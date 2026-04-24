start_time <- Sys.time()

library(bkmr)
# library(ggplot2)
library(dplyr)

# 1. Setup Simulation Parameters
#TODO set array number to seed if present
seed <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (is.na(seed)){
  seed <- 9999999
  mcmc_iter <- 100
} else {
  mcmc_iter <- 10000
}
set.seed(seed)

####### Control Parameters #######

# Defaults 
n <- 300
lod_quantile <- 0.75
exposure_dist <- "lnorm" # Options: lnorm unif gamma
h_dist <- "nonlinear" # Options: linear nonlinear

# Command line arguments (override defaults if included)
args <- commandArgs(TRUE)
if (length(args) >= 1) n <- as.numeric(args[1])
if (length(args) >= 2) lod_quantile <- as.numeric(args[2])
if (length(args) >= 3) exposure_dist <- args[3]
if (length(args) >= 4) h_dist <- args[4]

print("Simulation Settings:")
print(paste("seed:", seed))
print(paste("n:", n))
print(paste("lod_quantile:", lod_quantile))
print(paste("exposure_dist:", exposure_dist))
print(paste("h_dist:", h_dist))

####### Hyper Parameters #######
p <- 2
rho_corr <- 0 # Correlation between chemicals



####### Functions ########
mse <- function(true, pred) mean((true - apply(pred, 2, mean))^2)

mse_by_lod_count <- function(h_true, pred, group, p_val = NULL) {
  if (is.null(p_val)) p_val <- max(group)

  results <- lapply(0:p_val, function(k) {
    idx <- which(group == k)

    if (length(idx) > 0) {
      data.frame(
        group = k,
        n_obs = length(idx),
        mse = mse(h_true[idx], pred[, idx, drop = FALSE])
      )
    } else {
      data.frame(
        group = k,
        n_obs = 0L,
        mse = NA_real_
      )
    }
  }) |>
    dplyr::bind_rows()

  total_row <- data.frame(
    group = NA,
    n_obs = length(h_true),
    mse = mse(h_true, pred)
  )

  dplyr::bind_rows(results, total_row)
}


mse_by_first2_lod <- function(h_true, pred, Z_obs) {
  z1_below <- is.na(Z_obs[, 1])
  z2_below <- is.na(Z_obs[, 2])

  group <- dplyr::case_when(
    !z1_below & !z2_below ~ "++",
    z1_below & !z2_below ~ "-+",
    !z1_below & z2_below ~ "+-",
    z1_below & z2_below ~ "--"
  )

  pred_mean <- colMeans(pred)

  observed <- tibble::tibble(
    group = group,
    h_true = h_true,
    pred_mean = pred_mean
  ) |>
    dplyr::summarise(
      n_obs = dplyr::n(),
      mse = mean((h_true - pred_mean)^2),
      .by = group
    )

  template <- tibble::tibble(group = c("++", "+-", "-+", "--"))

  by_group <- template |>
    dplyr::left_join(observed, by = "group") |>
    dplyr::mutate(
      n_obs = dplyr::coalesce(n_obs, 0L),
      mse = ifelse(n_obs == 0L, NA_real_, mse)
    )

  overall <- tibble::tibble(
    group = "Overall",
    n_obs = length(h_true),
    mse = mean((h_true - pred_mean)^2)
  )

  dplyr::bind_rows(by_group, overall)
}
####### Simulation #######

#### Exposure
if (exposure_dist == "lnorm") {
  Z_true <- matrix(exp(rnorm(n * p,2,1)), ncol = p)
} else if (exposure_dist == "unif") {
  Z_true <- matrix(exp(runif(n * p, min = 0, max = 4)), ncol = p)
} else if (exposure_dist == "gamma") {
  Z_true <- matrix(exp(rgamma(n * p, shape = 2, scale = 1)), ncol = p)
} else {
  stop("exposure_dist must be 'lnorm' or 'unif' or 'gamma'")
}

#### LoD
lod <- apply(Z_true, 2, quantile, probs = lod_quantile)

# Create observed data
Z_obs <- Z_true
for(j in 1:p) {
  Z_obs[Z_true[,j] < lod[j], j] <- NA
}
complete_case_idx <- complete.cases(Z_obs)
if (!any(complete_case_idx)) {
  stop("No complete cases available for complete-case model fitting.")
}

##### Models
# A. Uncensored
# B. Imputation (LoD / sqrt(2))
# C. Augmented (Indicator + Continuous)
# D. Complete case
# TODO
# E. Truncated MI (tobit lognormal)
# F. Pseudo-Gibbs (Carli et al) ie. impute using gibbs in each MCMC iteration 
#this might not make sense for computational/time to implement


# A. Uncensored 
Z_uncensored <- Z_true

# B. Imputation (LoD / sqrt(2))
Z_impute <- Z_obs
for(j in 1:p) {
  Z_impute[is.na(Z_impute[,j]), j] <- lod[j] / sqrt(2)
}


# C. Augmented (Indicator + Continuous)
#Current functional form assumes everything is log-dist (normal, gamma, unif) so need this logic for all. 
#put this block in if statement and uncomment following block if this changes
  Z_aug_cont <- log(Z_obs)
  Z_aug_ind  <- matrix(1, nrow = n, ncol = p)
  for(j in 1:p) {
    Z_aug_cont[is.na(Z_aug_cont[,j]), j] <- log(lod[j])
    # Z_aug_cont[,j] <- Z_aug_cont[,j] - lod[j]
    Z_aug_ind[is.na(Z_obs[,j]), j] <- 0
  }
  Z_aug_cont <- sweep(Z_aug_cont, 2, log(lod), "-")
  Z_aug <- cbind(scale(Z_aug_cont), Z_aug_ind)

# if (exposure_dist == "lnorm") {}else {
#   Z_aug_cont <- Z_obs
#   Z_aug_ind  <- matrix(1, nrow = n, ncol = p)
#   for(j in 1:p) {
#     Z_aug_cont[is.na(Z_aug_cont[,j]), j] <- lod[j]
#     Z_aug_cont[,j] <- Z_aug_cont[,j] - lod[j]
#     Z_aug_ind[is.na(Z_obs[,j]), j] <- 0
#   }
#   Z_aug <- cbind(scale(Z_aug_cont), Z_aug_ind)
# }

#D. Complete case
Z_complete_case <- Z_uncensored[complete_case_idx, , drop = FALSE]


#### Scaling
#except for augmented continuous since that needs to be handled seperately due to the LoD adjustment
#should be in the code block directly above this 
if (exposure_dist == "lnorm") {
  Z_uncensored <- scale(log(Z_true))
  Z_impute <- scale(log(Z_impute))
  Z_complete_case <- scale(log(Z_complete_case))
} else {
  Z_uncensored <- scale(log(Z_true))
  Z_impute <- scale(log(Z_impute))
  Z_complete_case <- scale(log(Z_complete_case))
}



#### Response
if (exposure_dist == "lnorm") {
  Z_response <- log(Z_true)
  plogis_mean <- 2
} else {
  Z_response <- log(Z_true)
  plogis_mean <- 2
}

if (h_dist == "nonlinear") {
  h_true <- 4 * plogis(1/4 * (Z_response[,1] + Z_response[,2] + 1/2 * (Z_response[,1]) * (Z_response[,2])), location = plogis_mean, scale = 0.5)
} else if (h_dist == "linear") {
  h_true <- (Z_response[,1] + Z_response[,2] - 2*plogis_mean)^2
} else {
  stop("h_dist must be 'linear' or 'nonlinear'")
}

y <- h_true + rnorm(n, sd = 1)
y_complete_case <- y[complete_case_idx]

# Run Models

m_uncensored <- kmbayes(y = y, Z = Z_uncensored, iter = mcmc_iter)
m_impute     <- kmbayes(y = y, Z = Z_impute, iter = mcmc_iter)
m_augmented  <- kmbayes(y = y, Z = Z_aug, iter = mcmc_iter)
m_complete_case <- kmbayes(y = y_complete_case, Z = Z_complete_case, iter = mcmc_iter)

# Predict h for each model on the fitted data
pred_uncensored <- SamplePred(m_uncensored, Znew = Z_uncensored)
pred_impute     <- SamplePred(m_impute, Znew = Z_impute)
pred_augmented  <- SamplePred(m_augmented, Znew = Z_aug)
pred_complete_case <- SamplePred(m_complete_case, Znew = Z_complete_case)

# Determine how many values per observation are below LoD
group <- rowSums(is.na(Z_obs))
group_complete_case <- group[complete_case_idx]
results_uncens  <- mse_by_lod_count(h_true, pred_uncensored, group, p)
results_impute  <- mse_by_lod_count(h_true, pred_impute, group, p)
results_augment <- mse_by_lod_count(h_true, pred_augmented, group, p)
results_complete_case <- mse_by_lod_count(h_true[complete_case_idx], pred_complete_case, group_complete_case, p)

results_uncens_first2  <- mse_by_first2_lod(h_true, pred_uncensored, Z_obs)
results_imputes_first2  <- mse_by_first2_lod(h_true, pred_impute, Z_obs)
results_augments_first2 <- mse_by_first2_lod(h_true, pred_augmented, Z_obs)
results_complete_case_first2 <- mse_by_first2_lod(h_true[complete_case_idx], pred_complete_case, Z_obs[complete_case_idx, , drop = FALSE])

sim_results <- list(
  settings = list(
    seed = seed,
    n = n,
    p = p,
    lod_quantile = lod_quantile,
    exposure_dist = exposure_dist,
    h_dist = h_dist,
    mcmc_iter = mcmc_iter
  ),
  logistics = list(
    run_time = Sys.time() - start_time,
    run_mem = sum(gc()[,6])
  ),
  results = list(
    mse_by_lod_count = list(
      uncensored = results_uncens,
      impute = results_impute,
      augment = results_augment,
      complete_case = results_complete_case
    ),
    mse_by_first2_lod = list(
      uncensored = results_uncens_first2,
      impute = results_imputes_first2,
      augment = results_augments_first2,
      complete_case = results_complete_case_first2
    )
  )
)

name <- paste0(
  "SimResults_seed", seed,
  "_n", n,
  "_lod", lod_quantile,
  "_", exposure_dist,
  "_", h_dist,
  ".rds"
)

saveRDS(sim_results, file = name)
