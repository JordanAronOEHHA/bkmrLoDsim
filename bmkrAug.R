library(bkmr)
library(ggplot2)
library(dplyr)

# 1. Setup Simulation Parameters
seed <- 9
set.seed(seed)

####### Control Parameters #######

#TODO: take these as command line arguments
n <- 300
lod_quantile <- 0.75
exposure_dist <- "lnorm" # Options: lnorm unif 
h_dist <- "nonlinear" # Options: linear nonlinear

####### Hyper Parameters #######
p <- 3
rho_corr <- 0 # Correlation between chemicalse
offset <- 0
mcmc_iter <- 500


####### Functions ########
mse <- function(true, pred) mean((true - apply(pred, 2, mean))^2)

mse_by_lod_count <- function(h_true, pred, n_below_lod, p_val = NULL) {
  if(is.null(p_val)) p_val <- max(n_below_lod)
  results <- data.frame(n_below_lod = integer(), n_obs = integer(), mse = numeric())
  for(k in 0:p_val) {
    idx <- which(n_below_lod == k)
    if(length(idx) > 0) {
      mse_val <- mse(h_true[idx], pred[, idx, drop = FALSE])
      results <- rbind(results, data.frame(n_below_lod = k, n_obs = length(idx), mse = mse_val))
    }
  }
  total_mse <- mse(h_true, pred)
  results <- rbind(results, data.frame(n_below_lod = NA, n_obs = length(h_true), mse = total_mse))
  return(results)
}


mse_by_first2_lod <- function(h_true, pred, Z_obs) {
  z1_below <- is.na(Z_obs[, 1])
  z2_below <- is.na(Z_obs[, 2])

  group <- dplyr::case_when(
    !z1_below & !z2_below ~ "Z1 above, Z2 above",
    z1_below & !z2_below ~ "Z1 below, Z2 above",
    !z1_below & z2_below ~ "Z1 above, Z2 below",
    z1_below & z2_below ~ "Z1 below, Z2 below"
  )

  pred_mean <- colMeans(pred)

  dplyr::bind_rows(
    tibble::tibble(
      group = group,
      h_true = h_true,
      pred_mean = pred_mean
    ) |>
      dplyr::summarise(
        n_obs = dplyr::n(),
        mse = mean((h_true - pred_mean)^2),
        .by = group
      ),
    tibble::tibble(
      group = "Overall",
      n_obs = length(h_true),
      mse = mean((h_true - pred_mean)^2)
    )
  )
}

####### Simulation #######

#### Exposure
if (exposure_dist == "lnorm") {
  Z_true <- matrix(rlnorm(n * p,log(2) - .5,1), ncol = p)
} else if (exposure_dist == "unif") {
  Z_true <- matrix(runif(n * p, min = 0, max = 4), ncol = p) + offset
} else {
  stop("exposure_dist must be 'lnorm' or 'unif'")
}

#### Response
if (h_dist == "nonlinear") {
  h_true <- 4 * plogis(1/4 * (Z_true[,1] + Z_true[,2] + 1/2 * (Z_true[,1]) * (Z_true[,2])), location = 2, scale = 0.3)
} else if (h_dist == "linear") {
  h_true <- (Z_true[,1] + Z_true[,2])^2
}

y <- h_true + rnorm(n, sd = 1) 

#### LoD
lod <- apply(Z_true, 2, quantile, probs = lod_quantile)

# Create observed data
Z_obs <- Z_true
for(j in 1:p) {
  Z_obs[Z_true[,j] < lod[j], j] <- NA
}

# A. Uncensored (Oracle/Gold Standard)
Z_uncensored <- scale(Z_true)

# B. Imputation (LoD / sqrt(2))
Z_impute <- Z_obs
for(j in 1:p) {
  Z_impute[is.na(Z_impute[,j]), j] <- lod[j] / sqrt(2)
}
Z_impute <- scale(Z_impute)

# C. Augmented (Indicator + Continuous)
Z_aug_cont <- Z_obs
Z_aug_ind  <- matrix(1, nrow = n, ncol = p)
for(j in 1:p) {
  Z_aug_cont[is.na(Z_aug_cont[,j]), j] <- lod[j]
  Z_aug_cont[,j] <- Z_aug_cont[,j] - lod[j]
  Z_aug_ind[is.na(Z_obs[,j]), j] <- 0
}
Z_aug <- cbind(scale(Z_aug_cont), Z_aug_ind)


# Run Models

m_uncensored <- kmbayes(y = y, Z = Z_uncensored, iter = mcmc_iter)
m_impute     <- kmbayes(y = y, Z = Z_impute, iter = mcmc_iter)
m_augmented  <- kmbayes(y = y, Z = Z_aug, iter = mcmc_iter)

# Predict h for each model on the fitted data
pred_uncensored <- SamplePred(m_uncensored, Znew = Z_uncensored)
pred_impute     <- SamplePred(m_impute, Znew = Z_impute)
pred_augmented  <- SamplePred(m_augmented, Znew = Z_aug)

# Determine how many values per observation are below LoD
n_below_lod <- rowSums(is.na(Z_obs))
results_uncens  <- mse_by_lod_count(h_true, pred_uncensored, n_below_lod, p)
results_impute  <- mse_by_lod_count(h_true, pred_impute, n_below_lod, p)
results_augment <- mse_by_lod_count(h_true, pred_augmented, n_below_lod, p)

results_uncens_first2  <- mse_by_first2_lod(h_true, pred_uncensored, Z_obs)
results_imputes_first2  <- mse_by_first2_lod(h_true, pred_impute, Z_obs)
results_augments_first2 <- mse_by_first2_lod(h_true, pred_augmented, Z_obs)

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
  results_uncens = results_uncens,
  results_impute = results_impute,
  results_augment = results_augment,
  results_uncens_first2 = results_uncens_first2,
  results_imputes_first2 = results_imputes_first2,
  results_augments_first2 = results_augments_first2
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
