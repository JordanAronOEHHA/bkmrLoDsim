start_time <- Sys.time()

library(bkmr)
# library(ggplot2)
library(dplyr)
library(mice)
library(qgcomp)

# 1. Setup Simulation Parameters
#TODO set array number to seed if present
seed <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (is.na(seed)){
  seed <- 999999
  mcmc_iter <- 100
} else {
  mcmc_iter <- 10000
}
set.seed(seed)

####### Control Parameters #######

# Defaults 
n <- 300
n_te <- 1000
lod_quantile <- 0.9
exposure_dist <- "lnorm" # Options: lnorm unif gamma
h_dist <- "nonlinear" # Options: linear nonlinear

# Command line arguments (override defaults if included)
args <- commandArgs(TRUE)
if (length(args) >= 1) n <- as.numeric(args[1])
if (length(args) >= 2) lod_quantile <- as.numeric(args[2])
if (length(args) >= 3) exposure_dist <- args[3]
if (length(args) >= 4) h_dist <- args[4]

#prints out current settings for reference when looking at results
print("Simulation Settings:")
print(paste("seed:", seed))
print(paste("n:", n))
print(paste("lod_quantile:", lod_quantile))
print(paste("exposure_dist:", exposure_dist))
print(paste("h_dist:", h_dist))

####### Hyper Parameters #######
p <- 3

#number of completed datasets
m_imputations <- 5
#number of iterations for the imputation model
mi_maxit <- 10
mi_seed <- seed

####### Functions ########

##### MSE functions #####
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

#empty table makes for when no complete case data is available 
empty_mse_by_lod_count <- function(p) {
  tibble::tibble(
    group = c(0:p, NA),
    n_obs = 0L,
    mse = NA_real_
  )
}

empty_mse_by_first2_lod <- function() {
  tibble::tibble(
    group = c("++", "+-", "-+", "--", "Overall"),
    n_obs = 0L,
    mse = NA_real_
  )
}

##### Data functions #####

CensorData <- function(Z_true, lod) {
  Z_obs <- Z_true
  for(j in 1:ncol(Z_true)) {
    Z_obs[Z_true[,j] < lod[j], j] <- NA
  }
  return(Z_obs)
}

SingleImputation <- function(Z_obs, lod) {
  Z_impute <- Z_obs
  for(j in 1:ncol(Z_obs)) {
    Z_impute[is.na(Z_impute[,j]), j] <- lod[j] / sqrt(2)
  }
  return(Z_impute)
}
AugmentData <- function(Z_obs, lod, center = NULL, scale_vals = NULL) {
  p <- ncol(Z_obs)
  Z_aug_cont <- log(Z_obs)
  Z_aug_ind <- matrix(1, nrow = nrow(Z_obs), ncol = p)

  for (j in seq_len(p)) {
    is_censored <- is.na(Z_obs[, j])
    Z_aug_cont[is_censored, j] <- log(lod[j])
    Z_aug_ind[is_censored, j] <- 0
  }

  Z_aug_cont <- sweep(Z_aug_cont, 2, log(lod), "-")

  if (is.null(center)) {
    center <- colMeans(Z_aug_cont)
  }

  if (is.null(scale_vals)) {
    scale_vals <- apply(Z_aug_cont, 2, sd)
  }

  Z_aug_cont_scaled <- scale(
    Z_aug_cont,
    center = center,
    scale = scale_vals
  )

  list(
    Z = cbind(Z_aug_cont_scaled, Z_aug_ind),
    center = center,
    scale = scale_vals
  )
}

####### Simulation #######
#_tr for training, _te for testing


#### Exposure and LoD
#LoD is based on quantile of true distribution
#Alternative is to censor x% of data but not current use 
#lod <- apply(Z_true, 2, quantile, probs = lod_quantile)
if (exposure_dist == "lnorm") {
  norm1 <- 2
  norm2 <- 1
  Z_true_tr <- matrix(exp(rnorm(n * p,mean = norm1,sd = norm2)), ncol = p)
  Z_true_te <- matrix(exp(rnorm(n_te * p,mean = norm1,sd = norm2)), ncol = p)
  lod <- exp(qnorm(lod_quantile, mean = norm1, sd = norm2))
  lod <- rep(lod, p)

} else if (exposure_dist == "unif") {
  unif1 = 0
  unif2 = 4
  Z_true_tr <- matrix(exp(runif(n * p, min = 0, max = 4)), ncol = p)
  Z_true_te <- matrix(exp(runif(n_te * p, min = 0, max = 4)), ncol = p)
  lod <- exp(qunif(lod_quantile, min = unif1, max = unif2))
  lod <- rep(lod, p)

} else if (exposure_dist == "gamma") {
  gamma1 = 2
  gamma2 = 1
  Z_true_tr <- matrix(exp(rgamma(n * p, shape = 2, scale = 1)), ncol = p)
  Z_true_te <- matrix(exp(rgamma(n_te * p, shape = 2, scale = 1)), ncol = p)
  lod <- exp(qgamma(lod_quantile, shape = gamma1, scale = gamma2))
  lod <- rep(lod, p)
  
} else {
  stop("exposure_dist must be 'lnorm' or 'unif' or 'gamma'")
}

# Create observed data with censoring
Z_obs_tr <- CensorData(Z_true_tr, lod)
Z_obs_te <- CensorData(Z_true_te, lod)

#get indices for complete case, can be empty with higher % censoring
complete_case_idx_tr <- complete.cases(Z_obs_tr)
complete_case_idx_te <- complete.cases(Z_obs_te)

#### Response
Z_resp_tr <- log(Z_true_tr)
Z_resp_te <- log(Z_true_te)
plogis_mean <- 2

if (h_dist == "nonlinear") {
  h_true_tr <- 4 * plogis(1/4 * (Z_resp_tr[,1] + Z_resp_tr[,2] + 1/2 * (Z_resp_tr[,1]) * (Z_resp_tr[,2])), location = plogis_mean, scale = 0.5)
  h_true_te <- 4 * plogis(1/4 * (Z_resp_te[,1] + Z_resp_te[,2] + 1/2 * (Z_resp_te[,1]) * (Z_resp_te[,2])), location = plogis_mean, scale = 0.5)
} else if (h_dist == "linear") {
  h_true_tr <- Z_resp_tr[, 1] + Z_resp_tr[, 2] + 0.5 * Z_resp_tr[, 1] * Z_resp_tr[, 2]
  h_true_te <- Z_resp_te[, 1] + Z_resp_te[, 2] + 0.5 * Z_resp_te[, 1] * Z_resp_te[, 2]
} else {
  stop("h_dist must be 'linear' or 'nonlinear'")
}

y_tr <- h_true_tr + rnorm(n, sd = 1)
y_te <- h_true_te + rnorm(n_te, sd = 1)

y_complete_case_tr <- y_tr[complete_case_idx_tr]

##### Models #####
# A. Uncensored
# B. Imputation (LoD / sqrt(2))
# C. Augmented (Indicator + Continuous)
# D. Complete case
# E. Truncated MI (tobit lognormal?)
# TODO
# F. Pseudo-Gibbs (Carli et al) ie. impute using gibbs in each MCMC iteration 
#this might not make sense for computational/time to implement


##### A. Uncensored #####
Z_uncensored_tr <- Z_true_tr
Z_uncensored_te <- Z_true_te

##### B. Imputation (LoD / sqrt(2)) #####
Z_impute_tr <- SingleImputation(Z_obs_tr, lod)
Z_impute_te <- SingleImputation(Z_obs_te, lod)

##### C. Augmented (Indicator + Continuous) #####
aug_train <- AugmentData(Z_obs = Z_obs_tr,lod = lod)
Z_aug_tr <- aug_train$Z

aug_test <- AugmentData(Z_obs = Z_obs_te, lod = lod, center = aug_train$center, scale_vals = aug_train$scale)
Z_aug_te <- aug_test$Z

##### D. Complete case #####
Z_complete_case_tr <- Z_uncensored_tr[complete_case_idx_tr, , drop = FALSE]
Z_complete_case_te <- Z_uncensored_te[complete_case_idx_te, , drop = FALSE]

##### E. Truncated multiple imputation using qgcomp::mice.impute.leftcenslognorm #####
mice.impute.leftcenslognorm <- qgcomp::mice.impute.leftcenslognorm
mi_data_tr <- cbind(y_tr, Z_obs_tr) |> as.data.frame()
mi_data_te <- cbind(y_te, Z_obs_te) |> as.data.frame()
colnames(mi_data_tr) <- c("y", paste0("z", seq_len(p)))
colnames(mi_data_te) <- c("y", paste0("z", seq_len(p)))

#call with no iterations to get default settings
mi_init_tr <- mice(mi_data_tr, maxit = 0, printFlag = FALSE)

method_trunc <- mi_init_tr$method
method_trunc[2:(p+1)] <- "leftcenslognorm"

#only uses y to impute missing values, no covariates in current form and all z have some censoring
predictor_matrix_trunc <- matrix(0,nrow = p+1, ncol = p+1)
predictor_matrix_trunc[2:(p+1),1] <- 1
colnames(predictor_matrix_trunc) <- c("y", paste0("z", seq_len(p)))
rownames(predictor_matrix_trunc) <- c("y", paste0("z", seq_len(p)))

mids_trunc_tr <- mice(
  data = mi_data_tr,
  m = m_imputations,
  maxit = mi_maxit,
  method = method_trunc,
  predictorMatrix = predictor_matrix_trunc,
  lod = c(NA, lod),
  seed = mi_seed,
  printFlag = FALSE
)

mids_trunc_te <- mice(
  data = mi_data_te,
  m = m_imputations,
  maxit = mi_maxit,
  method = method_trunc,
  predictorMatrix = predictor_matrix_trunc,
  lod = c(NA, lod),
  seed = mi_seed,
  printFlag = FALSE
)

Z_trunc_mi_raw_list_tr <- complete(mids_trunc_tr, action = "all") |>
  lapply(function(dat) {
    dat[, paste0("z", seq_len(p)), drop = FALSE] |>
      as.matrix()
  })

Z_trunc_mi_raw_list_te <- complete(mids_trunc_te, action = "all") |>
  lapply(function(dat) {
    dat[, paste0("z", seq_len(p)), drop = FALSE] |>
      as.matrix()
  })

log_obs_for_scaling_tr <- log(Z_obs_tr)
trunc_mi_center_tr <- colMeans(log_obs_for_scaling_tr, na.rm = TRUE)
trunc_mi_scale_tr <- apply(log_obs_for_scaling_tr, 2, sd, na.rm = TRUE)

##### Scaling #####
#augmented continuous scalet seperately due to the LoD adjustment
#multiple imputation tobit also needs to be handled separately, done in prevous block
#training scale is applied to testing scale to ensure that scaling used in training is applied to testing data
#common scaling for all MI data sets is used to ensure that they are on the same scale for pooling results 

### uncensored scaling
uncens_center <- colMeans(log(Z_true_tr))
uncens_scale <- apply(log(Z_true_tr), 2, sd)

Z_uncensored_tr <- scale( log(Z_true_tr), center = uncens_center, scale = uncens_scale)
Z_uncensored_te <- scale( log(Z_true_te), center = uncens_center, scale = uncens_scale)

### SI scaling

impute_center <- colMeans(log(Z_impute_tr))
impute_scale <- apply(log(Z_impute_tr), 2, sd)

Z_impute_tr <- scale(log(Z_impute_tr), center = impute_center, scale = impute_scale)
Z_impute_te <- scale(log(Z_impute_te), center = impute_center, scale = impute_scale)


### CC scaling

if (nrow(Z_complete_case_tr) >= 2) {
  complete_center <- colMeans(log(Z_complete_case_tr))
  complete_scale <- apply(log(Z_complete_case_tr), 2, sd)

  Z_complete_case_tr <- scale(
    log(Z_complete_case_tr),
    center = complete_center,
    scale = complete_scale
  )

  if (nrow(Z_complete_case_te) > 0) {
    Z_complete_case_te <- scale(
      log(Z_complete_case_te),
      center = complete_center,
      scale = complete_scale
    )
  }
}


##### BKMR ##### 

#create empty lists to store results for MI 
Z_trunc_mi_list_tr <- vector("list", length = m_imputations)
fit_trunc_mi_list_tr <- vector("list", length = m_imputations)
pred_trunc_mi_list_tr <- vector("list", length = m_imputations)

Z_trunc_mi_list_te <- vector("list", length = m_imputations)
pred_trunc_mi_list_te <- vector("list", length = m_imputations)

#used to store results by how many values are below the LoD
#ie group 0 is all values above the LoD, which is our focus
group_tr <- rowSums(is.na(Z_obs_tr))
group_te <- rowSums(is.na(Z_obs_te))

group_complete_case_tr <- group_tr[complete_case_idx_tr]
group_complete_case_te <- group_te[complete_case_idx_te]

###############################################################################
#Model calls and MSE calculation


##### BKMR Uncensored #####
m_uncensored <- kmbayes(y = y_tr, Z = Z_uncensored_tr, iter = mcmc_iter)
pred_uncensored_tr <- SamplePred(m_uncensored, Znew = Z_uncensored_tr)
results_uncens_tr  <- mse_by_lod_count(h_true_tr, pred_uncensored_tr, group_tr, p)
results_uncens_first2_tr  <- mse_by_first2_lod(h_true_tr, pred_uncensored_tr, Z_obs_tr)

#remove cbind(0) when adding covariates
#need this when predicting on new data and dont have any new fixed covariates
pred_uncensored_te <- SamplePred(m_uncensored, Znew = Z_uncensored_te,Xnew = cbind(0))
results_uncens_te  <- mse_by_lod_count(h_true_te, pred_uncensored_te, group_te, p)
results_uncens_first2_te  <- mse_by_first2_lod(h_true_te, pred_uncensored_te, Z_obs_te)

##### BKMR Single Imputation #####

m_impute <- kmbayes(y = y_tr, Z = Z_impute_tr, iter = mcmc_iter)
pred_impute_tr <- SamplePred(m_impute, Znew = Z_impute_tr)
results_impute_tr <- mse_by_lod_count(h_true_tr, pred_impute_tr, group_tr, p)
results_imputes_first2_tr <- mse_by_first2_lod(h_true_tr, pred_impute_tr, Z_obs_tr)

pred_impute_te <- SamplePred(m_impute, Znew = Z_impute_te, Xnew = cbind(0))
results_impute_te <- mse_by_lod_count(h_true_te, pred_impute_te, group_te, p)
results_imputes_first2_te <- mse_by_first2_lod(h_true_te, pred_impute_te, Z_obs_te)


##### BKMR Missing Indicator Method #####

m_augmented <- kmbayes(y = y_tr, Z = Z_aug_tr, iter = mcmc_iter)
pred_augmented_tr <- SamplePred(m_augmented, Znew = Z_aug_tr, Xnew = cbind(0))
results_augment_tr <- mse_by_lod_count(h_true_tr, pred_augmented_tr, group_tr, p)
results_augments_first2_tr <- mse_by_first2_lod(h_true_tr, pred_augmented_tr, Z_obs_tr)

pred_augmented_te <- SamplePred(m_augmented, Znew = Z_aug_te, Xnew = cbind(0))
results_augment_te <- mse_by_lod_count(h_true_te, pred_augmented_te, group_te, p)
results_augments_first2_te <- mse_by_first2_lod(h_true_te, pred_augmented_te, Z_obs_te)


##### BKMR Complete Case #####

n_complete_tr <- nrow(Z_complete_case_tr)

if (n_complete_tr >= 2) {
  m_complete_case <- kmbayes( y = y_complete_case_tr, Z = Z_complete_case_tr, iter = mcmc_iter )
  pred_complete_case_tr <- SamplePred( m_complete_case, Znew = Z_complete_case_tr)
  results_complete_case_tr <- mse_by_lod_count(h_true_tr[complete_case_idx_tr],pred_complete_case_tr,group_complete_case_tr,p)
  results_complete_cases_first2_tr <- mse_by_first2_lod( h_true_tr[complete_case_idx_tr], pred_complete_case_tr, Z_obs_tr[complete_case_idx_tr, , drop = FALSE])

  if (nrow(Z_complete_case_te) > 0 && !is.null(m_complete_case)) {
    pred_complete_case_te <- SamplePred( m_complete_case, Znew = Z_complete_case_te, Xnew = cbind(0))
    results_complete_case_te <- mse_by_lod_count( h_true_te[complete_case_idx_te], pred_complete_case_te, group_complete_case_te, p)
    results_complete_cases_first2_te <- mse_by_first2_lod( h_true_te[complete_case_idx_te], pred_complete_case_te, Z_obs_te[complete_case_idx_te, , drop = FALSE])

  } else {
    pred_complete_case_te <- NULL
    results_complete_case_te <- empty_mse_by_lod_count(p)
    results_complete_cases_first2_te <- empty_mse_by_first2_lod()
  }

} else {
  m_complete_case <- NULL
  pred_complete_case_tr <- NULL
  pred_complete_case_te <- NULL
  results_complete_case_tr <- empty_mse_by_lod_count(p)
  results_complete_cases_first2_tr <- empty_mse_by_first2_lod()
  results_complete_case_te <- empty_mse_by_lod_count(p)
  results_complete_cases_first2_te <- empty_mse_by_first2_lod()
}


#MI model
for (m in seq_len(m_imputations)) {
  Z_trunc_mi_list_tr[[m]] <- scale( log(Z_trunc_mi_raw_list_tr[[m]]), center = trunc_mi_center_tr, scale = trunc_mi_scale_tr)
  Z_trunc_mi_list_te[[m]] <- scale( log(Z_trunc_mi_raw_list_te[[m]]), center = trunc_mi_center_tr, scale = trunc_mi_scale_tr)
  fit_trunc_mi_list_tr[[m]] <- kmbayes(y = y_tr,Z = Z_trunc_mi_list_tr[[m]],iter = mcmc_iter)
  pred_trunc_mi_list_tr[[m]] <- SamplePred(fit_trunc_mi_list_tr[[m]],Znew = Z_trunc_mi_list_tr[[m]])
  pred_trunc_mi_list_te[[m]] <- SamplePred(fit_trunc_mi_list_tr[[m]],Znew = Z_trunc_mi_list_te[[m]],Xnew = cbind(0))
}

pred_trunc_mi_tr <- do.call(rbind, pred_trunc_mi_list_tr)
pred_trunc_mi_te <- do.call(rbind, pred_trunc_mi_list_te)

results_trunc_mi_tr <- mse_by_lod_count( h_true_tr, pred_trunc_mi_tr, group_tr, p)
results_trunc_mi_first2_tr <- mse_by_first2_lod( h_true_tr, pred_trunc_mi_tr, Z_obs_tr)
results_trunc_mi_te <- mse_by_lod_count( h_true_te, pred_trunc_mi_te, group_te, p)
results_trunc_mi_first2_te <- mse_by_first2_lod( h_true_te, pred_trunc_mi_te, Z_obs_te)

########################### Results Compilation and Saving ###########################

sim_results <- list(
  settings = list(
    seed = seed,
    n = n,
    n_te = n_te,
    p = p,
    lod_quantile = lod_quantile,
    exposure_dist = exposure_dist,
    h_dist = h_dist,
    mcmc_iter = mcmc_iter,
    m_imputations = m_imputations,
    mi_maxit = mi_maxit
  ),
  logistics = list(
    run_time = Sys.time() - start_time,
    run_mem = sum(gc()[, 6])
  ),
  results = list(
    training = list(
      mse_by_lod_count = list(
        uncensored = results_uncens_tr,
        impute = results_impute_tr,
        augment = results_augment_tr,
        complete_case = results_complete_case_tr,
        trunc_mi = results_trunc_mi_tr
      ),
      mse_by_first2_lod = list(
        uncensored = results_uncens_first2_tr,
        impute = results_imputes_first2_tr,
        augment = results_augments_first2_tr,
        complete_case = results_complete_cases_first2_tr,
        trunc_mi = results_trunc_mi_first2_tr
      )
    ),
    testing = list(
      mse_by_lod_count = list(
        uncensored = results_uncens_te,
        impute = results_impute_te,
        augment = results_augment_te,
        complete_case = results_complete_case_te,
        trunc_mi = results_trunc_mi_te
      ),
      mse_by_first2_lod = list(
        uncensored = results_uncens_first2_te,
        impute = results_imputes_first2_te,
        augment = results_augments_first2_te,
        complete_case = results_complete_cases_first2_te,
        trunc_mi = results_trunc_mi_first2_te
      )
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
