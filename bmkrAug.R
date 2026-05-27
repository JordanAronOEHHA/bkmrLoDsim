q_fixed_effects <- 2L
fixed_effect_beta <- c(-.25, 0.25)

zero_fixed_effects <- function(n_rows) {
  matrix(0, nrow = n_rows, ncol = q_fixed_effects)
}


start_time <- Sys.time()

library(bkmr)
# library(ggplot2)
library(dplyr)
library(mice)
library(qgcomp)

#### Setup Parameters ####

# 1. Setup Simulation Parameters

RNGkind("L'Ecuyer-CMRG")

array_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))



if (is.na(array_id)){
  array_id <- 1
  seeds_per_job <- 1

  mcmc_iter <- 10
} else {
  seeds_per_job <- 4

  mcmc_iter <- 10000
}

seed_vec <- (array_id - 1) * seeds_per_job + seq_len(seeds_per_job)

##### Control Parameters #####
# Defaults 
n <- 400
lod_quantile <- 0.25
exposure_dist <- "norm" # Options: norm unif gamma
mean_offset <- 0
h_func <- 2
scale <- 0.5

# Command line arguments (override defaults if included)
args <- commandArgs(TRUE)
if (length(args) >= 1) n <- as.numeric(args[1])
if (length(args) >= 2) lod_quantile <- as.numeric(args[2])
if (length(args) >= 3) exposure_dist <- args[3]
if (length(args) >= 4) mean_offset <- as.numeric(args[4])
if (length(args) >= 5) h_func <- as.numeric(args[5])
if (length(args) >= 6) scale <- as.numeric(args[6])

#prints out current settings for reference when looking at results
print("Simulation Settings:")
print(paste("seed vector:", seed_vec))
print(paste("n:", n))
print(paste("lod_quantile:", lod_quantile))
print(paste("exposure_dist:", exposure_dist))
print(paste("mean_offset:", mean_offset))
print(paste("h_func:", h_func))
print(paste("scale:", scale))
print(paste("q_fixed_effects:", q_fixed_effects))
print(paste("fixed_effect_beta:", paste(fixed_effect_beta, collapse = ",")))

##### Hyper Parameters #####
p <- 4

#number of completed datasets
m_imputations <- 5
#number of iterations for the imputation model
mi_maxit <- 10


#### Functions #####

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

coverage_by_lod_count <- function(h_true, pred, group, p_val = NULL) {
  if (is.null(p_val)) p_val <- max(group)

  ci <- apply(
    pred,
    2,
    quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )

  covered <- h_true >= ci[1, ] & h_true <= ci[2, ]
  ci_width <- ci[2, ] - ci[1, ]

  observed <- tibble::tibble(
    group = group,
    covered = covered,
    ci_width = ci_width
  ) |>
    dplyr::summarise(
      n_obs = dplyr::n(),
      n_covered = sum(covered, na.rm = TRUE),
      empirical_coverage = mean(covered, na.rm = TRUE),
      mean_ci_width = mean(ci_width, na.rm = TRUE),
      .by = group
    )

  by_group <- tibble::tibble(group = 0:p_val) |>
    dplyr::left_join(observed, by = "group") |>
    dplyr::mutate(
      n_obs = dplyr::coalesce(n_obs, 0L),
      n_covered = dplyr::coalesce(n_covered, 0L),
      empirical_coverage = ifelse(n_obs == 0L, NA_real_, empirical_coverage),
      mean_ci_width = ifelse(n_obs == 0L, NA_real_, mean_ci_width)
    )

  overall <- tibble::tibble(
    group = NA,
    n_obs = length(h_true),
    n_covered = sum(covered, na.rm = TRUE),
    empirical_coverage = mean(covered, na.rm = TRUE),
    mean_ci_width = mean(ci_width, na.rm = TRUE)
  )

  dplyr::bind_rows(by_group, overall)
}

mse_by_first2_lod <- function(h_true, pred, Z_obs) {
  z1_below <- is.na(Z_obs[, 1])
  z2_below <- is.na(Z_obs[, 2])

  if (h_func == 2 | h_func == 3){
    group <- dplyr::case_when(
      !z1_below & !z2_below ~ "++",
      z1_below & !z2_below ~ "+-",
      !z1_below & z2_below ~ "+-",
      z1_below & z2_below ~ "--"
    )
  } else if (h_func == 1){
    group <- dplyr::case_when(
      !z1_below ~ "+",
      z1_below  ~ "-",
    )
  }
  

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

  if (h_func == 1){template <- tibble::tibble(group = c("+", "-"))}
  if (h_func == 2){template <- tibble::tibble(group = c("++", "+-", "--"))}
  if (h_func == 3){template <- tibble::tibble(group = c("++", "+-", "--"))}
  

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

coverage_by_first2_lod <- function(h_true, pred, Z_obs) {
  z1_below <- is.na(Z_obs[, 1])
  z2_below <- is.na(Z_obs[, 2])

  if (h_func == 2 | h_func == 3){
    group <- dplyr::case_when(
      !z1_below & !z2_below ~ "++",
      z1_below & !z2_below ~ "+-",
      !z1_below & z2_below ~ "+-",
      z1_below & z2_below ~ "--"
    )
  } else if (h_func == 1){
    group <- dplyr::case_when(
      !z1_below ~ "+",
      z1_below  ~ "-",
    )
  }

  ci <- apply(
    pred,
    2,
    quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )

  covered <- h_true >= ci[1, ] & h_true <= ci[2, ]
  ci_width <- ci[2, ] - ci[1, ]

  observed <- tibble::tibble(
    group = group,
    covered = covered,
    ci_width = ci_width
  ) |>
    dplyr::summarise(
      n_obs = dplyr::n(),
      n_covered = sum(covered, na.rm = TRUE),
      empirical_coverage = mean(covered, na.rm = TRUE),
      mean_ci_width = mean(ci_width, na.rm = TRUE),
      .by = group
    )

  if (h_func == 1){template <- tibble::tibble(group = c("+", "-"))}
  if (h_func == 2){template <- tibble::tibble(group = c("++", "+-", "--"))}
  if (h_func == 3){template <- tibble::tibble(group = c("++", "+-", "--"))}

  by_group <- template |>
    dplyr::left_join(observed, by = "group") |>
    dplyr::mutate(
      n_obs = dplyr::coalesce(n_obs, 0L),
      n_covered = dplyr::coalesce(n_covered, 0L),
      empirical_coverage = ifelse(n_obs == 0L, NA_real_, empirical_coverage),
      mean_ci_width = ifelse(n_obs == 0L, NA_real_, mean_ci_width)
    )

  overall <- tibble::tibble(
    group = "Overall",
    n_obs = length(h_true),
    n_covered = sum(covered, na.rm = TRUE),
    empirical_coverage = mean(covered, na.rm = TRUE),
    mean_ci_width = mean(ci_width, na.rm = TRUE)
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

empty_coverage_by_lod_count <- function(p) {
  tibble::tibble(
    group = c(0:p, NA),
    n_obs = 0L,
    n_covered = 0L,
    empirical_coverage = NA_real_,
    mean_ci_width = NA_real_
  )
}

empty_mse_by_first2_lod <- function() {
  if (h_func == 1){ grouping = c("+", "-", "Overall")}
  if (h_func == 2){ grouping = c("++", "+-", "--", "Overall")}
  if (h_func == 3){ grouping = c("++", "+-", "--", "Overall")}
  tibble::tibble(
    group = grouping,
    n_obs = 0L,
    mse = NA_real_
  )
}

empty_coverage_by_first2_lod <- function() {
  if (h_func == 1){ grouping = c("+", "-", "Overall")}
  if (h_func == 2){ grouping = c("++", "+-", "--", "Overall")}
  if (h_func == 3){ grouping = c("++", "+-", "--", "Overall")}
  tibble::tibble(
    group = grouping,
    n_obs = 0L,
    n_covered = 0L,
    empirical_coverage = NA_real_,
    mean_ci_width = NA_real_
  )
}

##### Data functions #####

true_h <- function(Z_raw, h_func, mean_offset, scale) {
  Z_log <- log(Z_raw)

  if (h_func == 1) {
    return(4 * plogis(Z_log[, 1], location = mean_offset, scale = scale))
  }

  if (h_func == 2) {
    return(
      2 * plogis(Z_log[, 1], location = mean_offset, scale = scale) +
        2 * plogis(Z_log[, 2], location = mean_offset, scale = scale)
    )
  }

  if (h_func == 3) {
    return(4 * plogis(0.5 * (Z_log[, 1] + Z_log[, 2]), location = mean_offset, scale = scale))
  }

  stop("Unsupported h_func")
}

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

  if (is.null(center)) {center <- colMeans(Z_aug_cont)}
  if (is.null(scale_vals)) {scale_vals <- apply(Z_aug_cont, 2, sd)}

  Z_aug_cont_scaled <- scale(Z_aug_cont,center = center,scale = scale_vals)

  list(
    Z = cbind(Z_aug_cont_scaled, Z_aug_ind),
    center = center,
    scale = scale_vals
  )
}


##### PIP functions #####
extract_augmented_chemical_pips <- function(fit, p, sel = NULL) {
  delta <- fit$delta

  if (is.null(delta)) {
    stop("fit$delta is not available")
  }

  if (is.null(sel)) {
    sel <- seq_len(nrow(delta)) > nrow(delta) / 2
  }

  delta <- delta[sel, , drop = FALSE]

  cont_idx <- seq_len(p)
  ind_idx <- p + seq_len(p)

  chem_pip_and <- vapply(seq_len(p), function(j) {
    mean(delta[, cont_idx[j]] == 1 & delta[, ind_idx[j]] == 1)
  }, numeric(1))

  chem_pip_or <- vapply(seq_len(p), function(j) {
    mean(delta[, cont_idx[j]] == 1 | delta[, ind_idx[j]] == 1)
  }, numeric(1))

  tibble::tibble(
    chemical = seq_len(p),
    pip_and = chem_pip_and,
    pip_or = chem_pip_or
  )
}

calc_sens_spec <- function(
  pip,
  thresholds = c(0.5, 0.75, 0.9)
) {

  if (h_func == 1){ true_active = c(TRUE, FALSE, FALSE, FALSE)}
  if (h_func == 2){ true_active = c(TRUE, TRUE, FALSE, FALSE)}
  if (h_func == 3){ true_active = c(TRUE, TRUE, FALSE, FALSE)}

  if (is.data.frame(pip)) {
    pip <- pip$PIP
  }

  tibble::tibble(threshold = thresholds) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      sensitivity = mean(pip[true_active] >= threshold),
      specificity = mean(pip[!true_active] < threshold)
    ) |>
    dplyr::ungroup()
}

##### Contrasts #####
get_fixed_values <- function(Z_true, lod, fixed, fixed_at) {
  fixed_above <- Z_true[Z_true[, fixed] >= lod[fixed], fixed]
  fixed_below <- Z_true[Z_true[, fixed] < lod[fixed], fixed]

  if (fixed_at == "below_lod") {
    return(sort(fixed_below))
  }

  if (fixed_at == "median_above_lod") {
    return(as.numeric(median(fixed_above, na.rm = TRUE)))
  }

  if (fixed_at == "lod") {
    return(lod[fixed])
  }

  if (fixed_at == "q25_above_lod") {
    return(as.numeric(quantile(fixed_above, 0.25, names = FALSE, type = 8)))
  }

  if (fixed_at == "q75_above_lod") {
    return(as.numeric(quantile(fixed_above, 0.75, names = FALSE, type = 8)))
  }

  stop("Unsupported fixed_at: ", fixed_at)
}

make_one_direction_contrast <- function(Z_true, lod, moving, fixed, fixed_at = "below_lod") {
  moving_above <- Z_true[Z_true[, moving] >= lod[moving], moving]
  fixed_values <- get_fixed_values(Z_true, lod, fixed, fixed_at)

  if (length(fixed_values) == 0) {
    stop("No fixed exposure values available for fixed_at = ", fixed_at)
  }

  z_low <- as.numeric(quantile(moving_above, 0.25, names = FALSE, type = 8))
  z_high <- as.numeric(quantile(moving_above, 0.75, names = FALSE, type = 8))

  baseline <- apply(Z_true, 2, median, na.rm = TRUE)

  Z_low <- matrix(rep(baseline, each = length(fixed_values)),
                  nrow = length(fixed_values))
  Z_high <- Z_low

  Z_low[, moving] <- z_low
  Z_high[, moving] <- z_high

  Z_low[, fixed] <- fixed_values
  Z_high[, fixed] <- fixed_values

  list(
    low = Z_low,
    high = Z_high,
    metadata = tibble::tibble(
      moving = moving,
      fixed = fixed,
      fixed_at = fixed_at,
      z_low = z_low,
      z_high = z_high,
      n_fixed_value = length(fixed_values),
      quantile_type = 8
    )
  )
}

make_symmetric_contrast <- function(Z_true, lod) {
  list(
    z1_given_z2_below = make_one_direction_contrast(Z_true, lod, moving = 1, fixed = 2, fixed_at = "below_lod"),
    z2_given_z1_below = make_one_direction_contrast(Z_true, lod, moving = 2, fixed = 1, fixed_at = "below_lod")
  )
}

prep_uncensored <- function(Z_raw) {
  scale(log(Z_raw), center = uncens_center, scale = uncens_scale)
}

prep_impute <- function(Z_raw) {
  Z_obs <- CensorData(Z_raw, lod)
  Z_imp <- SingleImputation(Z_obs, lod)
  scale(log(Z_imp), center = impute_center, scale = impute_scale)
}

prep_augmented <- function(Z_raw) {
  Z_obs <- CensorData(Z_raw, lod)
  AugmentData(
    Z_obs = Z_obs,
    lod = lod,
    center = aug_data$center,
    scale_vals = aug_data$scale
  )$Z
}

prep_complete_case <- function(Z_raw) {
  scale(log(Z_raw), center = complete_center, scale = complete_scale)
}

prep_trunc_mi <- function(Z_raw) {
  # Minimal version: represent censored contrast points by LoD/sqrt(2),
  # then use the common MI scaling already used in the script.
  Z_obs <- CensorData(Z_raw, lod)
  Z_imp <- SingleImputation(Z_obs, lod)
  scale(log(Z_imp), center = trunc_mi_center, scale = trunc_mi_scale)
}

estimate_contrast <- function(fit, contrast, prep_fun) {
  Z_low <- prep_fun(contrast$low)
  Z_high <- prep_fun(contrast$high)

  pred_low <- SamplePred(
    fit,
    Znew = Z_low,
    Xnew = zero_fixed_effects(nrow(Z_low))
  )

  pred_high <- SamplePred(
    fit,
    Znew = Z_high,
    Xnew = zero_fixed_effects(nrow(Z_high))
  )

  rowMeans(pred_high - pred_low)
}

estimate_contrast_mi <- function(fit_list, contrast, prep_fun) {
  draws <- lapply(fit_list, function(fit) {
    estimate_contrast(fit, contrast, prep_fun)
  })

  unlist(draws)
}

summarize_contrast <- function(draws, truth, method) {
  tibble::tibble(
    method = method,
    truth = truth,
    est = mean(draws),
    bias = mean(draws) - truth,
    ci_low = as.numeric(quantile(draws, 0.025)),
    ci_high = as.numeric(quantile(draws, 0.975)),
    covered = truth >= ci_low & truth <= ci_high
  )
}

beta_draw_matrix <- function(fit) {
  if (is.null(fit) || is.null(fit$beta)) {
    return(NULL)
  }

  beta <- fit$beta
  if (is.null(dim(beta))) {
    return(matrix(beta, ncol = q_fixed_effects))
  }

  beta <- as.matrix(beta)
  if (ncol(beta) == q_fixed_effects) {
    return(beta)
  }
  if (nrow(beta) == q_fixed_effects) {
    return(t(beta))
  }

  stop("Unexpected fixed-effect beta dimensions")
}

summarize_fixed_effect_draws <- function(beta_draws, method) {
  fixed_effect <- paste0("x", seq_len(q_fixed_effects))

  if (is.null(beta_draws) || nrow(beta_draws) == 0) {
    return(tibble::tibble(
      method = method,
      fixed_effect = fixed_effect,
      truth = fixed_effect_beta,
      est = NA_real_,
      bias = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      covered = NA
    ))
  }

  est <- colMeans(beta_draws)
  ci_low <- apply(beta_draws, 2, quantile, 0.025)
  ci_high <- apply(beta_draws, 2, quantile, 0.975)

  tibble::tibble(
    method = method,
    fixed_effect = fixed_effect,
    truth = fixed_effect_beta,
    est = est,
    bias = est - fixed_effect_beta,
    ci_low = as.numeric(ci_low),
    ci_high = as.numeric(ci_high),
    covered = fixed_effect_beta >= ci_low & fixed_effect_beta <= ci_high
  )
}

summarize_fixed_effects <- function(fit, method) {
  summarize_fixed_effect_draws(beta_draw_matrix(fit), method)
}

summarize_fixed_effects_mi <- function(fit_list, method = "trunc_mi") {
  draws <- lapply(fit_list, beta_draw_matrix)
  draws <- draws[!vapply(draws, is.null, logical(1))]

  summarize_fixed_effect_draws(
    if (length(draws) == 0) NULL else do.call(rbind, draws),
    method
  )
}

contrast_truth <- function(contrast) {
  mean(
    true_h(contrast$high, h_func, mean_offset, scale) -
      true_h(contrast$low, h_func, mean_offset, scale)
  )
}

estimate_all_methods_for_contrast <- function(contrast_obj, contrast_name) {
  truth <- contrast_truth(contrast_obj)
  metadata <- contrast_obj$metadata

  results <- list(
    summarize_contrast(
      estimate_contrast(m_uncensored, contrast_obj, prep_uncensored),
      truth,
      "uncensored"
    ),
    summarize_contrast(
      estimate_contrast(m_impute, contrast_obj, prep_impute),
      truth,
      "impute"
    ),
    summarize_contrast(
      estimate_contrast(m_augmented, contrast_obj, prep_augmented),
      truth,
      "augmented"
    ),
    summarize_contrast(
      estimate_contrast_mi(fit_trunc_mi_list, contrast_obj, prep_trunc_mi),
      truth,
      "trunc_mi"
    )
  )

  if (!is.null(m_complete_case)) {
    results <- c(
      results,
      list(
        summarize_contrast(
          estimate_contrast(m_complete_case, contrast_obj, prep_complete_case),
          truth,
          "complete_case"
        )
      )
    )
  }

  dplyr::bind_rows(results) |>
    dplyr::mutate(
      contrast = contrast_name,
      moving = metadata$moving,
      fixed = metadata$fixed,
      fixed_at = metadata$fixed_at,
      z_low = metadata$z_low,
      z_high = metadata$z_high,
      n_fixed_value = metadata$n_fixed_value,
      quantile_type = metadata$quantile_type,
      .before = method
    )
}

estimate_interaction_contrast <- function(fit, contrast_median, contrast_below, prep_fun) {
  estimate_contrast(fit, contrast_median, prep_fun) -
    estimate_contrast(fit, contrast_below, prep_fun)
}

estimate_interaction_contrast_mi <- function(fit_list, contrast_median, contrast_below, prep_fun) {
  draws <- lapply(fit_list, function(fit) {
    estimate_interaction_contrast(fit, contrast_median, contrast_below, prep_fun)
  })

  unlist(draws)
}

estimate_all_methods_for_interaction <- function(
  contrast_median,
  contrast_below,
  contrast_name
) {
  truth <- contrast_truth(contrast_median) - contrast_truth(contrast_below)
  metadata_median <- contrast_median$metadata
  metadata_below <- contrast_below$metadata

  results <- list(
    summarize_contrast(
      estimate_interaction_contrast(
        m_uncensored,
        contrast_median,
        contrast_below,
        prep_uncensored
      ),
      truth,
      "uncensored"
    ),
    summarize_contrast(
      estimate_interaction_contrast(
        m_impute,
        contrast_median,
        contrast_below,
        prep_impute
      ),
      truth,
      "impute"
    ),
    summarize_contrast(
      estimate_interaction_contrast(
        m_augmented,
        contrast_median,
        contrast_below,
        prep_augmented
      ),
      truth,
      "augmented"
    ),
    summarize_contrast(
      estimate_interaction_contrast_mi(
        fit_trunc_mi_list,
        contrast_median,
        contrast_below,
        prep_trunc_mi
      ),
      truth,
      "trunc_mi"
    )
  )

  if (!is.null(m_complete_case)) {
    results <- c(
      results,
      list(
        summarize_contrast(
          estimate_interaction_contrast(
            m_complete_case,
            contrast_median,
            contrast_below,
            prep_complete_case
          ),
          truth,
          "complete_case"
        )
      )
    )
  }

  dplyr::bind_rows(results) |>
    dplyr::mutate(
      contrast = contrast_name,
      moving = metadata_median$moving,
      fixed = metadata_median$fixed,
      fixed_at = paste0(metadata_median$fixed_at, "_minus_", metadata_below$fixed_at),
      z_low = metadata_median$z_low,
      z_high = metadata_median$z_high,
      n_fixed_value = metadata_median$n_fixed_value + metadata_below$n_fixed_value,
      quantile_type = metadata_median$quantile_type,
      .before = method
    )
}


#### Big Wrapper for multiple HPC runs ####
for (seed in seed_vec){
  #### Simulation ####
  print(seed)
  set.seed(seed)
  mi_seed <- seed

  ##### Exposure and LoD#####
  #LoD is based on quantile of true distribution
  #Alternative is to censor x% of data but not current use 
  #lod <- apply(Z_true, 2, quantile, probs = lod_quantile)
  if (exposure_dist == "norm") {
    norm1 <- 0
    norm2 <- 1
    Z_true <- matrix(exp(rnorm(n * p,mean = norm1,sd = norm2)), ncol = p)
    lod <- exp(qnorm(lod_quantile, mean = norm1, sd = norm2))
    lod <- rep(lod, p)

  } else if (exposure_dist == "unif") {
    unif1 = -3
    unif2 = 3
    Z_true <- matrix(exp(runif(n * p, min = unif1, max = unif2)), ncol = p)
    lod <- exp(qunif(lod_quantile, min = unif1, max = unif2))
    lod <- rep(lod, p)

  } else {
    stop("exposure_dist must be 'norm' or 'unif'")
  }

  # Create observed data with censoring
  Z_obs <- CensorData(Z_true, lod)

  #get indices for complete case, can be empty with higher % censoring
  complete_case_idx <- complete.cases(Z_obs)

  ##### Response #####
  h_true <- true_h(Z_true, h_func, mean_offset, scale)
  X <- matrix(rnorm(n * q_fixed_effects), nrow = n)
  X <- scale(X, center = TRUE, scale = FALSE)
  colnames(X) <- paste0("x", seq_len(q_fixed_effects))
  fixed_effect <- drop(X %*% fixed_effect_beta)

  y <- h_true + fixed_effect + rnorm(n, sd = 1)
  y_complete_case <- y[complete_case_idx]
  X_complete_case <- X[complete_case_idx, , drop = FALSE]

  ##### Models #####
  # A. Uncensored
  # B. Imputation (LoD / sqrt(2))
  # C. Augmented (Indicator + Continuous)
  # D. Complete case
  # E. Truncated MI (tobit lognormal)

  ###### A. Uncensored ######
  Z_uncensored <- Z_true

  ###### B. Imputation (LoD / sqrt(2)) ######
  Z_impute <- SingleImputation(Z_obs, lod)

  ###### C. Augmented (Indicator + Continuous) ######
  aug_data <- AugmentData(Z_obs = Z_obs,lod = lod)
  Z_aug <- aug_data$Z

  ###### D. Complete case ######
  Z_complete_case <- Z_uncensored[complete_case_idx, , drop = FALSE]

  ###### E. Truncated multiple imputation using qgcomp::mice.impute.leftcenslognorm ######
  mice.impute.leftcenslognorm <- qgcomp::mice.impute.leftcenslognorm
  x_names <- paste0("x", seq_len(q_fixed_effects))
  z_names <- paste0("z", seq_len(p))
  mi_data <- cbind(y, X, Z_obs) |> as.data.frame()
  colnames(mi_data) <- c("y", x_names, z_names)

  #call with no iterations to get default settings
  mi_init <- mice(mi_data, maxit = 0, printFlag = FALSE)

  method_trunc <- mi_init$method
  z_cols <- (q_fixed_effects + 2):(q_fixed_effects + p + 1)
  method_trunc[z_cols] <- "leftcenslognorm"

  predictor_matrix_trunc <- matrix(0, nrow = ncol(mi_data), ncol = ncol(mi_data))
  predictor_matrix_trunc[z_cols, c(1, seq_len(q_fixed_effects) + 1)] <- 1
  colnames(predictor_matrix_trunc) <- colnames(mi_data)
  rownames(predictor_matrix_trunc) <- colnames(mi_data)

  mids_trunc <- mice(
    data = mi_data,
    m = m_imputations,
    maxit = mi_maxit,
    method = method_trunc,
    predictorMatrix = predictor_matrix_trunc,
    lod = c(rep(NA, q_fixed_effects + 1), lod),
    seed = mi_seed,
    printFlag = FALSE
  )

  Z_trunc_mi_raw_list <- complete(mids_trunc, action = "all") |>
    lapply(function(dat) {
      dat[, z_names, drop = FALSE] |>
        as.matrix()
    })

  log_obs_for_scaling <- log(Z_obs)
  trunc_mi_center <- colMeans(log_obs_for_scaling, na.rm = TRUE)
  trunc_mi_scale <- apply(log_obs_for_scaling, 2, sd, na.rm = TRUE)

  ###### Scaling ######
  #augmented continuous scalet seperately due to the LoD adjustment
  #multiple imputation tobit also needs to be handled separately, done in prevous block
  #common scaling for all MI data sets is used to ensure that they are on the same scale for pooling results 

  ####### uncensored scaling#######
  uncens_center <- colMeans(log(Z_true))
  uncens_scale <- apply(log(Z_true), 2, sd)

  Z_uncensored <- scale( log(Z_true), center = uncens_center, scale = uncens_scale)

  ####### SI scaling #######

  impute_center <- colMeans(log(Z_impute))
  impute_scale <- apply(log(Z_impute), 2, sd)

  Z_impute <- scale(log(Z_impute), center = impute_center, scale = impute_scale)


  ####### CC scaling #######

  if (nrow(Z_complete_case) >= 2) {
    complete_center <- colMeans(log(Z_complete_case))
    complete_scale <- apply(log(Z_complete_case), 2, sd)

    Z_complete_case <- scale(
      log(Z_complete_case),
      center = complete_center,
      scale = complete_scale
    )
  }


  ##### BKMR ##### 
  #used to store results by how many values are below the LoD
  #ie group 0 is all values above the LoD, which is our focus
  group <- rowSums(is.na(Z_obs))

  group_complete_case <- group[complete_case_idx]

  ###############################################################################
  #Model calls and MSE calculation


  ###### BKMR Uncensored ######
  m_uncensored <- kmbayes(y = y, Z = Z_uncensored, X = X, iter = mcmc_iter,varsel = TRUE)
  pred_uncensored <- SamplePred(m_uncensored, Znew = Z_uncensored, Xnew = zero_fixed_effects(nrow(Z_uncensored)))
  results_uncens  <- mse_by_lod_count(h_true, pred_uncensored, group, p)
  results_uncens_first2  <- mse_by_first2_lod(h_true, pred_uncensored, Z_obs)
  coverage_uncens <- coverage_by_lod_count(h_true, pred_uncensored, group, p)
  coverage_uncens_first2 <- coverage_by_first2_lod(h_true, pred_uncensored, Z_obs)
  pip_uncensored <- ExtractPIPs(m_uncensored)
  sensspec_uncensored <- calc_sens_spec(pip_uncensored)

  ###### BKMR Single Imputation ######

  m_impute <- kmbayes(y = y, Z = Z_impute, X = X, iter = mcmc_iter,varsel = TRUE)
  pred_impute <- SamplePred(m_impute, Znew = Z_impute, Xnew = zero_fixed_effects(nrow(Z_impute)))
  results_impute <- mse_by_lod_count(h_true, pred_impute, group, p)
  results_imputes_first2 <- mse_by_first2_lod(h_true, pred_impute, Z_obs)
  coverage_impute <- coverage_by_lod_count(h_true, pred_impute, group, p)
  coverage_impute_first2 <- coverage_by_first2_lod(h_true, pred_impute, Z_obs)
  pip_impute <- ExtractPIPs(m_impute)
  sensspec_impute <- calc_sens_spec(pip_impute)


  ###### BKMR Missing Indicator Method ######



  m_augmented <- kmbayes(y = y, Z = Z_aug, X = X, iter = mcmc_iter,varsel = TRUE)
  pred_augmented <- SamplePred(m_augmented, Znew = Z_aug, Xnew = zero_fixed_effects(nrow(Z_aug)))
  results_augmented <- mse_by_lod_count(h_true, pred_augmented, group, p)
  results_augmented_first2 <- mse_by_first2_lod(h_true, pred_augmented, Z_obs)
  coverage_augmented <- coverage_by_lod_count(h_true, pred_augmented, group, p)
  coverage_augmented_first2 <- coverage_by_first2_lod(h_true, pred_augmented, Z_obs)


  augmented_pips <- extract_augmented_chemical_pips(m_augmented, p = 4)

  pip_augmented_and <- pip_uncensored
  pip_augmented_and[,2] <- augmented_pips$pip_and
  sensspec_augmented_and <- calc_sens_spec(pip_augmented_and)

  pip_augmented_or <- pip_uncensored
  pip_augmented_or[,2] <- augmented_pips$pip_or
  sensspec_augmented_or <- calc_sens_spec(pip_augmented_or)

  ###### BKMR Complete Case ######

  n_complete <- nrow(Z_complete_case)

  if (n_complete >= 2) {
    m_complete_case <- kmbayes( y = y_complete_case, Z = Z_complete_case, X = X_complete_case, iter = mcmc_iter,varsel = TRUE)
    pred_complete_case <- SamplePred( m_complete_case, Znew = Z_complete_case, Xnew = zero_fixed_effects(nrow(Z_complete_case)))
    results_complete_case <- mse_by_lod_count(h_true[complete_case_idx],pred_complete_case,group_complete_case,p)
    results_complete_cases_first2 <- mse_by_first2_lod( h_true[complete_case_idx], pred_complete_case, Z_obs[complete_case_idx, , drop = FALSE])
    coverage_complete_case <- coverage_by_lod_count(
      h_true[complete_case_idx],
      pred_complete_case,
      group_complete_case,
      p
    )
    coverage_complete_case_first2 <- coverage_by_first2_lod(
      h_true[complete_case_idx],
      pred_complete_case,
      Z_obs[complete_case_idx, , drop = FALSE]
    )
    pip_complete_case <- ExtractPIPs(m_complete_case)
    sensspec_complete_case <- calc_sens_spec(pip_complete_case)

  } else {
    m_complete_case <- NULL
    pred_complete_case <- NULL

    results_complete_case <- empty_mse_by_lod_count(p)
    results_complete_cases_first2 <- empty_mse_by_first2_lod()
    coverage_complete_case <- empty_coverage_by_lod_count(p)
    coverage_complete_case_first2 <- empty_coverage_by_first2_lod()

    pip_complete_case <- pip_uncensored
    pip_complete_case[,2] <- NA
    sensspec_complete_case <- sensspec_uncensored
    sensspec_complete_case[,2:3] <- NA
  }


  ###### BKMR MI Tobit ######
  #create empty lists to store results for MI 
  Z_trunc_mi_list <- vector("list", length = m_imputations)
  fit_trunc_mi_list <- vector("list", length = m_imputations)
  pred_trunc_mi_list <- vector("list", length = m_imputations)
  pip_trunc_mi_list <- vector("list", length = m_imputations)

  for (m in seq_len(m_imputations)) {
    Z_trunc_mi_list[[m]] <- scale( log(Z_trunc_mi_raw_list[[m]]), center = trunc_mi_center, scale = trunc_mi_scale)
    fit_trunc_mi_list[[m]] <- kmbayes(y = y,Z = Z_trunc_mi_list[[m]],X = X,iter = mcmc_iter,varsel = TRUE)
    pip_trunc_mi_list[[m]] <- ExtractPIPs(fit_trunc_mi_list[[m]])[,2]

    pred_trunc_mi_list[[m]] <- SamplePred(fit_trunc_mi_list[[m]],Znew = Z_trunc_mi_list[[m]],Xnew = zero_fixed_effects(nrow(Z_trunc_mi_list[[m]])))
  }

  #used to get formatting
  pip_trunc_mi <- pip_uncensored
  pip_trunc_mi[,2] <- Reduce("+", pip_trunc_mi_list) / m_imputations
  sensspec_trunc_mi <- calc_sens_spec(pip_trunc_mi)

  pred_trunc_mi <- do.call(rbind, pred_trunc_mi_list)

  results_trunc_mi <- mse_by_lod_count( h_true, pred_trunc_mi, group, p)
  results_trunc_mi_first2 <- mse_by_first2_lod( h_true, pred_trunc_mi, Z_obs)
  coverage_trunc_mi <- coverage_by_lod_count(h_true, pred_trunc_mi, group, p)
  coverage_trunc_mi_first2 <- coverage_by_first2_lod(h_true, pred_trunc_mi, Z_obs)

  fixed_effect_estimates <- dplyr::bind_rows(
    summarize_fixed_effects(m_uncensored, "uncensored"),
    summarize_fixed_effects(m_impute, "impute"),
    summarize_fixed_effects(m_augmented, "augmented"),
    summarize_fixed_effects(m_complete_case, "complete_case"),
    summarize_fixed_effects_mi(fit_trunc_mi_list, "trunc_mi")
  )


  ##### Contrast estimation #####
  contrast_z1_given_z2 <- make_one_direction_contrast(
    Z_true = Z_true,
    lod = lod,
    moving = 1,
    fixed = 2,
    fixed_at = "below_lod"
  )

  contrast_z2_given_z1 <- make_one_direction_contrast(
    Z_true = Z_true,
    lod = lod,
    moving = 2,
    fixed = 1,
    fixed_at = "below_lod"
  )

  contrast_z1_given_z2_median <- make_one_direction_contrast(
    Z_true = Z_true,
    lod = lod,
    moving = 1,
    fixed = 2,
    fixed_at = "median_above_lod"
  )

  contrast_z2_given_z1_median <- make_one_direction_contrast(
    Z_true = Z_true,
    lod = lod,
    moving = 2,
    fixed = 1,
    fixed_at = "median_above_lod"
  )

  contrast_results <- dplyr::bind_rows(
    estimate_all_methods_for_contrast(
      contrast_z1_given_z2,
      "z1_75_vs_25_given_z2_below"
    ),
    estimate_all_methods_for_contrast(
      contrast_z2_given_z1,
      "z2_75_vs_25_given_z1_below"
    ),
    estimate_all_methods_for_contrast(
      contrast_z1_given_z2_median,
      "z1_75_vs_25_given_z2_median_above_lod"
    ),
    estimate_all_methods_for_contrast(
      contrast_z2_given_z1_median,
      "z2_75_vs_25_given_z1_median_above_lod"
    ),
    estimate_all_methods_for_interaction(
      contrast_z1_given_z2_median,
      contrast_z1_given_z2,
      "interaction_z1_median_above_lod_minus_below_lod"
    ),
    estimate_all_methods_for_interaction(
      contrast_z2_given_z1_median,
      contrast_z2_given_z1,
      "interaction_z2_median_above_lod_minus_below_lod"
    )
  )

  contrast_results

  ##### Results Compilation and Saving #####

  sim_results <- list(
    settings = list(
      seed = seed,
      n = n,
      p = p,
      lod_quantile = lod_quantile,
      exposure_dist = exposure_dist,
      mean_offset = mean_offset,
      scale = scale,
      h_func = h_func,
      mcmc_iter = mcmc_iter,
      m_imputations = m_imputations,
      mi_maxit = mi_maxit,
      q_fixed_effects = q_fixed_effects,
      fixed_effect_beta = fixed_effect_beta
    ),
    logistics = list(
      run_time = Sys.time() - start_time,
      run_mem = sum(gc()[, 6])
    ),
    results = list(
      mse_by_lod_count = list(
        uncensored = results_uncens,
        impute = results_impute,
        augmented = results_augmented,
        complete_case = results_complete_case,
        trunc_mi = results_trunc_mi
      ),
      mse_by_first2_lod = list(
        uncensored = results_uncens_first2,
        impute = results_imputes_first2,
        augmented = results_augmented_first2,
        complete_case = results_complete_cases_first2,
        trunc_mi = results_trunc_mi_first2
      ),
      coverage_by_lod_count = list(
        uncensored = coverage_uncens,
        impute = coverage_impute,
        augmented = coverage_augmented,
        complete_case = coverage_complete_case,
        trunc_mi = coverage_trunc_mi
      ),
      coverage_by_first2_lod = list(
        uncensored = coverage_uncens_first2,
        impute = coverage_impute_first2,
        augmented = coverage_augmented_first2,
        complete_case = coverage_complete_case_first2,
        trunc_mi = coverage_trunc_mi_first2
      ),
      pips = list(
        uncensored = pip_uncensored,
        impute = pip_impute,
        augmented_and = pip_augmented_and,
        augmented_or = pip_augmented_or,
        complete_case = pip_complete_case,
        trunc_mi = pip_trunc_mi
      ),
      sensspec = list(
        uncensored = sensspec_uncensored,
        impute = sensspec_impute,
        augmented_and = sensspec_augmented_and,
        augmented_or = sensspec_augmented_or,
        complete_case = sensspec_complete_case,
        trunc_mi = sensspec_trunc_mi
      ),
      fixed_effects = fixed_effect_estimates,
      contrasts = contrast_results
    )
  )

  name <- paste0(
    "SimResults_seed", seed,
    "_n", n,
    "_lod", lod_quantile,
    "_", exposure_dist,
    "_mo", mean_offset,
    "_s", scale,
    "_h",h_func,
    ".rds"
  )

  saveRDS(sim_results, file = name)

}











