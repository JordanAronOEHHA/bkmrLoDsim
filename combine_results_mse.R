library(dplyr)
library(purrr)
library(tibble)
library(tidyr)

results_dir <- "Routputs"
result_files <- list.files(
  results_dir,
  pattern = "\\.rds$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(result_files) == 0) {
  stop("No .rds files found in Results/.")
}

counter  <- 0

prediction_methods <- c("uncensored", "impute", "augmented", "trunc_mi")
stored_result_methods <- c(prediction_methods, "augmented_and", "augmented_or")

filter_result_methods <- function(result_list) {
  if (is.null(result_list)) {
    return(NULL)
  }

  result_list[intersect(names(result_list), stored_result_methods)]
}

get_result_table <- function(results, name, fallback = NULL) {
  if (!is.null(results[[name]])) {
    return(results[[name]])
  }

  if (!is.null(fallback) && !is.null(results[[fallback]])) {
    return(results[[fallback]])
  }

  NULL
}

bind_method_tables <- function(result_list, metadata_tbl) {
  if (is.null(result_list)) {
    return(tibble())
  }

  result_list <- filter_result_methods(result_list)

  imap_dfr(
    result_list,
    \(result_tbl, method) {
      bind_cols(
        metadata_tbl,
        as_tibble(result_tbl) |>
          mutate(
            method = method
          )
      )
    }
  )
}

prediction_diagnostic_cols <- c(
  "group",
  "n_obs",
  "mse",
  "rmse",
  "mean_bias",
  "mae",
  "sd_error",
  "mean_posterior_sd",
  "calibration_ratio",
  "n_covered",
  "empirical_coverage",
  "n_lower_misses",
  "n_upper_misses",
  "lower_miss_rate",
  "upper_miss_rate",
  "mean_ci_width",
  "mean_interval_score"
)

fixed_effect_estimate_cols <- c(
  "method",
  "fixed_effect",
  "truth",
  "est",
  "bias",
  "ci_low",
  "ci_high",
  "covered"
)

standardize_prediction_diagnostics <- function(data) {
  if (nrow(data) == 0 && ncol(data) == 0) {
    return(tibble())
  }

  missing_cols <- setdiff(prediction_diagnostic_cols, names(data))
  for (col in missing_cols) {
    data[[col]] <- NA
  }

  data |>
    mutate(
      rmse = ifelse(is.na(rmse) & !is.na(mse), sqrt(mse), rmse)
    )
}

standardize_fixed_effect_estimates <- function(data) {
  data <- as_tibble(data)

  fixed_effect_col_defaults <- list(
    method = NA_character_,
    fixed_effect = NA_character_,
    truth = NA_real_,
    est = NA_real_,
    bias = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    covered = NA
  )

  missing_cols <- setdiff(fixed_effect_estimate_cols, names(data))
  for (col in missing_cols) {
    data[[col]] <- rep(fixed_effect_col_defaults[[col]], nrow(data))
  }

  data
}

legacy_prediction_diagnostics <- function(mse_list, coverage_list) {
  if (is.null(mse_list) && is.null(coverage_list)) {
    return(NULL)
  }

  method_names <- union(names(mse_list), names(coverage_list))

  stats::setNames(
    lapply(method_names, function(method) {
      mse_tbl <- if (!is.null(mse_list[[method]])) {
        as_tibble(mse_list[[method]])
      } else {
        tibble()
      }

      coverage_tbl <- if (!is.null(coverage_list[[method]])) {
        as_tibble(coverage_list[[method]])
      } else {
        tibble()
      }

      out <- if (nrow(mse_tbl) == 0) {
        coverage_tbl
      } else if (nrow(coverage_tbl) == 0) {
        mse_tbl
      } else {
        full_join(
          mse_tbl,
          coverage_tbl,
          by = "group",
          suffix = c("_mse", "_coverage")
        ) |>
          mutate(
            n_obs = coalesce(n_obs_mse, n_obs_coverage)
          ) |>
          select(-any_of(c("n_obs_mse", "n_obs_coverage")))
      }

      standardize_prediction_diagnostics(out)
    }),
    method_names
  )
}

read_sim_result <- function(path) {
  counter  <<- counter  +1
  if (counter %% 500 ==0){print(counter)}
  sim_result <- readRDS(path)

  vector_setting_names <- c("active_idx", "fixed_effect_beta")
  settings <- purrr::imap(sim_result$settings, function(x, name) {
    if (length(x) == 0) {
      return(NA)
    }

    if (name %in% vector_setting_names) {
      return(paste(x, collapse = ","))
    }

    if (length(x) == 1) {
      return(x)
    }

    paste(x, collapse = ",")
  })

  if (is.null(settings$correlation)) {
    settings$correlation <- "ind"
  }

  settings_tbl <- as_tibble(settings)

  logistics_tbl <- tibble(
    run_time = as.numeric(sim_result$logistics$run_time, units = "mins"),
    run_mem = sim_result$logistics$run_mem,
    file = basename(path)
  )

  metadata_tbl <- bind_cols(settings_tbl, logistics_tbl)

  active_prediction_diagnostics <- get_result_table(
    sim_result$results,
    "prediction_diagnostics_by_active_lod_burden"
  )

  if (is.null(active_prediction_diagnostics)) {
    active_mse <- get_result_table(
      sim_result$results,
      "mse_by_active_lod_burden",
      fallback = "mse_by_first2_lod"
    )

    active_coverage <- get_result_table(
      sim_result$results,
      "coverage_by_active_lod_burden",
      fallback = "coverage_by_first2_lod"
    )

    active_prediction_diagnostics <- legacy_prediction_diagnostics(
      active_mse,
      active_coverage
    )
  }

  prediction_diagnostics_by_active_lod_burden <- bind_method_tables(
    active_prediction_diagnostics,
    metadata_tbl
  ) |>
    standardize_prediction_diagnostics()

  mse_by_active_lod_burden <- prediction_diagnostics_by_active_lod_burden |>
    select(any_of(c(
      names(metadata_tbl),
      "method",
      "group",
      "n_obs",
      "mse"
    )))

  coverage_by_active_lod_burden <- prediction_diagnostics_by_active_lod_burden |>
    select(any_of(c(
      names(metadata_tbl),
      "method",
      "group",
      "n_obs",
      "n_covered",
      "empirical_coverage",
      "mean_ci_width"
    )))

  sensspec <- imap_dfr(
    filter_result_methods(sim_result$results$sensspec),
    \(result_tbl, method) {
      bind_cols(
        metadata_tbl,
        as_tibble(result_tbl) |>
          mutate(method = method)
      )
    }
  )

  pips <- imap_dfr(
    filter_result_methods(sim_result$results$pips),
    \(result_tbl, method) {
      out <- as_tibble(result_tbl)

      if ("chemical" %in% names(out) && "pip" %in% names(out)) {
        out <- out |>
          rename(PIP = pip)
      }

      bind_cols(
        metadata_tbl,
        out |>
          mutate(method = method)
      )
    }
  )

  contrasts <- if (is.null(sim_result$results$contrasts)) {
    tibble()
  } else {
    bind_cols(
      metadata_tbl,
      as_tibble(sim_result$results$contrasts)
    )
  }

  fixed_effect_estimates <- get_result_table(
    sim_result$results,
    "fixed_effects",
    fallback = "fixed_effect_estimates"
  )

  fixed_effect_estimates <- if (is.null(fixed_effect_estimates)) {
    tibble()
  } else {
    bind_cols(
      metadata_tbl,
      as_tibble(fixed_effect_estimates)
    )
  } |>
    standardize_fixed_effect_estimates()

  list(
    file_metadata = metadata_tbl,
    prediction_diagnostics_by_active_lod_burden = prediction_diagnostics_by_active_lod_burden,
    mse_by_active_lod_burden = mse_by_active_lod_burden,
    coverage_by_active_lod_burden = coverage_by_active_lod_burden,
    sensspec = sensspec,
    pips = pips,
    fixed_effect_estimates = fixed_effect_estimates,
    contrasts = contrasts
  )
}

pool_mse <- function(data, group_cols) {
  data |>
    mutate(sse = mse * n_obs) |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      runs = n_distinct(seed),
      total_n_obs = sum(n_obs, na.rm = TRUE),
      pooled_mse = sqrt(sum(sse, na.rm = TRUE) / total_n_obs),
      .groups = "drop"
    ) |>
    arrange(across(all_of(group_cols)))
}

pool_prediction_diagnostics <- function(data, group_cols) {
  if (nrow(data) == 0) {
    return(tibble())
  }

  weighted_mean <- function(value, weight) {
    valid <- !is.na(value) & !is.na(weight) & weight > 0
    if (!any(valid)) {
      return(NA_real_)
    }

    sum(value[valid] * weight[valid]) / sum(weight[valid])
  }

  sum_count <- function(count) {
    valid <- !is.na(count)
    if (!any(valid)) {
      return(NA_real_)
    }

    sum(count[valid])
  }

  count_rate <- function(count, weight) {
    valid <- !is.na(count) & !is.na(weight) & weight > 0
    if (!any(valid)) {
      return(NA_real_)
    }

    sum(count[valid]) / sum(weight[valid])
  }

  pooled_sd_error <- function(mean_bias, mse, weight) {
    valid <- !is.na(mean_bias) & !is.na(mse) & !is.na(weight) & weight > 0
    total_n <- sum(weight[valid])
    if (total_n <= 1) {
      return(NA_real_)
    }

    sum_error <- sum(mean_bias[valid] * weight[valid])
    sum_sq_error <- sum(mse[valid] * weight[valid])
    variance <- (sum_sq_error - (sum_error^2 / total_n)) / (total_n - 1)

    sqrt(max(variance, 0))
  }

  data |>
    mutate(
      source_n_obs = n_obs,
      source_mse = mse,
      source_mean_bias = mean_bias,
      source_mae = mae,
      source_mean_posterior_sd = mean_posterior_sd,
      source_mean_ci_width = mean_ci_width,
      source_mean_interval_score = mean_interval_score
    ) |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      runs = n_distinct(seed),
      total_n_obs = sum(source_n_obs, na.rm = TRUE),
      total_covered = sum_count(n_covered),
      total_lower_misses = sum_count(n_lower_misses),
      total_upper_misses = sum_count(n_upper_misses),
      mse = weighted_mean(source_mse, source_n_obs),
      mean_bias = weighted_mean(source_mean_bias, source_n_obs),
      mae = weighted_mean(source_mae, source_n_obs),
      sd_error = pooled_sd_error(source_mean_bias, source_mse, source_n_obs),
      mean_posterior_sd = weighted_mean(source_mean_posterior_sd, source_n_obs),
      empirical_coverage = count_rate(n_covered, source_n_obs),
      empirical_coverage_pct = 100 * empirical_coverage,
      lower_miss_rate = count_rate(n_lower_misses, source_n_obs),
      upper_miss_rate = count_rate(n_upper_misses, source_n_obs),
      mean_ci_width = weighted_mean(source_mean_ci_width, source_n_obs),
      mean_interval_score = weighted_mean(source_mean_interval_score, source_n_obs),
      .groups = "drop"
    ) |>
    mutate(
      rmse = sqrt(mse),
      calibration_ratio = ifelse(
        mean_posterior_sd > 0,
        sd_error / mean_posterior_sd,
        NA_real_
      ),
      .after = mean_posterior_sd
    ) |>
    arrange(across(all_of(group_cols)))
}

pool_logistics <- function(data, group_cols) {
  data |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      runs = n_distinct(seed),
      max_run_time = max(run_time, na.rm = TRUE),
      max_run_mem = max(run_mem, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(across(all_of(group_cols)))
}

pool_sensspec <- function(data, group_cols) {
  data |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      runs = n_distinct(seed),
      mean_sensitivity = mean(sensitivity, na.rm = TRUE),
      mean_specificity = mean(specificity, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(across(all_of(group_cols)))
}

pool_coverage <- function(data, group_cols) {
  if (nrow(data) == 0) {
    return(tibble())
  }

  data |>
    mutate(
      ci_width_weight = mean_ci_width * n_obs,
      ci_width_n = ifelse(is.na(mean_ci_width), 0L, n_obs)
    ) |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      runs = n_distinct(seed),
      total_n_obs = sum(n_obs, na.rm = TRUE),
      total_covered = sum(n_covered, na.rm = TRUE),
      empirical_coverage = ifelse(
        total_n_obs == 0,
        NA_real_,
        total_covered / total_n_obs
      ),
      empirical_coverage_pct = 100 * empirical_coverage,
      mean_ci_width = ifelse(
        sum(ci_width_n, na.rm = TRUE) == 0,
        NA_real_,
        sum(ci_width_weight, na.rm = TRUE) / sum(ci_width_n, na.rm = TRUE)
      ),
      .groups = "drop"
    ) |>
    arrange(across(all_of(group_cols)))
}

pool_contrasts <- function(data, group_cols) {
  if (nrow(data) == 0) {
    return(tibble())
  }

  data |>
    mutate(ci_width = ci_high - ci_low) |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      runs = n_distinct(seed),
      n_contrast_estimates = dplyr::n(),
      mean_n_fixed_value = mean(n_fixed_value, na.rm = TRUE),
      mean_truth = mean(truth, na.rm = TRUE),
      mean_est = mean(est, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      rmse = sqrt(mean(bias^2, na.rm = TRUE)),
      coverage = mean(covered, na.rm = TRUE),
      mean_ci_width = mean(ci_width, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(across(all_of(group_cols)))
}

pool_fixed_effect_estimates <- function(data, group_cols) {
  if (nrow(data) == 0) {
    return(tibble())
  }

  mean_or_na <- function(value) {
    valid <- !is.na(value)
    if (!any(valid)) {
      return(NA_real_)
    }

    mean(value[valid])
  }

  rate_or_na <- function(value) {
    valid <- !is.na(value)
    if (!any(valid)) {
      return(NA_real_)
    }

    mean(value[valid])
  }

  data |>
    mutate(
      ci_width = ci_high - ci_low,
      covered = as.logical(covered)
    ) |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      runs = n_distinct(seed),
      n_fixed_effect_estimates = dplyr::n(),
      truth = mean_or_na(truth),
      mean_est = mean_or_na(est),
      mean_bias = mean_or_na(bias),
      rmse = sqrt(mean_or_na(bias^2)),
      coverage = rate_or_na(covered),
      coverage_pct = 100 * coverage,
      mean_ci_width = mean_or_na(ci_width),
      .groups = "drop"
    ) |>
    arrange(across(all_of(group_cols)))
}

add_contrast_family <- function(data) {
  if (nrow(data) == 0) {
    return(data)
  }

  data |>
    mutate(
      contrast_family = dplyr::case_when(
        contrast %in% c(
          "z1_75_vs_25_given_z2_below",
          "z2_75_vs_25_given_z1_below"
        ) ~ "75_vs_25_given_other_below",
        contrast %in% c(
          "z1_75_vs_25_given_z2_median_above_lod",
          "z2_75_vs_25_given_z1_median_above_lod"
        ) ~ "75_vs_25_given_other_median_above_lod",
        contrast %in% c(
          "interaction_z1_median_above_lod_minus_below_lod",
          "interaction_z2_median_above_lod_minus_below_lod"
        ) ~ "interaction_median_above_lod_minus_below_lod",
        TRUE ~ contrast
      )
    )
}

pivot_pooled_mse_wider <- function(data) {
  data |>
    pivot_wider(
      id_cols = c(n, p, lod_quantile, exposure_dist, correlation, mean_offset, h_func, scale, mcmc_iter, group),
      names_from = method,
      values_from = pooled_mse,
      names_prefix = "pooled_mse_"
    )
}

pivot_sensspec_wider <- function(data) {
  data |>
    pivot_wider(
      id_cols = c(n, p, lod_quantile, exposure_dist, correlation, mean_offset, h_func, scale, mcmc_iter, threshold),
      names_from = method,
      values_from = c(mean_sensitivity, mean_specificity),
      names_sep = "_"
    )
}

combined_raw <- map(result_files, read_sim_result)
combined_file_metadata <- map_dfr(combined_raw, "file_metadata")
combined_prediction_diagnostics_by_active_lod_burden <- map_dfr(
  combined_raw,
  "prediction_diagnostics_by_active_lod_burden"
)
combined_mse_by_active_lod_burden <- map_dfr(combined_raw, "mse_by_active_lod_burden")
combined_coverage_by_active_lod_burden <- map_dfr(combined_raw, "coverage_by_active_lod_burden")
combined_sensspec <- map_dfr(combined_raw, "sensspec")
combined_pips <- map_dfr(combined_raw, "pips")
combined_fixed_effect_estimates <- map_dfr(combined_raw, "fixed_effect_estimates") |>
  standardize_fixed_effect_estimates() |>
  filter(method %in% prediction_methods)
combined_contrasts <- map_dfr(combined_raw, "contrasts") |>
  filter(method %in% prediction_methods) |>
  add_contrast_family()

scenario_cols <- c(
  "n",
  "p",
  "lod_quantile",
  "exposure_dist",
  "correlation",
  "mean_offset",
  "h_func",
  "scale",
  "mcmc_iter"
)

prediction_diagnostics_by_active_lod_burden_summary <- pool_prediction_diagnostics(
  combined_prediction_diagnostics_by_active_lod_burden,
  c(scenario_cols, "method", "group")
)

sensspec_summary_long <- pool_sensspec(
  combined_sensspec,
  c(scenario_cols, "method", "threshold")
)

contrast_summary <- pool_contrasts(
  combined_contrasts,
  c(scenario_cols, "contrast_family", "fixed_at", "method")
)

fixed_effect_estimates_summary <- pool_fixed_effect_estimates(
  combined_fixed_effect_estimates,
  c(scenario_cols, "method", "fixed_effect")
)

logistics_summary <- pool_logistics(
  combined_file_metadata,
  scenario_cols
)

logistics_summary |> 
  group_by(n) |> 
  summarize(
    max_run_time = max(max_run_time, na.rm = TRUE),
    max_run_mem = max(max_run_mem, na.rm = TRUE)
)

runs_by_method <- pool_mse(
  combined_mse_by_active_lod_burden |>
    filter(group == "Overall"),
  c(scenario_cols)
) |>
  select(
    n,
    lod_quantile,
    exposure_dist,
    correlation,
    mean_offset,
    h_func,
    scale,
    runs
  )

combined_results <- list(
  files = combined_file_metadata,
  combined_prediction_diagnostics_by_active_lod_burden = combined_prediction_diagnostics_by_active_lod_burden,
  combined_mse_by_active_lod_burden = combined_mse_by_active_lod_burden,
  combined_coverage_by_active_lod_burden = combined_coverage_by_active_lod_burden,
  combined_sensspec = combined_sensspec,
  combined_pips = combined_pips,
  combined_fixed_effect_estimates = combined_fixed_effect_estimates,
  combined_contrasts = combined_contrasts,
  prediction_diagnostics_by_active_lod_burden_summary = prediction_diagnostics_by_active_lod_burden_summary,
  sensspec_summary_long = sensspec_summary_long,
  fixed_effect_estimates_summary = fixed_effect_estimates_summary,
  contrast_summary = contrast_summary
)

saveRDS(combined_results, file = "combined_results.rds")
