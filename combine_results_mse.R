library(dplyr)
library(purrr)
library(tibble)
library(tidyr)

results_dir <- "Routputs"
result_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(result_files) == 0) {
  stop("No .rds files found in Results/.")
}

counter  <- 0
read_sim_result <- function(path) {
  counter  <<- counter  +1
  if (counter %% 500 ==0){print(counter)}
  sim_result <- readRDS(path)

  settings_tbl <- as_tibble(sim_result$settings)

  logistics_tbl <- tibble(
    run_time = as.numeric(sim_result$logistics$run_time, units = "mins"),
    run_mem = sim_result$logistics$run_mem,
    file = basename(path)
  )

  metadata_tbl <- bind_cols(settings_tbl, logistics_tbl)

  mse_by_lod_count <- imap_dfr(
    sim_result$results$mse_by_lod_count,
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

  mse_by_first2_lod <- imap_dfr(
    sim_result$results$mse_by_first2_lod,
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

  coverage_by_lod_count <- if (is.null(sim_result$results$coverage_by_lod_count)) {
    tibble()
  } else {
    imap_dfr(
      sim_result$results$coverage_by_lod_count,
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

  coverage_by_first2_lod <- if (is.null(sim_result$results$coverage_by_first2_lod)) {
    tibble()
  } else {
    imap_dfr(
      sim_result$results$coverage_by_first2_lod,
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

  sensspec <- imap_dfr(
    sim_result$results$sensspec,
    \(result_tbl, method) {
      bind_cols(
        metadata_tbl,
        as_tibble(result_tbl) |>
          mutate(method = method)
      )
    }
  )

  pips <- imap_dfr(
    sim_result$results$pips,
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

  list(
    file_metadata = metadata_tbl,
    mse_by_lod_count = mse_by_lod_count,
    mse_by_first2_lod = mse_by_first2_lod,
    coverage_by_lod_count = coverage_by_lod_count,
    coverage_by_first2_lod = coverage_by_first2_lod,
    sensspec = sensspec,
    pips = pips,
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
      id_cols = c(n, p, lod_quantile, exposure_dist, mean_offset, h_func, scale, mcmc_iter, group),
      names_from = method,
      values_from = pooled_mse,
      names_prefix = "pooled_mse_"
    )
}

pivot_sensspec_wider <- function(data) {
  data |>
    pivot_wider(
      id_cols = c(n, p, lod_quantile, exposure_dist, mean_offset, h_func, scale, mcmc_iter, threshold),
      names_from = method,
      values_from = c(mean_sensitivity, mean_specificity),
      names_sep = "_"
    )
}

combined_raw <- map(result_files, read_sim_result)
combined_file_metadata <- map_dfr(combined_raw, "file_metadata")
combined_mse_by_lod_count <- map_dfr(combined_raw, "mse_by_lod_count")
combined_mse_by_first2_lod <- map_dfr(combined_raw, "mse_by_first2_lod")
combined_coverage_by_lod_count <- map_dfr(combined_raw, "coverage_by_lod_count")
combined_coverage_by_first2_lod <- map_dfr(combined_raw, "coverage_by_first2_lod")
combined_sensspec <- map_dfr(combined_raw, "sensspec")
combined_pips <- map_dfr(combined_raw, "pips")
combined_contrasts <- map_dfr(combined_raw, "contrasts") |>
  add_contrast_family()

scenario_cols <- c(
  "n",
  "p",
  "lod_quantile",
  "exposure_dist",
  "mean_offset",
  "h_func",
  "scale",
  "mcmc_iter"
)

mse_by_lod_count_summary <- pool_mse(
  combined_mse_by_lod_count,
  c(scenario_cols, "method", "group")
) |>
  pivot_pooled_mse_wider()

mse_by_first2_lod_summary <- pool_mse(
  combined_mse_by_first2_lod,
  c(scenario_cols, "method", "group")
) |>
  pivot_pooled_mse_wider()

coverage_by_lod_count_summary <- pool_coverage(
  combined_coverage_by_lod_count,
  c(scenario_cols, "method", "group")
)

coverage_by_first2_lod_summary <- pool_coverage(
  combined_coverage_by_first2_lod,
  c(scenario_cols, "method", "group")
)

ci_width_by_lod_count_summary <- if (nrow(coverage_by_lod_count_summary) == 0) {
  tibble()
} else {
  coverage_by_lod_count_summary |>
    select(
      all_of(scenario_cols),
      method,
      group,
      runs,
      total_n_obs,
      total_covered,
      pooled_mean_ci_width = mean_ci_width,
      empirical_coverage,
      empirical_coverage_pct
    )
}

ci_width_by_first2_lod_summary <- if (nrow(coverage_by_first2_lod_summary) == 0) {
  tibble()
} else {
  coverage_by_first2_lod_summary |>
    select(
      all_of(scenario_cols),
      method,
      group,
      runs,
      total_n_obs,
      total_covered ,
      pooled_mean_ci_width = mean_ci_width,
      empirical_coverage,
      empirical_coverage_pct
    )
}

sensspec_summary_long <- pool_sensspec(
  combined_sensspec,
  c(scenario_cols, "method", "threshold")
)

sensspec_summary <- sensspec_summary_long |>
  pivot_sensspec_wider()

contrast_summary <- pool_contrasts(
  combined_contrasts,
  c(scenario_cols, "contrast_family", "fixed_at", "method")
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
  combined_mse_by_lod_count |>
    filter(method != "complete_case"),
  c(scenario_cols)
) |>
  select(
    n,
    lod_quantile,
    exposure_dist,
    mean_offset,
    h_func,
    scale,
    runs
  )

combined_results <- list(
  files = combined_file_metadata,
  combined_mse_by_lod_count = combined_mse_by_lod_count,
  combined_mse_by_first2_lod = combined_mse_by_first2_lod,
  combined_coverage_by_lod_count = combined_coverage_by_lod_count,
  combined_coverage_by_first2_lod = combined_coverage_by_first2_lod,
  combined_sensspec = combined_sensspec,
  combined_pips = combined_pips,
  combined_contrasts = combined_contrasts,
  mse_by_lod_count_summary = mse_by_lod_count_summary,
  mse_by_first2_lod_summary = mse_by_first2_lod_summary,
  coverage_by_lod_count_summary = coverage_by_lod_count_summary,
  coverage_by_first2_lod_summary = coverage_by_first2_lod_summary,
  ci_width_by_lod_count_summary = ci_width_by_lod_count_summary,
  ci_width_by_first2_lod_summary = ci_width_by_first2_lod_summary,
  sensspec_summary_long = sensspec_summary_long,
  sensspec_summary = sensspec_summary,
  contrast_summary = contrast_summary
)


process_df <- combined_results$mse_by_first2_lod_summary |> 
  filter(
    # h_func == 3
    # group == "-+" | group == "+-",
    group == "++",
    # group == '0',
    # is.na(group),
  ) |>
  select(
    n,
    lod_quantile,
    exposure_dist,
    mean_offset,
    scale,
    # group,
    h_func,
    any_of(c(
      "pooled_mse_uncensored",
      "pooled_mse_impute",
      "pooled_mse_augmented",
      "pooled_mse_trunc_mi"
      # "pooled_mse_complete_case"
    ))
  )

break
library(looplot)
colnames(process_df)[(ncol(process_df) - 3):ncol(process_df)] <- c("Oracle", "Single Imputation", "Augmented", "Truncated MI")


# plot_data = nested_loop_base_data(
#     process_df, 
#     x = "n", steps = c("group","exposure_dist","scale"),
#     grid_cols = "mean_offset", grid_rows = "lod_quantile",
#     spu_x_shift = .2 * 1000
# )


plot_data = nested_loop_base_data(
    process_df, 
    x = "lod_quantile", steps = c("mean_offset","scale","n"),
    grid_cols = "exposure_dist", grid_rows = "h_func",
    spu_x_shift = .2
)

plot_data = nested_loop_paramsteps_data(
    plot_data,
    steps_y_base = -.1,
    steps_y_height = .05,
    steps_y_shift = .1
)

p = nested_loop_base_plot(
    plot_data,
    x_name = "Detection Limit Quantile",
    y_name = "RMSE", 
    colors = scales::viridis_pal(end = .85, option = "A"),
    grid_scales = "free_y"
)

p = nested_loop_paramsteps_plot(
    p, plot_data, 
    steps_values_annotate = TRUE, 
    steps_annotation_size = 3
)

p = add_processing(
    p, 
    list(
        # set limits
        adjust_ylim = list(
            y_expand_add = c(.25, NULL)
        ),
        # adjust theme
        add_custom_theme = list(
            axis.text.x = element_text(angle = 90, 
                                       vjust = 0.5, 
                                       size = 7)
        ), 
        # add horizontal lines
        add_abline = list(
            intercept = 0
        )
    )
)
print(p)
ggsave("RMSE++.png",width = 12,height = 12)

##############################

# R
sens_plot_df <- combined_results$sensspec_summary_long |>
  select(
    n,
    lod_quantile,
    exposure_dist,
    mean_offset,
    h_func,
    scale,
    threshold,
    method,
    mean_sensitivity,
    mean_specificity
  ) |>
  pivot_longer(
    cols = c(mean_sensitivity, mean_specificity),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = recode(
      metric,
      mean_sensitivity = "sensitivity",
      mean_specificity = "specificity"
    )
  ) |>
  pivot_wider(
    names_from = method,
    values_from = value
  )


sens_plot_df <- sens_plot_df |> 
  filter(metric=="sensitivity") |> 
  # filter(metric=="specificity") |> 
  filter(threshold==0.75 | threshold==0.9) |> 
  # filter(h_func==3) |> 
  select(
  n,
  lod_quantile,
  exposure_dist,
  mean_offset,
  scale,
  h_func,
  # metric,
  threshold,
  uncensored,
  impute,
  # augmented_and,
  augmented_or,
  trunc_mi
) 

# colnames(sens_plot_df)[8:12] <- c("Oracle", "Single Imputation", "Augmented and","Augmented or", "Truncated MI")

colnames(sens_plot_df)[(length(colnames(sens_plot_df))-3):length(colnames(sens_plot_df))] <- c("Oracle", "Single Imputation","Augmented or", "Truncated MI")


# plot_data = nested_loop_base_data(
#     sens_plot_df, 
#     x = "n", steps = c("threshold","exposure_dist","metric","scale"),
#     grid_cols = "mean_offset", grid_rows = "lod_quantile",
#     spu_x_shift = .2 * 1000
# )

plot_data = nested_loop_base_data(
    sens_plot_df, 
    x = "lod_quantile", steps = c("threshold","n","scale","mean_offset"),
    grid_cols = "exposure_dist", grid_rows = "h_func",
    spu_x_shift = .3
)


p = nested_loop_base_plot(
    plot_data,
    x_name = "Detection Limit Quantile",
    y_name = "Sensitivity", 
    colors = scales::viridis_pal(end = .85, option = "A")
)

plot_data = nested_loop_paramsteps_data(
    plot_data,
    steps_y_base = 0.6,
    steps_y_height = .05,
    steps_y_shift = .05
)

p = nested_loop_paramsteps_plot(
    p, plot_data, 
    steps_values_annotate = TRUE, 
    steps_annotation_size = 3
)

p = add_processing(
    p, 
    list(
        # set limits
        adjust_ylim = list(
            y_expand_add = c(.25, NULL)
        ),
        # adjust theme
        add_custom_theme = list(
            axis.text.x = element_text(angle = 90, 
                                       vjust = 0.5, 
                                       size = 6)
        )
        # add horizontal lines
        # add_abline = list(
        #     intercept = 0
        # )
    )
)
print(p)
ggsave("Sensitivity.png",width = 16,height = 12)


