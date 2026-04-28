library(dplyr)
library(purrr)
library(tibble)
library(tidyr)

results_dir <- "Routputs"
result_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(result_files) == 0) {
  stop("No .rds files found in Results/.")
}

read_sim_result <- function(path) {
  sim_result <- readRDS(path)

  settings_tbl <- as_tibble(sim_result$settings)

  logistics_tbl <- tibble(
    run_time = as.numeric(sim_result$logistics$run_time, units = "mins"),
    run_mem = sim_result$logistics$run_mem,
    file = basename(path)
  )

  metadata_tbl <- bind_cols(settings_tbl, logistics_tbl)

  mse_by_lod_count <- imap_dfr(
    sim_result$results[c("training", "testing")],
    \(split_obj, split_name) {
      imap_dfr(
        split_obj$mse_by_lod_count,
        \(result_tbl, method) {
          bind_cols(
            metadata_tbl,
            as_tibble(result_tbl) |>
              mutate(
                split = split_name,
                method = method
              )
          )
        }
      )
    }
  )

  mse_by_first2_lod <- imap_dfr(
    sim_result$results[c("training", "testing")],
    \(split_obj, split_name) {
      imap_dfr(
        split_obj$mse_by_first2_lod,
        \(result_tbl, method) {
          bind_cols(
            metadata_tbl,
            as_tibble(result_tbl) |>
              mutate(
                split = split_name,
                method = method
              )
          )
        }
      )
    }
  )

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

  list(
    file_metadata = metadata_tbl,
    mse_by_lod_count = mse_by_lod_count,
    mse_by_first2_lod = mse_by_first2_lod,
    sensspec = sensspec,
    pips = pips
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

pivot_pooled_mse_wider <- function(data) {
  data |>
    pivot_wider(
      id_cols = c(n, n_te, p, lod_quantile, exposure_dist, h_dist, mcmc_iter, split, group),
      names_from = method,
      values_from = pooled_mse,
      names_prefix = "pooled_mse_"
    )
}

pivot_sensspec_wider <- function(data) {
  data |>
    pivot_wider(
      id_cols = c(n, n_te, p, lod_quantile, exposure_dist, h_dist, mcmc_iter, threshold),
      names_from = method,
      values_from = c(mean_sensitivity, mean_specificity),
      names_sep = "_"
    )
}

combined_raw <- map(result_files, read_sim_result)

combined_file_metadata <- map_dfr(combined_raw, "file_metadata")
combined_mse_by_lod_count <- map_dfr(combined_raw, "mse_by_lod_count")
combined_mse_by_first2_lod <- map_dfr(combined_raw, "mse_by_first2_lod")
combined_sensspec <- map_dfr(combined_raw, "sensspec")
combined_pips <- map_dfr(combined_raw, "pips")

scenario_cols <- c(
  "n",
  "n_te",
  "p",
  "lod_quantile",
  "exposure_dist",
  "h_dist",
  "mcmc_iter"
)

mse_by_lod_count_summary <- pool_mse(
  combined_mse_by_lod_count,
  c(scenario_cols, "split", "method", "group")
) |>
  pivot_pooled_mse_wider()

mse_by_first2_lod_summary <- pool_mse(
  combined_mse_by_first2_lod,
  c(scenario_cols, "split", "method", "group")
) |>
  pivot_pooled_mse_wider()

sensspec_summary_long <- pool_sensspec(
  combined_sensspec,
  c(scenario_cols, "method", "threshold")
)

sensspec_summary <- sensspec_summary_long |>
  pivot_sensspec_wider()

logistics_summary <- pool_logistics(
  combined_file_metadata,
  scenario_cols
)

runs_by_method <- pool_mse(
  combined_mse_by_lod_count |>
    filter(method != "complete_case"),
  c(scenario_cols)
) |>
  select(
    n,
    n_te,
    lod_quantile,
    exposure_dist,
    h_dist,
    runs
  )

combined_results <- list(
  files = combined_file_metadata,
  combined_mse_by_lod_count = combined_mse_by_lod_count,
  combined_mse_by_first2_lod = combined_mse_by_first2_lod,
  combined_sensspec = combined_sensspec,
  combined_pips = combined_pips,
  mse_by_lod_count_summary = mse_by_lod_count_summary,
  mse_by_first2_lod_summary = mse_by_first2_lod_summary,
  sensspec_summary_long = sensspec_summary_long,
  sensspec_summary = sensspec_summary,
  logistics_summary = logistics_summary,
  runs_by_method = runs_by_method
)


process_df <- combined_results$mse_by_first2_lod_summary |> 
  filter(
    split == "training",
    # group == "-+" | group == "+-"
    group == "++",
    lod_quantile < 0.9
  ) |>
  select(
    n,
    lod_quantile,
    exposure_dist,
    h_dist,
    # group,
    any_of(c(
      "pooled_mse_uncensored",
      "pooled_mse_impute",
      "pooled_mse_augmented",
      "pooled_mse_trunc_mi"
      # "pooled_mse_complete_case"
    ))
  )


library(looplot)



plot_data = nested_loop_base_data(
    process_df, 
    x = "lod_quantile", steps = c("n"),
    # x = "lod_quantile", steps = c("n","group"),
    grid_cols = "exposure_dist", grid_rows = "h_dist",
    spu_x_shift = .2
)

p = nested_loop_base_plot(
    plot_data,
    x_name = "Detection Limit Quantile",
    y_name = "RMSE", 
    colors = scales::viridis_pal(end = .85, option = "A"),
    grid_scales = "free_y"
)

plot_data = nested_loop_paramsteps_data(
    plot_data,
    steps_y_base = -.2,
    steps_y_height = .2
)

p = nested_loop_paramsteps_plot(
    p, plot_data, 
    steps_values_annotate = TRUE, 
    steps_annotation_size = 5
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
            axis.text.x = element_text(angle = 0, 
                                       vjust = 0.5, 
                                       size = 8)
        ), 
        # add horizontal lines
        add_abline = list(
            intercept = 0
        )
    )
)
print(p)
ggsave("RMSE.png",width = 10,height = 10)

##############################

# R
sens_plot_df <- combined_results$sensspec_summary_long |>
  select(
    n,
    lod_quantile,
    exposure_dist,
    h_dist,
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


sens_plot_df <- sens_plot_df |> filter(metric=="specificity") |> select(
  # n,
  lod_quantile,
  exposure_dist,
  h_dist,
  # metric,
  threshold,
  uncensored,
  impute,
  augmented,
  trunc_mi
) 

plot_data = nested_loop_base_data(
    sens_plot_df, 
    x = "lod_quantile", steps = c("threshold"),
    grid_cols = "exposure_dist", grid_rows = "h_dist",
    spu_x_shift = .2
)

p = nested_loop_base_plot(
    plot_data,
    x_name = "Detection Limit Quantile",
    y_name = "Specificity", 
    colors = scales::viridis_pal(end = .85, option = "A")
)

plot_data = nested_loop_paramsteps_data(
    plot_data,
    steps_y_base = -.2,
    steps_y_height = .2
)

p = nested_loop_paramsteps_plot(
    p, plot_data, 
    steps_values_annotate = TRUE, 
    steps_annotation_size = 5
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
            axis.text.x = element_text(angle = 0, 
                                       vjust = 0.5, 
                                       size = 8)
        )
        # add horizontal lines
        # add_abline = list(
        #     intercept = 0
        # )
    )
)
print(p)
ggsave("Specificity.png",width = 10,height = 10)


