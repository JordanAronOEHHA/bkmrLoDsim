library(looplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)


combined_results <- readRDS("combined_results.rds")

h_func_levels <- c("1", "2", "3", "4")
h_func_labels <- c("Single Active", "No Interaction", "Interaction", "Four Active")
label_h_func <- function(x) {
  factor(as.character(x), levels = h_func_levels, labels = h_func_labels)
}

draw_nested_metric_plot <- function(
  data,
  y_name,
  steps_y_base,
  steps_y_shift,
  spu_x_shift,
  axis_text_size,
  abline_intercept = NULL,
  save_path = NULL,
  save_width = 10,
  save_height = 10,
  x = "lod_quantile",
  steps = c("mean_offset", "scale"),
  grid_rows = "group",
  grid_cols = NULL,
  x_name = "Detection Limit Quantile",
  grid_scales = "free_y",
  steps_y_height = 0.025,
  steps_annotation_size = 3
) {
  if (!is.null(grid_cols)) {
    plot_data <- nested_loop_base_data(
      data,
      x = x,
      steps = steps,
      grid_rows = grid_rows,
      grid_cols = grid_cols,
      spu_x_shift = spu_x_shift
    )
  } else {
    plot_data <- nested_loop_base_data(
      data,
      x = x,
      steps = steps,
      grid_rows = grid_rows,
      spu_x_shift = spu_x_shift
    )
  }

  plot_data <- nested_loop_paramsteps_data(
    plot_data,
    steps_y_base = steps_y_base,
    steps_y_height = steps_y_height,
    steps_y_shift = steps_y_shift
  )

  p <- nested_loop_base_plot(
    plot_data,
    x_name = x_name,
    y_name = y_name,
    colors = scales::viridis_pal(end = .85, option = "A"),
    grid_scales = grid_scales
  )

  p <- nested_loop_paramsteps_plot(
    p,
    plot_data,
    steps_values_annotate = TRUE,
    steps_annotation_size = steps_annotation_size
  )

  processing <- list(
    add_custom_theme = list(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        size = axis_text_size
      )
    )
  )

  if (!is.null(abline_intercept)) {
    processing$add_abline <- list(intercept = abline_intercept)
  }

  p <- add_processing(p, processing)
  print(p)

  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = save_width, height = save_height)
  }

  invisible(p)
}

############################## RMSE ##################################

mse_by_active_lod_burden_summary <- combined_results$prediction_diagnostics_by_active_lod_burden_summary |> select(
    n,
    p,
    lod_quantile,
    exposure_dist,
    correlation,
    mean_offset,
    h_func,
    scale,
    mcmc_iter, 
    group,
    method,
    mean_bias
  ) |>
  pivot_wider(
    names_from = method,
    values_from = mean_bias
  )

process_df <- mse_by_active_lod_burden_summary |>
  select(
    # n,
    lod_quantile,
    # exposure_dist,
    # mean_offset,
    # scale,
    group,
    # h_func,
    correlation, 
    uncensored,
    impute,
    augmented,
    trunc_mi
  )


# colnames(process_df)[(ncol(process_df) - 3):ncol(process_df)] <- c("Oracle", "Single Imputation", "Augmented", "Truncated MI")

# process_df$h_func <- label_h_func(process_df$h_func)


draw_nested_metric_plot(
  process_df,
  y_name = "Bias",
  steps = c("correlation"),
  steps_y_base = -.025,
  steps_y_shift = .075,
  spu_x_shift = .25,
  axis_text_size = 7,
  abline_intercept = 0
)
ggsave("Biasp8.png",width = 10,height = 10)

##############################

############################## Sensitivity/Specificity ##################################

sens_plot_df <- combined_results$sensspec_summary_long |>
  select(
    # n,
    lod_quantile,
    # exposure_dist,
    # mean_offset,
    # h_func,
    # scale,
    threshold,
    method,
    mean_sensitivity,
    mean_specificity,
    correlation
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


################################################ Specificity ################################################

specificity_df <- sens_plot_df |> 
  # filter(metric=="sensitivity") |> 
  filter(metric=="specificity") |> 
  # filter(n==800) |>
  select(
  # n,
  lod_quantile,
  correlation,
  # exposure_dist,
  # mean_offset,
  # scale,
  # h_func,
  # metric,
  threshold,
  uncensored,
  impute,
  # augmented_and,
  augmented_or,
  trunc_mi
) 

# colnames(sens_plot_df)[8:12] <- c("Oracle", "Single Imputation", "Augmented and","Augmented or", "Truncated MI")

colnames(specificity_df)[(length(colnames(specificity_df))-3):length(colnames(specificity_df))] <- c("Oracle", "Single Imputation","Augmented or", "Truncated MI")


specificity_df$h_func <- label_h_func(specificity_df$h_func)

draw_nested_metric_plot(
  specificity_df,
  y_name = "Specificity",
  steps = c("correlation"),
  steps_y_base = 1.175,
  steps_y_shift = .025,
  spu_x_shift = .3,
  axis_text_size = 6,
  grid_rows = "threshold",
  save_path = "Specificity.png",
  save_width = 9,
  save_height = 9
)


################################################ Sensitivity ################################################

sensitivity_df <- sens_plot_df |> 
  filter(metric=="sensitivity") |> 
  # filter(n==800) |>
  select(
  # n,
  lod_quantile,
  correlation,
  # exposure_dist,
  # mean_offset,
  # scale,
  # h_func,
  # metric,
  threshold,
  uncensored,
  impute,
  # augmented_and,
  augmented_or,
  trunc_mi
) 

# colnames(sens_plot_df)[8:12] <- c("Oracle", "Single Imputation", "Augmented and","Augmented or", "Truncated MI")

colnames(sensitivity_df)[(length(colnames(sensitivity_df))-3):length(colnames(sensitivity_df))] <- c("Oracle", "Single Imputation","Augmented or", "Truncated MI")


sensitivity_df$h_func <- label_h_func(sensitivity_df$h_func)

draw_nested_metric_plot(
  sensitivity_df,
  y_name = "Sensitivity",
  # steps = c("mean_offset", "exposure_dist", "scale"),
  steps = c("correlation"),
  grid_rows = "threshold",
  grid_cols = "h_func",
  steps_y_base = 1.175,
  steps_y_shift = .025,
  spu_x_shift = .3,
  axis_text_size = 6,
  save_path = "Sensitivity.png",
  save_width = 9,
  save_height = 9
)


############################## Coverage ##################################

coverage_df <- combined_results$prediction_diagnostics_by_active_lod_burden_summary |> select(
    # n,
    # p,
    lod_quantile,
    # exposure_dist,
    correlation,
    # mean_offset,
    # h_func,
    # scale,
    # mcmc_iter, 
    group,
    method,
    empirical_coverage
  ) |>
  pivot_wider(
    names_from = method,
    values_from = empirical_coverage
  ) |> 
  relocate(uncensored, impute, augmented, trunc_mi, .after = group)
  

# coverage_df <- coverage_df |> 
#   # filter(metric=="Coverage") |>
#   # filter(n==800) |>
#   # select(lod_quantile, group, mean_offset, scale, h_func, uncensored, impute, augmented, trunc_mi) |>
#   relocate(uncensored, impute, augmented, trunc_mi, .after = h_func )

# col_len <- ncol(coverage_df)
# colnames(coverage_df)[(col_len-3):col_len] <- c("Oracle", "Single Imputation","Augmented", "Truncated MI")

coverage_df$h_func <- label_h_func(coverage_df$h_func)


draw_nested_metric_plot(
  coverage_df,
  y_name = "CI Coverage",
  steps = c("correlation"),
  steps_y_base = 1.1,
  steps_y_shift = .05,
  spu_x_shift = .3,
  axis_text_size = 6,
  abline_intercept = 0.95
  # save_path = "Coverage.png",
  # save_width = 10,
  # save_height = 10
)




############################## Contrasts ##################################


contrast_df <- combined_results$contrast_summary |>
  select(
    n,
    lod_quantile,
    exposure_dist,
    mean_offset,
    scale,
    h_func,
    correlation,
    contrast_family,
    method,
    rmse,
    coverage,
  ) |>
  pivot_longer(
    cols = c(rmse, coverage),
    names_to = "metric",
    values_to = "value"
  ) |>
  pivot_wider(
    names_from = method,
    values_from = value
  ) |> 
  filter(metric == "rmse") |>
  select(-'metric')

contrast_df <- contrast_df |> 
  select(lod_quantile, 
    # exposure_dist, 
    # mean_offset, 
    # scale, 
    # h_func, 
    correlation,
    contrast_family, 
    uncensored, 
    impute, 
    augmented, 
    trunc_mi) |>
  relocate(uncensored, impute, augmented, trunc_mi, .after = contrast_family )

col_len <- ncol(contrast_df)
colnames(contrast_df)[(col_len-3):col_len] <- c("Oracle", "Single Imputation","Augmented", "Truncated MI")

contrast_df$h_func <- label_h_func(contrast_df$h_func)


draw_nested_metric_plot(
  contrast_df,
  y_name = "Contrast RMSE",
  # steps = c("mean_offset", "exposure_dist", "scale"),
  steps = c("correlation"),
  grid_rows = "contrast_family",
  # grid_cols = "h_func",
  steps_y_base = -.05,
  steps_y_shift = .05,
  spu_x_shift = .3,
  axis_text_size = 6,
  abline_intercept = 0.95,
  save_path = "Contrast.png",
  save_width = 10,
  save_height = 10
)
