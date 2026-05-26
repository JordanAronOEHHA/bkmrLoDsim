library(looplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)


combined_results <- readRDS("combined_results.rds")

############################## RMSE ##################################

process_df <- combined_results$mse_by_first2_lod_summary |>
  filter(n==800)|>
  filter(group == '++' | group == '+-')|>
  select(
    # n,
    lod_quantile,
    exposure_dist,
    mean_offset,
    scale,
    group,
    h_func,
    any_of(c(
      "pooled_mse_uncensored",
      "pooled_mse_impute",
      "pooled_mse_augmented",
      "pooled_mse_trunc_mi"
      # "pooled_mse_complete_case"
    ))
  )


colnames(process_df)[(ncol(process_df) - 3):ncol(process_df)] <- c("Oracle", "Single Imputation", "Augmented", "Truncated MI")

process_df$h_func <- factor(process_df$h_func, 
                      levels = c("2", "3"), 
                      labels = c("No Interaction", "Interaction"))

# plot_data = nested_loop_base_data(
#     process_df, 
#     x = "n", steps = c("group","exposure_dist","scale"),
#     grid_cols = "mean_offset", grid_rows = "lod_quantile",
#     spu_x_shift = .2 * 1000
# )


plot_data = nested_loop_base_data(
    process_df, 
    x = "lod_quantile", steps = c("mean_offset","exposure_dist","scale"),
    grid_cols = "h_func", grid_rows = "group",
    spu_x_shift = .25
)

plot_data = nested_loop_paramsteps_data(
    plot_data,
    steps_y_base = -.025,
    steps_y_height = .025,
    steps_y_shift = .075
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
        # adjust_ylim = list(
        #     y_expand_add = c(.25, NULL)
        # ),
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
ggsave("RMSE.png",width = 10,height = 10)

##############################

############################## Sensitivity/Specificity ##################################

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


################################################ Specificity ################################################

specificity_df <- sens_plot_df |> 
  # filter(metric=="sensitivity") |> 
  filter(metric=="specificity") |> 
  filter(n==800) |>
  select(
  # n,
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

colnames(specificity_df)[(length(colnames(specificity_df))-3):length(colnames(specificity_df))] <- c("Oracle", "Single Imputation","Augmented or", "Truncated MI")


specificity_df$h_func <- factor(specificity_df$h_func, 
                      levels = c("2", "3"), 
                      labels = c("No Interaction", "Interaction"))

plot_data = nested_loop_base_data(
    specificity_df, 
    x = "lod_quantile", steps = c("mean_offset","exposure_dist","scale"),
    grid_rows = "threshold", grid_cols = "h_func",
    spu_x_shift = .3
)


p = nested_loop_base_plot(
    plot_data,
    x_name = "Detection Limit Quantile",
    y_name = "Specificity", 
    colors = scales::viridis_pal(end = .85, option = "A"),
    grid_scales = "free_y"
)

plot_data = nested_loop_paramsteps_data(
    plot_data,
    steps_y_base = 1.175,
    steps_y_height = .025,
    steps_y_shift = .025
)

p = nested_loop_paramsteps_plot(
    p, plot_data, 
    steps_values_annotate = TRUE, 
    steps_annotation_size = 3
)

p = add_processing(
    p, 
    list(
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
ggsave("Specificity.png",width = 9,height = 9)


################################################ Sensitivity ################################################

sensitivity_df <- sens_plot_df |> 
  filter(metric=="sensitivity") |> 
  filter(n==800) |>
  select(
  # n,
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

colnames(sensitivity_df)[(length(colnames(sensitivity_df))-3):length(colnames(sensitivity_df))] <- c("Oracle", "Single Imputation","Augmented or", "Truncated MI")


sensitivity_df$h_func <- factor(sensitivity_df$h_func, 
                      levels = c("2", "3"), 
                      labels = c("No Interaction", "Interaction"))

plot_data = nested_loop_base_data(
    sensitivity_df, 
    x = "lod_quantile", steps = c("mean_offset","exposure_dist","scale"),
    grid_rows = "threshold", grid_cols = "h_func",
    spu_x_shift = .3
)


p = nested_loop_base_plot(
    plot_data,
    x_name = "Detection Limit Quantile",
    y_name = "Sensitivity", 
    colors = scales::viridis_pal(end = .85, option = "A"),
    grid_scales = "free_y"
)

plot_data = nested_loop_paramsteps_data(
    plot_data,
    steps_y_base = 1.175,
    steps_y_height = .025,
    steps_y_shift = .025
)

p = nested_loop_paramsteps_plot(
    p, plot_data, 
    steps_values_annotate = TRUE, 
    steps_annotation_size = 3
)

p = add_processing(
    p, 
    list(
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
ggsave("Sensitivity.png",width = 9,height = 9)


############################## Coverage ##################################


coverage_df <- combined_results$coverage_by_first2_lod_summary |>
  filter(group == '++' | group == '+-' ) |>
  select(
    n,
    lod_quantile,
    exposure_dist,
    group,
    mean_offset,
    scale,
    h_func,
    method,
    empirical_coverage,
    mean_ci_width
  ) |>
  pivot_longer(
    cols = c(empirical_coverage, mean_ci_width),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = recode(
      metric,
      empirical_coverage = "Coverage",
      mean_ci_width = "Width"
    )
  ) |>
  pivot_wider(
    names_from = method,
    values_from = value
  )

coverage_df <- coverage_df |> 
  filter(metric=="Coverage") |>
  filter(n==800) |>
  select(-'metric',-'n',-'complete_case') |>
  relocate(uncensored, impute, augmented, trunc_mi, .after = h_func )

col_len <- ncol(coverage_df)
colnames(coverage_df)[(col_len-3):col_len] <- c("Oracle", "Single Imputation","Augmented", "Truncated MI")

coverage_df$h_func <- factor(coverage_df$h_func, 
                      levels = c("2", "3"), 
                      labels = c("No Interaction", "Interaction"))


plot_data = nested_loop_base_data(
    coverage_df, 
    x = "lod_quantile", steps = c("mean_offset","exposure_dist","scale"),
    grid_cols = "h_func", grid_rows = "group",
    spu_x_shift = .3
)


p = nested_loop_base_plot(
    plot_data,
    x_name = "Detection Limit Quantile",
    y_name = "CI Coverage", 
    colors = scales::viridis_pal(end = .85, option = "A"),
    grid_scales = "free_y"
)

plot_data = nested_loop_paramsteps_data(
    plot_data,
    steps_y_base = 1.3,
    steps_y_height = .025,
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
        # adjust theme
        add_custom_theme = list(
            axis.text.x = element_text(angle = 90, 
                                       vjust = 0.5, 
                                       size = 6)
          ),
        add_abline = list(
            intercept = 0.95
        )
      )
)
print(p)
ggsave("Coverage.png",width = 10,height = 10)




############################## Contrasts ##################################


contrast_df <- combined_results$contrast_summary |>
  select(
    n,
    lod_quantile,
    exposure_dist,
    mean_offset,
    scale,
    h_func,
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
  filter(n==800) |>
  select(-'n',-'complete_case') |>
  relocate(uncensored, impute, augmented, trunc_mi, .after = contrast_family )

col_len <- ncol(contrast_df)
colnames(contrast_df)[(col_len-3):col_len] <- c("Oracle", "Single Imputation","Augmented", "Truncated MI")

contrast_df$h_func <- factor(contrast_df$h_func, 
                      levels = c("2", "3"), 
                      labels = c("No Interaction", "Interaction"))


plot_data = nested_loop_base_data(
    contrast_df, 
    x = "lod_quantile", steps = c("mean_offset","exposure_dist","scale"),
    grid_cols = "h_func", grid_rows = "contrast_family",
    spu_x_shift = .3
)


p = nested_loop_base_plot(
    plot_data,
    x_name = "Detection Limit Quantile",
    y_name = "Contrast RMSE", 
    colors = scales::viridis_pal(end = .85, option = "A"),
    grid_scales = "free_y"
)

plot_data = nested_loop_paramsteps_data(
    plot_data,
    steps_y_base = -.05,
    steps_y_height = .025,
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
        # adjust theme
        add_custom_theme = list(
            axis.text.x = element_text(angle = 90, 
                                       vjust = 0.5, 
                                       size = 6)
          ),
        add_abline = list(
            intercept = 0.95
        )
      )
)
print(p)
ggsave("Contrast.png",width = 10,height = 10)
