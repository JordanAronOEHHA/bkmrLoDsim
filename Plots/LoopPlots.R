library(looplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)


combined_results <- readRDS("combined_results.rds")

############################## RMSE ##################################

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
# ggsave("RMSE++.png",width = 12,height = 12)

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
# ggsave("Sensitivity.png",width = 16,height = 12)


############################## RMSE ##################################


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
