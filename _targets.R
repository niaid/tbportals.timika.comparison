source("./packages.R")

lapply(list.files("./R", full.names = TRUE), source)

tar_option_set(
  format = "qs",
  resources = tar_resources(
    qs = tar_resources_qs(preset = "fast")
  )
)

tar_plan(
  
  tar_target(endpoints, list_endpoints()$endpoint),
  
  tar_target(endpoints_dfs, get_endpoints_dfs(endpoints), pattern = map(endpoints), iteration = "list",
             cue = tar_cue_age(name = endpoints_dfs, age = as.difftime(14, units = "days"))),
  
  tar_target(cohort_df, generate_cohort_df(endpoints_dfs)),
  
  tar_target(sensitive_predictions, make_sensitive_predictions(cohort_df)),
  
  tar_target(mdr_predictions, make_mdr_predictions(cohort_df)),
  
  tar_target(matched_sensitive_predictions, make_matched_sensitive_predictions(cohort_df)),
  
  tar_target(matched_mdr_predictions, make_matched_mdr_predictions(cohort_df)),
  
  tar_target(best_comparisons, make_best_comparisons(sensitive_predictions, mdr_predictions)),
  
  tar_target(best_matched_comparisons, make_best_matched_comparisons(matched_sensitive_predictions, matched_mdr_predictions)),
  
  tar_target(visualizations_tables, make_visualizations_tables(sensitive_predictions, mdr_predictions,
                                                               matched_sensitive_predictions, matched_mdr_predictions,
                                                               best_comparisons, best_matched_comparisons,
                                                               cohort_df)),
  
  tar_target(sensitive_itt_predictions, make_itt_sensitive_predictions(cohort_df)),
  
  tar_target(mdr_itt_predictions, make_itt_mdr_predictions(cohort_df)),
  
  tar_target(best_itt_comparisons, make_best_itt_comparisons(sensitive_itt_predictions, mdr_itt_predictions)),
  
  tar_target(flextable_images, make_flextable_images(cohort_df, visualizations_tables, best_itt_comparisons), format = "file"),
  
  tar_target(tables_editable, make_tables_editable(cohort_df, visualizations_tables, best_itt_comparisons), format = "file"),
  
  tar_quarto(name = article, path = "article.qmd",
             execute_params = list(num_itt = abs(cohort_df$cxr_microscopy_outcome_qure_counts$loss[grepl("Has treatment outcome",
                                                                                                         cohort_df$cxr_microscopy_outcome_qure_counts$step)]),
                                   num_matching_drop = abs(cohort_df$cxr_microscopy_outcome_qure_counts$loss[grepl("matched",
                                                                                                                   cohort_df$cxr_microscopy_outcome_qure_counts$step)]),
                                   flow_table = cohort_df$cxr_microscopy_outcome_qure_counts %>% dplyr::filter(!grepl("qure", step))
             )
  )
  
)
