#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param cohort_df
make_matched_mdr_predictions <- function(cohort_df) {

  # Positive or negative microscopy prediction-------------
  df <- cohort_df$cxr_microscopy_matched_df %>% ungroup() %>%
    mutate(across(where(is.logical), ~as.integer(.x)))
  df %<>% filter(type_of_resistance == "MDR non XDR")
  df %<>% select(overallpercentofabnormalvolume:aremediastinallymphnodespresent,
                 matches("_sextant"), microscopyresults, timika_score:nodule,
                 starts_with("qure")) %>%
    select(-ends_with("_num"))
  df %<>% mutate(microscopyresults = factor(ifelse(microscopyresults == "Negative", "neg", "pos"), levels = c("neg", "pos")))
  df %<>% filter()
  
  set.seed(-10); cv_folds <- vfold_cv(df, strata = "microscopyresults")
  
  compute_roc <- function(split, wflow, out_col) {
    
    out_col <- ensym(out_col)
    
    # Fit the model to 90% of the data
    mod <- fit(wflow, data = analysis(split))
    # Predict the other 10%
    pred <- predict(mod, new_data = assessment(split), type = "prob")
    # Compute the area under the ROC curve
    pred %>% 
      bind_cols(assessment(split)) %>% 
      roc_auc(!!out_col, .pred_neg) %>% 
      pull(.estimate)
  }
  
  logit_mod <-
    logistic_reg() %>%
    set_engine("glm")
  
  rand_mod <- 
    rand_forest(trees = 2000, mode = "classification", ) %>%
    set_engine("randomForest")
  
  all_recipe <- recipe(microscopyresults ~ ., data = df) %>%
    step_rm(timika_score | starts_with("qure")) %>%
    step_upsample()
  
  all_corr_recipe <- all_recipe %>%
    step_nzv(all_predictors()) %>%
    step_select_mrmr(all_predictors(), top_p = 10, outcome = "microscopyresults")
  
  timika_recipe <- recipe(microscopyresults ~ timika_score, data = df) %>%
    step_upsample()
  
  timika_plus_recipe <- recipe(microscopyresults ~ timika_score + upper + middle + lower + nodule, data = df) %>%
    step_upsample()
  
  glm_workflow <- 
    workflow() %>%
    add_model(logit_mod)
  
  rand_workflow <- 
    workflow() %>%
    add_model(rand_mod)
  
  all_glm <- glm_workflow %>%
    add_recipe(all_corr_recipe)
  
  all_rand <- rand_workflow %>%
    add_recipe(all_recipe)
  
  timika_glm <- glm_workflow %>%
    add_recipe(timika_recipe)
  
  timika_rand <- rand_workflow %>%
    add_recipe(timika_recipe)
  
  timika_plus_glm <- glm_workflow %>%
    add_recipe(timika_plus_recipe)
  
  timika_plus_rand <- rand_workflow %>%
    add_recipe(timika_plus_recipe)
  
  # Get model interpretability information
  all_rand_fit <- all_rand %>% fit(object = ., data = df %>% mutate(overallpercentofabnormalvolume = as.double(overallpercentofabnormalvolume)))
  e_microscopy <- explain_tidymodels(model = all_rand_fit, data = df %>% mutate(overallpercentofabnormalvolume = as.double(overallpercentofabnormalvolume)), y = df$microscopyresults)
  microscopy_pdp_overall <- model_profile(
    e_microscopy,
    variables = "overallpercentofabnormalvolume",
    N = NULL,
    groups = "cavity"
  )
  microscopy_pdp_overall_plot <- as_tibble(microscopy_pdp_overall$agr_profiles) %>%
    mutate(`_label_` = str_remove(`_label_`, "workflow_")) %>%
    ggplot(aes(`_x_`, `_yhat_`, color = `_label_`)) +
    geom_line(size = 1.2, alpha = 0.8) +
    labs(
      x = "Overall percent of abnormal lung volume",
      y = "Predicted probability of positive baseline microscopy",
      color = NULL,
      title = "Partial dependence plot for baseline microscopy",
      subtitle = "Predictions from a random forest model"
    ) + theme(legend.position = "none")
  microscopy_imp_vars <- all_rand_fit %>%
    extract_fit_parsnip() %>%
    vip::vip(num_features = 10)
  
  model_interpretability_microscopy <- list("all_rand_fit" = all_rand_fit, "pdp_overallabrnormal_cavity" = microscopy_pdp_overall_plot, "imp_vars" = microscopy_imp_vars)
  
  roc_values <- 
    cv_folds %>% 
    mutate(
      all_log = map_dbl(splits, ~compute_roc(split = .x, out_col = "microscopyresults", wflow = all_glm)),
      timika_log      = map_dbl(splits, ~compute_roc(split = .x, out_col = "microscopyresults", wflow = timika_glm)),
      timika_plus_log = map_dbl(splits, ~compute_roc(split = .x, out_col = "microscopyresults", wflow = timika_plus_glm)),
      all_rand = map_dbl(splits, ~compute_roc(split = .x, out_col = "microscopyresults", wflow = all_rand)),
      timika_rand = map_dbl(splits, ~compute_roc(split = .x, out_col = "microscopyresults", wflow = timika_rand)),
      timika_plus_rand = map_dbl(splits, ~compute_roc(split = .x, out_col = "microscopyresults", wflow = timika_plus_rand))
    )
  
  rset_mod <- perf_mod(roc_values, seed = 10, iter = 10000, chains = 5, refresh = 1, hetero_var = T)
  
  preproc_diff <- contrast_models(rset_mod, seed = -10)
  
  roc_longer <- 
    roc_values %>% 
    select(-splits) %>% 
    pivot_longer(cols = c(-id), names_to = "preprocessor", values_to = "roc")
  
  roc_fit <- glm(roc ~ preprocessor, roc_longer, family = gaussian)
  
  plot_roc_fit <- roc_fit %>% 
    augment() %>% 
    ggplot(aes(sample = .resid)) + 
    geom_qq() + 
    geom_qq_line(lty = 2) +
    coord_fixed(ratio  = 20)
  
  # Positive or negative outcome prediction-------------
  df <- cohort_df$cxr_microscopy_matched_df %>% ungroup() %>%
    mutate(across(where(is.logical), ~as.integer(.x)))
  df %<>% filter(type_of_resistance == "MDR non XDR")
  df %<>% select(overallpercentofabnormalvolume:aremediastinallymphnodespresent, matches("_sextant"), outcome, timika_score:nodule, starts_with("qure")) %>%
    select(-ends_with("_num"))
  df %<>% mutate(outcome = factor(ifelse(outcome %in% c("Completed", "Cured"), "neg", "pos"), levels = c("neg", "pos")))
  set.seed(-10); cv_folds <- vfold_cv(df, strata = "outcome")
  
  logit_mod <-
    logistic_reg() %>%
    set_engine("glm")
  
  rand_mod <- 
    rand_forest(trees = 2000, mode = "classification") %>%
    set_engine("randomForest")
  
  all_recipe <- recipe(outcome ~ ., data = df) %>%
    step_rm(timika_score | starts_with("qure")) %>%
    step_upsample()
  
  all_corr_recipe <- all_recipe %>%
    step_nzv(all_predictors()) %>%
    step_select_mrmr(all_predictors(), top_p = 10, outcome = "outcome")
  
  timika_recipe <- recipe(outcome ~ timika_score, data = df) %>%
    step_upsample()
  
  timika_plus_recipe <- recipe(outcome ~ timika_score + upper + middle + lower + nodule, data = df) %>%
    step_upsample()
  
  glm_workflow <- 
    workflow() %>%
    add_model(logit_mod)
  
  rand_workflow <- 
    workflow() %>%
    add_model(rand_mod)
  
  all_glm <- glm_workflow %>%
    add_recipe(all_corr_recipe)
  
  all_rand <- rand_workflow %>%
    add_recipe(all_recipe)
  
  timika_glm <- glm_workflow %>%
    add_recipe(timika_recipe)
  
  timika_rand <- rand_workflow %>%
    add_recipe(timika_recipe)
  
  timika_plus_glm <- glm_workflow %>%
    add_recipe(timika_plus_recipe)
  
  timika_plus_rand <- rand_workflow %>%
    add_recipe(timika_plus_recipe)
  
  # Get model interpretability information
  all_rand_fit <- all_rand %>% fit(object = ., data = df %>% mutate(overallpercentofabnormalvolume = as.double(overallpercentofabnormalvolume)))
  e_treatment <- explain_tidymodels(model = all_rand_fit, data = df %>% mutate(overallpercentofabnormalvolume = as.double(overallpercentofabnormalvolume)), y = df$outcome)
  treatment_pdp_overall <- model_profile(
    e_treatment,
    variables = "overallpercentofabnormalvolume",
    N = NULL,
    groups = "cavity"
  )
  treatment_pdp_overall_plot <- as_tibble(treatment_pdp_overall$agr_profiles) %>%
    mutate(`_label_` = str_remove(`_label_`, "workflow_")) %>%
    ggplot(aes(`_x_`, `_yhat_`, color = `_label_`)) +
    geom_line(size = 1.2, alpha = 0.8) +
    labs(
      x = "Overall percent of abnormal lung volume",
      y = "Predicted probability of treatment failure",
      color = NULL,
      title = "Partial dependence plot for treatment failure",
      subtitle = "Predictions from a random forest model"
    ) + theme(legend.position = "none")
  treatment_imp_vars <- all_rand_fit %>%
    extract_fit_parsnip() %>%
    vip::vip(num_features = 10)
  
  model_interpretability_treatment <- list("all_rand_fit" = all_rand_fit, "pdp_overallabrnormal_cavity" = treatment_pdp_overall_plot, "imp_vars" = treatment_imp_vars)
  
  roc_values_outcome <- 
    cv_folds %>% 
    mutate(
      all_log_outcome = map_dbl(splits, ~compute_roc(split = .x, out_col = "outcome", wflow = all_glm)),
      timika_log_outcome      = map_dbl(splits, ~compute_roc(split = .x, out_col = "outcome", wflow = timika_glm)),
      timika_plus_log_outcome = map_dbl(splits, ~compute_roc(split = .x, out_col = "outcome", wflow = timika_plus_glm)),
      all_rand_outcome = map_dbl(splits, ~compute_roc(split = .x, out_col = "outcome", wflow = all_rand)),
      timika_rand_outcome = map_dbl(splits, ~compute_roc(split = .x, out_col = "outcome", wflow = timika_rand)),
      timika_plus_rand_outcome = map_dbl(splits, ~compute_roc(split = .x, out_col = "outcome", wflow = timika_plus_rand))
    )
  
  rset_mod_outcome <- perf_mod(roc_values_outcome, seed = 10, iter = 10000, chains = 5, refresh = 1, hetero_var = T)
  
  preproc_diff_outcome <- contrast_models(rset_mod_outcome, seed = -10)
  
  roc_longer <- 
    roc_values_outcome %>% 
    select(-splits) %>% 
    pivot_longer(cols = c(-id), names_to = "preprocessor", values_to = "roc")
  
  roc_fit_outcome <- glm(roc ~ preprocessor, roc_longer, family = gaussian)
  
  plot_roc_fit_outcome <- roc_fit_outcome %>% 
    augment() %>% 
    ggplot(aes(sample = .resid)) + 
    geom_qq() + 
    geom_qq_line(lty = 2) +
    coord_fixed(ratio  = 20)
  
  #### Comparison of both ------------
  
  roc_values_both <- roc_values %>% left_join(roc_values_outcome %>% select(-splits))
  
  rset_mod_both <- perf_mod(roc_values_both, seed = 10, iter = 10000, chains = 5, refresh = 1, hetero_var = T)
  
  preproc_diff_both <- contrast_models(rset_mod_both, seed = -10)
  
  roc_longer <- 
    roc_values_both %>% 
    select(-splits) %>% 
    pivot_longer(cols = c(-id), names_to = "preprocessor", values_to = "roc")
  
  roc_fit_both <- glm(roc ~ preprocessor, roc_longer, family = gaussian)
  
  plot_roc_fit_both <- roc_fit_outcome %>% 
    augment() %>% 
    ggplot(aes(sample = .resid)) + 
    geom_qq() + 
    geom_qq_line(lty = 2) +
    coord_fixed(ratio  = 20)
  
  
  ### Return results -----------
  
  list("rset_mod" = rset_mod,
       "rset_mod_outcome" = rset_mod_outcome,
       "rset_mod_both" = rset_mod_both,
       "preproc_diff" = preproc_diff,
       "preproc_diff_outcome" = preproc_diff_outcome,
       "preproc_diff_both" = preproc_diff_both,
       "roc_fit" = roc_fit,
       "roc_fit_outcome" = roc_fit_outcome,
       "roc_fit_both" = roc_fit_both,
       "roc_values" = roc_values,
       "roc_values_outcome" = roc_values_outcome,
       "roc_values_both" = roc_values_both,
       "model_interpretability_microscopy" = model_interpretability_microscopy,
       "model_interpretability_treatment" = model_interpretability_treatment)
  

}
