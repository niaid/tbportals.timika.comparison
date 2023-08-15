#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param sensitive_predictions
#' @param mdr_predictions
#' @param matched_sensitive_predictions
#' @param matched_mdr_predictions
#' @param best_comparisons
#' @param best_matched_comparisons
#' @param cohort_df
make_visualizations_tables <- function(sensitive_predictions, mdr_predictions,
                                       matched_sensitive_predictions,
                                       matched_mdr_predictions,
                                       best_comparisons,
                                       best_matched_comparisons, cohort_df) {

  raw <- c("sensitive_predictions", "mdr_predictions",
           "matched_sensitive_predictions",
           "matched_mdr_predictions",
           "best_comparisons",
           "best_matched_comparisons", "cohort_df")
  
  #Cohort description
  smd_plot <- plot(summary(cohort_df$matched_summary))
  
  unmatched_cohort_desc <- tableby(type_of_resistance ~ microscopyresults + outcome + timika_score + upper + middle + lower + cavity +
                                     nodule + registration_date + age_of_onset + gender + country + education + employment + case_definition +
                                     bmi + hiv + risk_smoker + risk_alcohol + risk_drug, data =  cohort_df$cxr_microscopy_outcome_df %>% slice_head() %>% ungroup()) %>%
    summary() %>%
    data.frame() %>%
    mutate_all(~gsub("&nbsp;", "", .x))
  
  matched_cohort_desc <- tableby(type_of_resistance ~ microscopyresults + outcome + timika_score + upper + middle + lower + cavity +
                                     nodule + registration_date + age_of_onset + gender + country + education + employment + case_definition +
                                     bmi + hiv + risk_smoker + risk_alcohol + risk_drug, data =  cohort_df$cxr_microscopy_matched_df %>% slice_head() %>% ungroup()) %>%
    summary() %>%
    data.frame() %>%
    mutate_all(~gsub("&nbsp;", "", .x))
  
  # Cross-validation summary tables
  sensitive_cv <- sensitive_predictions$roc_values_both %>%
    select(-c(1:2)) %>%
    pivot_longer(everything(), names_to = "model", values_to = "roc") %>%
    group_by(model) %>% 
    summarise(across(.fns = list("mean" = ~mean(.x), "median" = ~median(.x), "sd" = ~sd(.x), "mad" = ~mad(.x)))) %>%
    mutate(resistance = "sensitive",
           matching = "unmatched")
  
  mdr_cv <- mdr_predictions$roc_values_both %>%
    select(-c(1:2)) %>%
    pivot_longer(everything(), names_to = "model", values_to = "roc") %>%
    group_by(model) %>% 
    summarise(across(.fns = list("mean" = ~mean(.x), "median" = ~median(.x), "sd" = ~sd(.x), "mad" = ~mad(.x)))) %>%
    mutate(resistance = "mdr",
           matching = "unmatched")
  
  sensitive_matched_cv <- matched_sensitive_predictions$roc_values_both %>%
    select(-c(1:2)) %>%
    pivot_longer(everything(), names_to = "model", values_to = "roc") %>%
    group_by(model) %>% 
    summarise(across(.fns = list("mean" = ~mean(.x), "median" = ~median(.x), "sd" = ~sd(.x), "mad" = ~mad(.x)))) %>%
    mutate(resistance = "sensitive",
           matching = "matched")
  
  mdr_matched_cv <- matched_mdr_predictions$roc_values_both %>%
    select(-c(1:2)) %>%
    pivot_longer(everything(), names_to = "model", values_to = "roc") %>%
    group_by(model) %>% 
    summarise(across(.fns = list("mean" = ~mean(.x), "median" = ~median(.x), "sd" = ~sd(.x), "mad" = ~mad(.x)))) %>%
    mutate(resistance = "mdr",
           matching = "matched")
  
  cv_summary <- bind_rows(sensitive_cv, mdr_cv, sensitive_matched_cv, mdr_matched_cv)
  cv_summary %<>% mutate(outcome = ifelse(grepl("outcome", model), "treatment", "baseline_sputum"))
  
  rm(sensitive_cv, mdr_cv, sensitive_matched_cv, mdr_matched_cv)
  
  # Bayesian Anova model performance
  sen_anova_model_metrics = sensitive_predictions$rset_mod_both %>% summary()
  sen_mcmc_traceplot <- mcmc_trace(x = sensitive_predictions$rset_mod_both$stan)
  sen_mcmc_density <- mcmc_dens_overlay(x = sensitive_predictions$rset_mod_both$stan)
  sen_rhat <- mcmc_rhat(rhat(sensitive_predictions$rset_mod_both$stan))
  sen_neff <- mcmc_neff(neff_ratio(object = sensitive_predictions$rset_mod_both$stan))
  
  mdr_anova_model_metrics = mdr_predictions$rset_mod_both %>% summary()
  mdr_mcmc_traceplot <- mcmc_trace(x = mdr_predictions$rset_mod_both$stan)
  mdr_mcmc_density <- mcmc_dens_overlay(x = mdr_predictions$rset_mod_both$stan)
  mdr_rhat <- mcmc_rhat(rhat(mdr_predictions$rset_mod_both$stan))
  mdr_neff <- mcmc_neff(neff_ratio(object = mdr_predictions$rset_mod_both$stan))
  
  matched_sen_anova_model_metrics = matched_sensitive_predictions$rset_mod_both %>% summary()
  matched_sen_mcmc_traceplot <- mcmc_trace(x = matched_sensitive_predictions$rset_mod_both$stan)
  matched_sen_mcmc_density <- mcmc_dens_overlay(x = matched_sensitive_predictions$rset_mod_both$stan)
  matched_sen_rhat <- mcmc_rhat(rhat(matched_sensitive_predictions$rset_mod_both$stan))
  matched_sen_neff <- mcmc_neff(neff_ratio(object = matched_sensitive_predictions$rset_mod_both$stan))
  
  matched_mdr_anova_model_metrics = matched_mdr_predictions$rset_mod_both %>% summary()
  matched_mdr_mcmc_traceplot <- mcmc_trace(x = matched_mdr_predictions$rset_mod_both$stan)
  matched_mdr_mcmc_density <- mcmc_dens_overlay(x = matched_mdr_predictions$rset_mod_both$stan)
  matched_mdr_rhat <- mcmc_rhat(rhat(matched_mdr_predictions$rset_mod_both$stan))
  matched_mdr_neff <- mcmc_neff(neff_ratio(object = matched_mdr_predictions$rset_mod_both$stan))
  
  # Bayesian Anova Results
  sen_microscopy <- autoplot(sensitive_predictions$rset_mod) + ggtitle("Sensitive baseline sputum microscopy + (unmatched)")
  mdr_microscopy <- autoplot(mdr_predictions$rset_mod) + ggtitle("MDR baseline sputum microscopy + (unmatched)")
  sen_outcome <- autoplot(sensitive_predictions$rset_mod_outcome) + ggtitle("Sensitive treatment outcome (unmatched)")
  mdr_outcome <- autoplot(mdr_predictions$rset_mod_outcome) + ggtitle("MDR treatment outcome (unmatched)")

  matched_sen_microscopy <- autoplot(matched_sensitive_predictions$rset_mod) + ggtitle("Sensitive baseline sputum microscopy + (matched)")
  matched_mdr_microscopy <- autoplot(matched_mdr_predictions$rset_mod) + ggtitle("MDR baseline sputum microscopy + (matched)")
  matched_sen_outcome <- autoplot(matched_sensitive_predictions$rset_mod_outcome) + ggtitle("Sensitive treatment outcome (matched)")
  matched_mdr_outcome <- autoplot(matched_mdr_predictions$rset_mod_outcome) + ggtitle("MDR treatment outcome (matched)")
  
  all_best <- bind_rows(best_comparisons, best_matched_comparisons)
  all_best %<>%
    mutate(n = str_count(pattern = "outcome", string = contrast)) %>%
    mutate(outcome = case_when(n == 0 ~ "baseline_sputum",
                               n == 1 ~ "timika outcome comparison",
                               n == 2 ~ "treatment")) %>%
    select(-n)
  
  sen <- summary(sensitive_predictions$preproc_diff_both, size = 0.02) %>%
    mutate(resistance = "sensitive", matching = "unmatched")
  
  sen <- bind_rows(sen, summary(matched_sensitive_predictions$preproc_diff_both, size = 0.02) %>%
                     mutate(resistance = "sensitive", matching = "matched"))
  
  mdr <- summary(mdr_predictions$preproc_diff_both, size = 0.02) %>%
    mutate(resistance = "mdr", matching = "unmatched")
  
  mdr <- bind_rows(mdr, summary(matched_mdr_predictions$preproc_diff_both, size = 0.02) %>%
                     mutate(resistance = "mdr", matching = "matched"))
  
  all_comparisons <- bind_rows(sen, mdr)
  
  rm(sen, mdr)
  
  l <- ls()
  l <- setdiff(l, c(raw, "raw", "l"))
  
  names(l) <- l
  
  l <- map(l, ~get0(.x))
  
  l
}
