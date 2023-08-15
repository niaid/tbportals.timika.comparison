#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param matched_sensitive_predictions
#' @param matched_mdr_predictions
make_best_matched_comparisons <- function(matched_sensitive_predictions,
                                          matched_mdr_predictions) {

  sen_summary <- matched_sensitive_predictions$rset_mod_both %>% tidy() %>% group_by(model) %>% summarise(mean_posterior = mean(posterior))
  best_timika_sen <- sen_summary %>%
    filter(model %in% c("timika_log", "timika_rand")) %>%
    filter(mean_posterior == max(mean_posterior))
  best_timika_sen <- best_timika_sen$model %>% unlist() %>% paste("\\b", ., "\\b", sep = "")
  best_model_sen <- sen_summary %>%
    filter(!(model %in% c("timika_log", "timika_rand")) & !grepl("outcome", model)) %>%
    filter(mean_posterior == max(mean_posterior))
  best_model_sen <- best_model_sen$model %>% unlist() %>% paste("\\b", ., "\\b", sep = "")
  
  mdr_summary <- matched_mdr_predictions$rset_mod_both %>% tidy() %>% group_by(model) %>% summarise(mean_posterior = mean(posterior))
  best_timika_mdr <- mdr_summary %>%
    filter(model %in% c("timika_log", "timika_rand")) %>%
    filter(mean_posterior == max(mean_posterior))
  best_timika_mdr <- best_timika_mdr$model %>% unlist() %>% paste("\\b", ., "\\b", sep = "")
  best_model_mdr <- mdr_summary %>%
    filter(!(model %in% c("timika_log", "timika_rand")) & !grepl("outcome", model)) %>%
    filter(mean_posterior == max(mean_posterior))
  best_model_mdr <- best_model_mdr$model %>% unlist() %>% paste("\\b", ., "\\b", sep = "")
  
  sen_comparison = matched_sensitive_predictions$preproc_diff_both %>%
    summary(size = 0.02) %>%
    filter(grepl(best_timika_sen, contrast) & grepl(best_model_sen, contrast))
  mdr_comparison = matched_mdr_predictions$preproc_diff_both %>%
    summary(size = 0.02) %>%
    filter(grepl(best_timika_mdr, contrast) & grepl(best_model_mdr, contrast))
  
  best_timika_outcome_sen <- sen_summary %>%
    filter(model %in% c("timika_log_outcome", "timika_rand_outcome")) %>%
    filter(mean_posterior == max(mean_posterior))
  best_timika_outcome_sen <- best_timika_outcome_sen$model %>% unlist() %>% paste("\\b", ., "\\b", sep = "")
  best_model_outcome_sen <- sen_summary %>%
    filter(!(model %in% c("timika_log_outcome", "timika_rand_outcome")) & grepl("outcome", model)) %>%
    filter(mean_posterior == max(mean_posterior))
  best_model_outcome_sen <- best_model_outcome_sen$model %>% unlist() %>% paste("\\b", ., "\\b", sep = "")
  
  best_timika_outcome_mdr <- mdr_summary %>%
    filter(model %in% c("timika_log_outcome", "timika_rand_outcome")) %>%
    filter(mean_posterior == max(mean_posterior))
  best_timika_outcome_mdr <- best_timika_outcome_mdr$model %>% unlist() %>% paste("\\b", ., "\\b", sep = "")
  best_model_outcome_mdr <- mdr_summary %>%
    filter(!(model %in% c("timika_log_outcome", "timika_rand_outcome")) & grepl("outcome", model)) %>%
    filter(mean_posterior == max(mean_posterior))
  best_model_outcome_mdr <- best_model_outcome_mdr$model %>% unlist() %>% paste("\\b", ., "\\b", sep = "")
  
  sen_outcome_comparison = matched_sensitive_predictions$preproc_diff_both %>%
    summary(size = 0.02) %>%
    filter(grepl(best_timika_outcome_sen, contrast) & grepl(best_model_outcome_sen, contrast))
  mdr_outcome_comparison = matched_mdr_predictions$preproc_diff_both %>%
    summary(size = 0.02) %>%
    filter(grepl(best_timika_outcome_mdr, contrast) & grepl(best_model_outcome_mdr, contrast))
  
  sen_both_comparison = matched_sensitive_predictions$preproc_diff_both %>%
    summary(size = 0.02) %>%
    filter(grepl(best_timika_sen, contrast) & grepl(best_timika_outcome_sen, contrast))
  mdr_both_comparison = matched_mdr_predictions$preproc_diff_both %>%
    summary(size = 0.02) %>%
    filter(grepl(best_timika_mdr, contrast) & grepl(best_timika_outcome_mdr, contrast))
  
  sensitive_comparison <- bind_rows(sen_comparison, sen_outcome_comparison, sen_both_comparison) %>%
    mutate(resistance = "sensitive", matching = "matched")
  mdr_comparison <- bind_rows(mdr_comparison, mdr_outcome_comparison, mdr_both_comparison) %>%
    mutate(resistance = "mdr", matching = "matched")
  
  all_comparison = bind_rows(sensitive_comparison, mdr_comparison)
  
  all_comparison

}
