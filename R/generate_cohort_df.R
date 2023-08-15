#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param endpoints_dfs
generate_cohort_df <- function(endpoints_dfs) {

  # Baseline sputum microscopy cohort
  case <- endpoints_dfs[[9]]
  
  treat <- endpoints_dfs[[11]]
  
  treat_starts <- treat %>%
    group_by(condition_id) %>%
    summarise(treat_start = min(c_across(ends_with("start"))))
  
  dst <- endpoints_dfs[[7]]
  
  counts_df <- data.frame(
    "step" = "all cases",
    "type" = "inclusion",
    "num_cases" = n_distinct(case$condition_id),
    "loss" = as.integer(NA)
  )
  
  # Updated cohort selection criteria
  cohort_df <- case
  
  cohort_df %<>% mutate(hiv = grepl("HIV", comorbidity, ignore.case = F),
                        risk_smoker = ifelse(grepl("smoker", social_risk_factors), 1, 0),
                        risk_drug = ifelse(grepl("drug", social_risk_factors), 1, 0),
                        risk_alcohol = ifelse(grepl("alcohol", social_risk_factors), 1, 0)) %>%
    filter(age_of_onset >= 18)
  
  counts_df %<>%
    add_row(step = "Age of onset >= 18",
            type = "inclusion",
            num_cases = n_distinct(cohort_df$condition_id)) %>%
    mutate(loss = num_cases - lag(num_cases))
  
  microscopy <- dst %>% filter(!is.na(microscopyresults) & specimen_collection_site == "sputum") %>%
    select(condition_id, specimen_collection_site, specimen_collection_date, microscopyresults) %>%
    distinct() %>%
    mutate(microscopyresults = factor(x = microscopyresults, levels = c("Negative", "1 to 9 in 100 (1-9/100)",
                                                                        "10 to 99 in 100 (1+)", "1 to 9 in 1 (2+)", 
                                                                        "10 to 99 in 1 (3+)", "More than 99 in 1 (4+)"), ordered = T)) %>%
    filter(!is.na(microscopyresults))
  
  # Get last and worst sputum microscopy
  microscopy %<>%
    left_join(treat_starts) %>%
    group_by(condition_id) %>%
    filter(specimen_collection_date <= treat_start) %>%
    slice_max(specimen_collection_date) %>%
    slice_max(microscopyresults)
  
  cohort_df <- cohort_df[cohort_df$condition_id %in% microscopy$condition_id,]
  
  counts_df %<>%
    add_row(step = "Has positive or negative sputum microscopy before treatment start",
            type = "inclusion",
            num_cases = n_distinct(cohort_df$condition_id)) %>%
    mutate(loss = num_cases - lag(num_cases))
  
  cxr_manual <- endpoints_dfs[[5]]

  cxr_manual %<>% mutate(across(everything(), ~as.character(.x))) %>%
    pivot_longer(cols = c(collapse:highdensitycalcifiedtpicallysequella), 
                 names_to = "measure", values_to = "value") %>%
    pivot_wider(names_from = c(sextant, measure), values_from = value) %>%
    rename_with(., ~tolower(gsub(pattern = " ", replacement = "_", .x)), everything()) %>%
    type.convert() %>%
    select(-starts_with("none_")) %>%
    mutate(across(where(is.character) & !matches("_id"), ~as.factor(.x))) %>%
    mutate(across(where(is.factor) & !matches("rater|modality_cd"), ~droplevels(replace_na(data = .x, replace = "No")))) %>%
    mutate(across(where(is.numeric), ~replace_na(data = .x, replace = 0)))
  
  cxr_manual %<>% left_join(treat_starts)
  cxr_manual %<>% rowwise() %>% filter(between(imaging_date, treat_start - 90, treat_start))
  
  cohort_df <- cohort_df[cohort_df$condition_id %in% cxr_manual$condition_id,]
  
  counts_df %<>%
    add_row(step = "Has annotated CXR within 90 days of treatment start",
            type = "inclusion",
            num_cases = n_distinct(cohort_df$condition_id)) %>%
    mutate(loss = num_cases - lag(num_cases))
  
  cxr_microscopy <- cxr_manual %>% inner_join(microscopy %>% select(condition_id, specimen_collection_date, microscopyresults)) %>%
    rowwise() %>%
    mutate(timika_score = ifelse(if_any(.cols = ends_with("cavities")), overallpercentofabnormalvolume + 40, overallpercentofabnormalvolume)) %>%
    filter(between(imaging_date, specimen_collection_date - 30, specimen_collection_date)) %>%
    group_by(condition_id) %>%
    filter(imaging_date == max(imaging_date)) %>%
    slice_max(timika_score)
  
  cohort_df <- cohort_df[cohort_df$condition_id %in% cxr_microscopy$condition_id, ]
  
  counts_df %<>%
    add_row(step = "Has annotated CXR within 30 days prior to sputum microscopy",
            type = "inclusion",
            num_cases = n_distinct(cohort_df$condition_id)) %>%
    mutate(loss = num_cases - lag(num_cases))
  
  cohort_df %<>%
    filter(type_of_resistance %in% c("Sensitive", "MDR non XDR"))
  
  counts_df %<>%
    add_row(step = "Type of resistance = Sensitive, MDR non XDR",
            type = "inclusion",
            num_cases = n_distinct(cohort_df$condition_id)) %>%
    mutate(loss = num_cases - lag(num_cases))
  
  cxr_microscopy %<>%
    mutate(upper = if_any(starts_with("upper") & where(is.numeric), ~.x > 0) | if_any(starts_with("upper") & where(is.factor), ~.x != "No"),
           middle = if_any(starts_with("middle") & where(is.numeric), ~.x > 0) | if_any(starts_with("middle") & where(is.factor), ~.x != "No"),
           lower = if_any(starts_with("lower") & where(is.numeric), ~.x > 0) | if_any(starts_with("lower") & where(is.factor), ~.x != "No"),
           cavity = if_any(ends_with("cavities"), ~.x > 0),
           nodule = if_any(ends_with("nodules"), ~.x > 0))

  cxr_microscopy %<>% inner_join(cohort_df %>% select(-one_of(c("ispleuraleffusionbilateral", "rater")))) %>%
    select(-matches("qure"))
  
  cxr_microscopy_counts_df <- counts_df

  cxr_qure <- endpoints_dfs[[6]] %>%
    mutate(across(where(is.character) & starts_with("qure"), ~factor(.x, levels = c("No", "Yes"))))
  
  cxr_microscopy %<>% left_join(cxr_qure)
  
  # Outcome prediction
  
  cxr_microscopy_outcome <- cxr_microscopy %>%
    filter(outcome %in% c("Completed", "Cured", "Died", "Failure"))
  
  counts_df %<>%
    add_row(step = "Has treatment outcome of Completed, Cured, Died, or Failure",
            type = "inclusion",
            num_cases = n_distinct(cxr_microscopy_outcome$condition_id)) %>%
    mutate(loss = num_cases - lag(num_cases))
  
  # Matching of age, gender across type of resistance
  
  m <- cxr_microscopy_outcome %>% ungroup() %>% 
    select(condition_id, imagingstudy_id, microscopyresults, outcome, age_of_onset,
           gender, type_of_resistance, hiv, case_definition,
           timika_score, bmi, risk_smoker, risk_drug) %>%
    mutate(bmi = case_when(bmi < 18.5 ~ "underweight",
                           is.na(bmi) ~ "missing",
                           between(bmi, 18.5, 25) ~ "healthy",
                           T ~ "overweight/obese")) %>%
    mutate(across(where(is.character) & -ends_with("_id"), ~as.factor(.x)))
  
  # After matching
  m.balanced <- matchit(type_of_resistance ~ age_of_onset + gender + hiv + case_definition + timika_score + bmi + risk_smoker + risk_drug,
                        data = m, method = "cem", cutpoints = list(age_of_onset = 5,
                                                         timika_score = 5))
  
  cxr_microscopy_matched <- match.data(m.balanced)
  
  cxr_microscopy_matched <- cxr_microscopy_outcome %>% filter(condition_id %in% cxr_microscopy_matched$condition_id)
  
  counts_df %<>%
    add_row(step = "Age, gender, hiv, case definition, timika score, bmi, risk smoker, risk drug matched",
            type = "inclusion",
            num_cases = n_distinct(cxr_microscopy_matched$condition_id)) %>%
    mutate(loss = num_cases - lag(num_cases))
  
  cxr_microscopy_outcome_qure <- cxr_microscopy_matched %>% filter(!is.na(qure_abnormal))
  
  counts_df %<>%
    add_row(step = "Has qure AI annotations",
            type = "inclusion",
            num_cases = n_distinct(cxr_microscopy_outcome_qure$condition_id)) %>%
    mutate(loss = num_cases - lag(num_cases))
  
  return(list("cxr_microscopy_df" = cxr_microscopy, "cxr_microscopy_counts" = cxr_microscopy_counts_df,
              "cxr_microscopy_outcome_df" = cxr_microscopy_outcome,
              "cxr_microscopy_matched_df" = cxr_microscopy_matched,
              "cxr_microscopy_outcome_qure_df" = cxr_microscopy_outcome_qure,
              "cxr_microscopy_outcome_qure_counts" = counts_df,
              "matched_summary" = m.balanced))

}
