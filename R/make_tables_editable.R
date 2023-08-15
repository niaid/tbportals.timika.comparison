#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param cohort_df
#' @param visualizations_tables
#' @param best_itt_comparisons
make_tables_editable <- function(cohort_df, visualizations_tables, best_itt_comparisons) {
  
  table1 = cohort_df$cxr_microscopy_outcome_qure_counts %>% dplyr::filter(!grepl("qure", step))
  table2 = visualizations_tables$unmatched_cohort_desc %>% select(-p.value) %>% rename(covariate=1) %>%
    rename_all(~gsub("MDR.non.XDR", "MDR", .x)) %>% rename_all(~gsub("\\.\\.", " \\(", .x)) %>%
    rename_all(~gsub("N\\.", "N = ", .x)) %>% rename_all(~gsub("\\.", "\\)", .x))
  table3 = visualizations_tables$matched_cohort_desc %>% select(-p.value) %>% rename(covariate=1) %>%
    rename_all(~gsub("MDR.non.XDR", "MDR", .x)) %>% rename_all(~gsub("\\.\\.", " \\(", .x)) %>%
    rename_all(~gsub("N\\.", "N = ", .x)) %>% rename_all(~gsub("\\.", "\\)", .x))
  sup_table1 = visualizations_tables$cv_summary %>% mutate(model = gsub("_outcome", "", model))
  sup_table2 = visualizations_tables$all_best %>% select(-c(probability:size))
  sup_table3 = best_itt_comparisons %>% mutate(outcome = rep(c("baseline_sputum", "treatment", "timika outcome comparison"), 2)) %>% select(-c(probability:size))
  
  l <- ls()
  l <- l[!grepl("cohort_df|visualizations|best_itt", l)]
  names(l) <- l
  l <- map(l ,~get0(.x))
  
  l <- map(l, ~mutate(.x, across(where(is.numeric), function(x) round(x, digits = 2))))
  
  # Create save outputs
  if(dir.exists("editable_tables")){
    for(n in names(l)){
      write_csv(x = l[[n]], file = glue::glue("editable_tables/{n}.csv"))
    }
  }else{
    usethis::use_directory("editable_tables")
    for(n in names(l)){
      write_csv(x = l[[n]], file = glue::glue("editable_tables/{n}.csv"))
    }
  }
  
  l <- list.files("editable_tables", full.names = T)
  names(l) <- basename(l)
  
  return(l)
  
}