#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param cohort_df
#' @param visualizations_tables
#' @param best_itt_comparisons
make_flextable_images <- function(cohort_df, visualizations_tables, best_itt_comparisons) {

  table1 = cohort_df$cxr_microscopy_outcome_qure_counts %>% dplyr::filter(!grepl("qure", step))%>% flextable() %>% theme_booktabs() %>% fontsize(size = 8, part = "all")
  table2 = visualizations_tables$unmatched_cohort_desc %>% select(-p.value) %>% rename(covariate=1) %>%
    rename_all(~gsub("MDR.non.XDR", "MDR", .x)) %>% rename_all(~gsub("\\.\\.", " \\(", .x)) %>%
    rename_all(~gsub("N\\.", "N = ", .x)) %>% rename_all(~gsub("\\.", "\\)", .x))%>% flextable() %>% theme_booktabs() %>% fontsize(size = 8, part = "all")
  table3 = visualizations_tables$matched_cohort_desc %>% select(-p.value) %>% rename(covariate=1) %>%
    rename_all(~gsub("MDR.non.XDR", "MDR", .x)) %>% rename_all(~gsub("\\.\\.", " \\(", .x)) %>%
    rename_all(~gsub("N\\.", "N = ", .x)) %>% rename_all(~gsub("\\.", "\\)", .x))%>% flextable() %>% theme_booktabs() %>% fontsize(size = 8, part = "all")
  sup_table1 = visualizations_tables$cv_summary %>% mutate(model = gsub("_outcome", "", model))%>% flextable() %>% theme_booktabs() %>% fontsize(size = 8, part = "all")
  sup_table2 = visualizations_tables$all_best %>% select(-c(probability:size))%>% flextable() %>% theme_booktabs() %>% fontsize(size = 8, part = "all")
  sup_table3 = best_itt_comparisons %>% mutate(outcome = rep(c("baseline_sputum", "treatment", "timika outcome comparison"), 2)) %>% select(-c(probability:size))%>% flextable() %>% theme_booktabs() %>% fontsize(size = 8, part = "all")
  
  l <- ls()
  l <- l[!grepl("cohort_df|visualizations|best_itt", l)]
  names(l) <- l
  l <- map(l ,~get0(.x))
  
  # Create save outputs
  if(dir.exists("plots")){
    for(n in names(l)){
      save_as_image(x = l[[n]], path = glue::glue("plots/{n}.png"))
    }
  }else{
    usethis::use_directory("plots")
    for(n in names(l)){
      save_as_image(x = l[[n]], path = glue::glue("plots/{n}.png"))
    }
  }
  
  l <- list.files("plots", full.names = T)
  names(l) <- basename(l)
  
  return(l)

}
