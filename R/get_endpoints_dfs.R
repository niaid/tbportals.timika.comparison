#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param endpoints
get_endpoints_dfs <- function(endpoints) {
  
  t <- get_token()
  
  r <- tidy_depot_api(path = endpoints, token = t)
  
  df <- r$content %>%
    select(-in_requested_cohort)
  
  df
  
}
