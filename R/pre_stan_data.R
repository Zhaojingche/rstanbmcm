
#' pre_stan_data
#'
#' @param input_data
#' @param event_type
#' @param tx_name
#' @param centre_age
#' @param mean_cf
#' @param var_cf
#'
#' @return
#' @export
#'
pre_stan_data <- function(input_data,
                          event_type,
                          tx_name,
                          centre_age,
                          mean_cf,
                          var_cf) {

  event_type <- match.arg(arg = event_type, c("PFS", "OS"))
  tx_name <- match.arg(arg = tx_name,
                       c("IPILIMUMAB", "NIVOLUMAB", "NIVOLUMAB+IPILIMUMAB"))

  beta_params <-
    if (!is.na(mean_cf) && !is.na(var_cf)) {
      MoM_beta(mean_cf, var_cf)
    } else {
      list(a = 1, b = 1)}

  if (event_type == "PFS") {

    tx_dat <-
      input_data %>%
      select(TRTA, pfs, pfs_event, PFSage, PFS_rate) %>%
      mutate(PFS_rate =
               ifelse(PFS_rate == 0, 0.00001, PFS_rate)) %>% # replace 0
      split(input_data$TRTA)
  } else if (event_type == "OS") {

    tx_dat <-
      input_data %>%
      select(TRTA, os, os_event, OSage, OS_rate) %>%
      mutate(OS_rate =
               ifelse(OS_rate == 0, 0.00001, OS_rate)) %>%
      split(input_data$TRTA)
  }

  # centring
  age_adj <- ifelse(centre_age, mean(tx_dat[[tx_name]][[4]]), 0)

  list(
    n = nrow(tx_dat[[tx_name]]),
    t = tx_dat[[tx_name]][[2]],
    d = tx_dat[[tx_name]][[3]],
    H = 2,
    X = matrix(c(rep(1, nrow(tx_dat[[tx_name]])),
                 tx_dat[[tx_name]][[4]] - age_adj),
               byrow = FALSE,
               ncol = 2),
    mu_beta = c(0,0),
    sigma_beta = c(1,1),
    mu_bg = c(-8.25, 0.066),
    sigma_bg = c(1, 1),
    a_cf = beta_params$a,
    b_cf = beta_params$b#,
    # h_bg = tx_dat[[tx_name]][[5]]
  )
}

