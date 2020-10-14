
#' Run stan mixture cure model
#'
#' @param input_data Dataframe with all variables
#' @param model Distribution name; exp, weibull, gompertz
#' @param event_type Overall survival or progression-free survival; OS, PFS
#' @param tx_name Treatment name; IPILIMUMAB, NIVOLUMAB, NIVOLUMAB+IPILIMUMAB
#' @param n_iter int
#' @param n_warmup int
#' @param n_chains Number of chains; int
#' @param mean_cf Mean cure fraction; default to U[0,1]
#' @param var_cf Variance of cure fraction
#' @param ... Additional arguments
#'
#' @import rstan
#' @import dplyr
#'
bmcm_stan <- function(input_data,
                      model = "exp",
                      event_type = "PFS",
                      tx_name = "IPILIMUMAB",
                      n_iter = 2000,
                      n_warmup = 1000,
                      n_chains = 2,
                      mean_cf = NA,
                      var_cf = NA,
                      ...) {

  event_type <- match.arg(arg = event_type, c("PFS", "OS"))
  tx_name <- match.arg(arg = tx_name,
                       c("IPILIMUMAB", "NIVOLUMAB", "NIVOLUMAB+IPILIMUMAB"))

  beta_params <-
    if (!is.na(mean_cf) && !is.na(var_cf)) {
      MoM_beta(mean_cf, var_cf)
    } else {
      list(a = 1, b = 1)}

  stan_model <-
    switch(model,
           exp      = stanmodels$exp_relative_mix,
           weibull  = stanmodels$weibull_relative_mix,
           gompertz = stanmodels$gompertz_relative_mix)

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

  data_list <-
    list(
      n = nrow(tx_dat[[tx_name]]),
      t = tx_dat[[tx_name]][[2]],
      d = tx_dat[[tx_name]][[3]],
      H = 2,
      X = matrix(c(rep(1, nrow(tx_dat[[tx_name]])),
                   tx_dat[[tx_name]][[4]]),
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

  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  # stan_rdump(c("n_obs", "y"), file = "mix.data.R")

  rstan::sampling(
    stan_model,
    data = data_list,
    warmup = n_warmup,
    control = list(adapt_delta = 0.99,
                   max_treedepth = 20),
    iter = n_iter,
    chains = n_chains, ...)
}

