
#' Run stan mixture cure model from stan file
#'
#' Does not use pre-compiled C code.
#'
#' @param input_data Dataframe with all variables
#' @param model Distribution name; exp, weibull, gompertz
#' @param event_type Overall survival or progression-free survival; OS, PFS
#' @param tx_name Treatment name; IPILIMUMAB, NIVOLUMAB, NIVOLUMAB+IPILIMUMAB
#' @param iter int
#' @param warmup Number of iterations burn-in; int
#' @param chains Number of chains; int
#' @param mean_cf Mean cure fraction; default to U[0,1]
#' @param var_cf Variance of cure fraction
#' @param centre_age Logical to centre regression covariate
#' @param ... Additional arguments
#'
#' @import rstan
#' @import dplyr
#'
bmcm_stan_file <- function(input_data,
                           model = "exp",
                           event_type = "PFS",
                           tx_name = "IPILIMUMAB",
                           iter = 3000,
                           warmup = 1000,
                           chains = 3,
                           thin = 10,
                           mean_cf = NA,
                           var_cf = NA,
                           centre_age = TRUE,
                           ...) {
  data_list <-
    c(prep_stan_params(model),
      prep_stan_data(input_data,
                     event_type,
                     tx_name,
                     centre_age,
                     mean_cf,
                     var_cf))

  stan_file <-
    switch(model,
           exp      = here::here("inst", "stan", "exp_relative_mix.stan"),
           weibull  = here::here("inst", "stan", "weibull_relative_mix.stan"),
           gompertz = here::here("inst", "stan", "gompertz_relative_mix.stan"))

  rstan_options(auto_write = TRUE)
  options(mc.cores = min(n_chains, parallel::detectCores() - 1))
  # stan_rdump(c("n_obs", "y"), file = "mix.data.R")

  res <-
    rstan::stan(
      file = stan_file,
      data = data_list,
      warmup = warmup,
      control = list(adapt_delta = 0.99,
                     max_treedepth = 20),
      iter = iter,
      chains = chains, ...)

  return(res)
}

