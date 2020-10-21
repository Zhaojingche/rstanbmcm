
#'
prep_stan_params <- function(model) {

  if (model == "exp") {
    return(
      list(mu_beta = c(0.01, 0.01),
           sigma_beta = c(1,1)))
  }
  if (model == "weibull") {
    return(
      list(mu_alpha = 0.1,
           sigma_alpha = 0.1,
           mu_beta = c(0.01, 0.01),
           sigma_beta = c(1,1)))
  }

  stop("distribution not found.")
}

