
#' Generate time to event sample with censoring
#'
rsurv <- function(n = 100,
                  distn = "exp",
                  prop_cens = 0.1,
                  params = list(mu = c(-5, 0.005)),
                  X = rep(c(10, 25, 50, 100), each = 25)) {

  times <-
    if (distn == "exp") {
      rexp_rgn(n, params$mu, X)
    } else if (distn == "weibull") {
      rweibull_rgn(n, params$alpha, param$mu, X)
    } else if (distn == "bi-weibull") {
      rbiweibull_rgn(n, params$alpha, params$mu_exp, params$mu_w, X)
    }

  # uniform sample then censor selected
  cens_idx <- sample(1:n,
                     size = n*prop_cens,
                     replace = FALSE)

  t_cens <- times
  t_cens[cens_idx] <-
    map_dbl(t_cens[cens_idx], function(x) runif(1, 0, x))

  status <- as.numeric(!(1:n) %in% cens_idx)

  list(times = times,
       t_cens = t_cens,
       status = status,
       X = X)
}


#'
rsurv_mix <- function(cf = 0.2,
                      n = 200,
                      distn = c("exp", "exp"),
                      prop_cens = 0.1, # c(0.1, 0.2)
                      X = X,
                      params =
                        list(
                          list(
                            mu = c(2.5, 0.005)),
                          list(
                            mu = c(-8, 0.005))),
                      X = rep(c(10, 25, 50, 100), each = 25)) {

  z <- rbinom(n, 1, cf) + 1
  s <- sum(z == 1)
  s[2] <- sum(z == 2)
  m <- X[z == 1]
  m[2] <- X[z == 2]

  res <- list()

  for (i in seq_along(distn)) {

    res[[i]] <-
      rsurv(n = s[i],
            X = m[i],
            params[[i]],
            distn = distn[i],
            prop_cens = prop_cens)
  }

  list(
    times = c(res[[1]]$times, res[[2]]$times),
    times_cens = c(res[[1]]$times_cens, res[[2]]$times_cens),
    status = c(res[[1]]$status, res[[2]]$status)
  )
}

