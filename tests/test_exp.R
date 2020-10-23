
# https://betanalpha.github.io/assets/case_studies/identifying_mixture_models.html


library(purrr)
library(reshape2)
library(dplyr)
library(rstan)
library(shinystan)
library(dplyr)
library(ggplot2)
# library(rstanbmcm)
devtools::load_all()


cf <- 0.2
# n <- 200
n <- 20                     #smaller amount of data
z <- rbinom(n, 1, cf) + 1
n_bg <- sum(z == 1)
n_0 <- sum(z == 2)

age <- rep(c(10,25,50,100), each = n/4)

beta_bg <- 0.005
rates_bg <- exp(-8 + beta_bg*age[z == 1])
times_bg <- rexp(n_bg, rate = rates_bg)

beta_0 <- 0.005
rates_0 <- exp(-5 + beta_0*age[z == 2])
times_0 <- rexp(n_0, rate = rates_0)

# why signs apparently wrong way around??
lm(log(times_bg) ~ age[z == 1])
lm(log(times_0) ~ age[z == 2])

data_list <-
  list(mu_0 = c(-5, beta_0),
       sigma_0 = c(1,1),
       mu_bg = c(-8, beta_bg),
       sigma_bg = c(1,1),
       a_cf = 3,
       b_cf = 12,
       n = n,
       t = c(times_bg, times_0),
       d = rep(1, n),
       H = 2,
       X = matrix(c(rep(1, n),
                    age - mean(age)),
                  byrow = FALSE,
                  ncol = 2),
       t_max = 60)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)

res <-
  rstan::sampling(
    stanmodels$exp_relative_mix,
    # stanmodels$exp_cf_fixed,
    data = data_list,
    warmup = 1000,
    iter = 10000,
    thin= 10,
    control = list(adapt_delta = 0.99,
                   max_treedepth = 20),
    chains = 2)

rstan::check_divergences(res)

## output

fit_stan <- extract(res)

mean(fit_stan$curefrac)
quantile(fit_stan$curefrac, c(0.025, 0.5, 0.975))


par(mfrow = c(2,2))
hist(fit_stan$beta0[,1])
hist(fit_stan$beta0[,2])
hist(fit_stan$beta_bg[,1])
hist(fit_stan$beta_bg[,2])

## survival curves

par(mfrow = c(1,3))

# posteriors
plot(colMeans(fit_stan$S_pred), type = "l", ylim = c(0,1))
lines(colMeans(fit_stan$S_bg), type = "l", col = "red")
lines(colMeans(fit_stan$S0), type = "l", col = "blue")

# lines(apply(fit_stan$S0, 2, function(x) quantile(x, probs = 0.025)),
# type = "l", col = "blue")

# kaplan-meier
library(survival)
plot(
  survfit(Surv(times_bg, rep(1, n_bg)) ~ 1),
  col = "red",
  xlim = c(0, 60))
lines(
  survfit(Surv(times_0, rep(1, n_0)) ~ 1),
  col = "blue")
lines(
  survfit(Surv(c(times_bg, times_0), rep(1, n)) ~ 1))

# for each ages in sample
plot(NA, type = 'n', xlim = c(0,60), ylim = c(0,1))
lapply(rates_0, function(x) lines(times, exp(-times * x), col = "blue"))
lapply(rates_bg, function(x) lines(times, exp(-times * x), col = "red"))

# means of generation model
times <- 0:60
lambda_0 <- exp(-5 + beta_0*mean(age))
lines(times, exp(-times * lambda_0), col = "blue", lty = 2)

lambda_bg <- exp(-8 + beta_bg*mean(age))
lines(times, exp(-times * lambda_bg), col = "red", lty = 2)

## chains ----
par(mfrow=c(1,1))
xx <- as.data.frame(extract(res, permuted = FALSE)[,1,])
plot(exp(xx$`beta0[1]`), exp(xx$`beta_bg[1]`))

# highly correlated...
plot(xx$`beta0[1]`, xx$`beta0[2]`)

plot(xx$`beta_bg[1]`, xx$`beta_bg[2]`)
plot(xx$`beta_bg[1]`, xx$curefrac)
plot(xx$`beta0[1]`, xx$curefrac)
