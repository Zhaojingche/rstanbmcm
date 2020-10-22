
# run stan mixture cure model
# CheckMate 067 dataset


# kill chains with:
# system("killall R")
# system("taskkill /F /IM R.exe /T")

library(purrr)
library(reshape2)
library(dplyr)
library(rstan)
library(shinystan)
library(dplyr)
library(ggplot2)
# library(rstanbmcm)
devtools::load_all()

# surv_input_data
# load("C:/Users/Nathan/Documents/R/mixture_cure_model/data/surv_input_data.RData")
data("surv_input_data")

all_tx_names <- c("IPILIMUMAB", "NIVOLUMAB", "NIVOLUMAB+IPILIMUMAB")
all_event_types <- c("PFS", "OS")
model_names <- c("exp")#, "weibull", "gompertz")

## choose compiled stan?
# stan_fn <- bmcm_stan_file
stan_fn <- bmcm_stan

stan_out <- list()

for (k in model_names) {
  for (i in all_event_types) {
    for (j in all_tx_names) {
      # tryCatch({

      out <- stan_fn(input_data = surv_input_data,
                     model = k,
                     event_type = i,
                     tx_name = j,
                     warmup = 100,
                     iter = 10000,#50000,
                     thin = 10)#50)

      file_name <-
        here::here(paste("data/stan", k, i, j, ".Rds", sep = "_"))

      stan_files[[k]][[i]][[j]] <- file_name
      saveRDS(out, file = file_name)

      # },
      # error = function(e) e)
    }
  }
}

# save(stan_files, file = "data/stan_filenames.RData")

plot_S_event_type(stan_files$exp)
plot_S_event_type(stan_files$weibull)


####################
# 10% cure fraction
####################

for (i in all_event_types) {
  for (j in all_tx_names) {
    tryCatch({
      stan_10[[i]][[j]] <-
        bmcm_stan(input_data = surv_input_data,
                  model = "exp",
                  event_type = i,
                  tx_name = j,
                  mean_cf = 0.1,
                  var_cf = 0.1)},
      error = function(e) e)
  }
}

save(stan_10, file = "data/stan_10.RData")

plot_S_event_type(stan_10)


####################
# 40% cure fraction
####################

for (i in all_event_types) {
  for (j in all_tx_names) {
    tryCatch({
      stan_40[[i]][[j]] <-
        bmcm_stan(input_data = surv_input_data,
                  model = "exp",
                  event_type = i,
                  tx_name = j,
                  mean_cf = 0.4,
                  var_cf = 0.1)},
      error = function(e) e)
  }
}

save(stan_40, file = "data/stan_40.RData")

plot_S_event_type(stan_40)


# -------------------------------------------------------------------------

# priors for cure fraction

beta_nivo <- c(-0.2, 0.001) # mean 45%
beta_ipi <- c(-1.2, 0.001)  # mean 23%
beta_both <- c(0.16, 0.001) # mean 54%

