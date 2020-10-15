
# run stan mixture cure model
# CheckMate 067 dataset


library(rstan)
library(shinystan)
library(dplyr)
# library(rstanbmcm)

# surv_input_data
load("C:/Users/Nathan/Documents/R/mixture_cure_model/data/surv_input_data.RData")

all_tx_names <- c("IPILIMUMAB", "NIVOLUMAB", "NIVOLUMAB+IPILIMUMAB")
all_event_types <- c("PFS", "OS")
model_names <- "exp" #c("exp", "weibull", "gompertz")

# stan_fn <- bmcm_stan
stan_fn <- bmcm_stan_file

stan_out <- list()

for (k in model_names) {
  for (i in all_event_types) {
    for (j in all_tx_names) {
      tryCatch({
        stan_out[[k]][[i]][[j]] <-
          stan_fn(input_data = surv_input_data,
                  model = k,
                  event_type = i,
                  tx_name = j)},
        error = function(e) e)
    }
  }
}

# save(stan_out, file = "data/stan_out.RData")

plot_S_event_type(stan_out[[1]])


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

