// exponential mixture cure model
// relative survival

functions {
#include /include/distributions.stan
}

// input data ----
data {
  int<lower=0> n;             // number of observations
  vector[n] t;                // observed times
  vector[n] d;                // censoring indicator (1 = observed, 0 = censored)
  int<lower = 0> H;           // number of covariates
  matrix[n,H] X;              // matrix of covariates (with n rows and H columns)

  vector[H] mu_0;
  vector<lower=0> [H] sigma_0;
  vector[H] mu_bg;
  vector<lower=0> [H] sigma_bg;

  real<lower=0, upper=1> curefrac;

  int<lower=0> t_max;
}

parameters {
  vector[H] beta0;         // coefficients in linear predictor (including intercept)
  vector[H] beta_bg;
}

transformed parameters {
  vector[n] linpred0;
  vector[n] linpred_bg;
  vector[n] lambda0;
  vector[n] lambda_bg;

  linpred0 = X*beta0;
  linpred_bg = X*beta_bg;

  // rate parameters
  lambda0 = exp(linpred0);
  lambda_bg = exp(linpred_bg); // background survival with uncertainty
}

model {
  beta0 ~ normal(mu_0, sigma_0);
  beta_bg ~ normal(mu_bg, sigma_bg);

  for (i in 1:n) {
    target += log_sum_exp(curefrac +
                surv_exp_lpdf(t[i] | d[i], lambda_bg[i]),
                (1 - curefrac) +
                surv_exp_lpdf(t[i] | d[i], lambda_bg[i] + lambda0[i]));
  }
}

generated quantities {
  real rate0;
  real rate_bg;
  vector[t_max] S_bg;
  vector[t_max] S0;
  vector[t_max] S_pred;

  # intercept
  rate0 = exp(beta0[1]);
  rate_bg = exp(beta_bg[1]);

  for (i in 1:t_max) {
    S_bg[i] = exp_Surv(i, rate_bg);
    S0[i] = exp_Surv(i, rate_bg + rate0);
    S_pred[i] = curefrac*S_bg[i] + (1 - curefrac)*S0[i];
  }
}

