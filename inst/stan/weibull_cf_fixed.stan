// weibull mixture cure model
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

  // intercept only -
  // real mu_beta;	              // means of the covariates coefficients
  // real mu_bg;
  // real<lower=0> sigma_beta;    // sds of the covariates coefficients
  // real<lower=0> sigma_bg;
  // intercept and gradient -
  real mu_alpha;
  real<lower=0> sigma_alpha;
  vector[H] mu_0;
  vector<lower=0> [H] sigma_0;
  vector[H] mu_bg;
  vector<lower=0> [H] sigma_bg;

  real<lower=0, upper=1> curefrac;
  int<lower=0> t_max;
}

parameters {
  real alpha0;
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
  lambda_bg = exp(linpred_bg);
}

model {
  beta0 ~ normal(mu_0, sigma_0);
  beta_bg ~ normal(mu_bg, sigma_bg);

  // not dependent on X or tranformed
  alpha0 ~ gamma(mu_alpha, sigma_alpha);

  for (i in 1:n) {

    target += log_sum_exp(log(curefrac) +
                surv_exp_lpdf(t[i] | d[i], lambda_bg[i]),
                log1m(curefrac) +
                joint_exp_weibull_lpdf(t[i] | d[i], alpha0, lambda0[i], lambda_bg[i]));
  }
}

generated quantities {
  real rate0;
  real rate_bg;
  vector[t_max] S_bg;
  vector[t_max] S0;
  vector[t_max] S_pred;

  rate0 = exp(beta0[1]);
  rate_bg = exp(beta_bg[1]);

  for (i in 1:t_max) {
    S_bg[i] = exp_Surv(i, rate_bg);
    S0[i] = weibull_Surv(i, alpha0, rate0);
    S_pred[i] = curefrac*S_bg[i] + (1 - curefrac)*S0[i]*S_bg[i];
  }
}

