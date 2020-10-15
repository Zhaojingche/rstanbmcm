// gompertz mixture cure model
// relative survival

// user-defined functions ----
functions {
#include /include/distributions.stan
}

// input data ----
data {
  int<lower=0> n;             // number of observations
  vector[n] t;                // observed times
  vector[n] d;                // censoring indicator (1 = observed, 0 = censored)
  int H;                      // number of covariates
  matrix[n,H] X;              // matrix of covariates (with n rows and H columns)

  // intercept only -
  // real mu_beta;	              // means of the covariates coefficients
  // real mu_bg;
  // real<lower=0> sigma_beta;    // sds of the covariates coefficients
  // real<lower=0> sigma_bg;
  // intercept and gradient -
  vector[H] mu_beta;
  vector<lower=0> [H] sigma_beta;
  vector[H] mu_bg;
  vector<lower=0> [H] sigma_bg;

  vector[H] mu_alpha;
  vector<lower=0> [H] sigma_alpha;

  real a_cf;                  // cure fraction ~ Beta(a,b)
  real b_cf;
}

parameters {
  vector[H] alpha0;
  vector[H] beta0;         // coefficients in linear predictor (including intercept)
  vector[H] beta_bg;
  real<lower=0, upper=1> curefrac;  //TODO: define as simplex?
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
  beta0 ~ normal(mu_beta, sigma_beta);
  beta_bg ~ normal(mu_bg, sigma_bg);
  alpha0 ~ normal(mu_alpha, sigma_alpha);

  curefrac ~ beta(a_cf, b_cf);

  for (i in 1:n) {

    target += log_sum_exp(log(curefrac) +
                surv_exp_lpdf(t[i] | d[i], lambda_bg[i]),
                log1m(curefrac) +
                surv_exp_lpdf(t[i] | d[i], lambda_bg[i]) +
                surv_gompertz_lpdf(t[i] | d[i], alpha0[i], lambda0[i]));
  }
}

generated quantities {
  real rate0;
  real rate_bg;
  vector[60] S_bg;
  vector[60] S0;
  vector[60] S_pred;

  rate0 = exp(beta0[1]);
  rate_bg = exp(beta_bg[1]);

  for (i in 1:60) {
    S_bg[i] = exp_Surv(i, rate_bg);
    S0[i] = gompertz_Surv(i, mu_alpha, rate0);
    S_pred[i] = curefrac*S_bg[i] + (1 - curefrac)*S0[i]*S_bg[i];
  }

  //// for each individual
  // matrix[60,n] S_bg;
  // matrix[60,n] S0;
  // matrix[60,n] S_pred;
  //
  // for (j in 1:n) {
    //   for (i in 1:60) {
      //     S_bg[i,j] = exp_Surv(i, lambda_bg[j]);
      //     S0[i,j] = gompertz_Surv(i, lambda_bg[j] + lambda0[j]);
      //     S_pred[i,j] = curefrac*S_bg[i,j] + (1 - curefrac)*S0[i,j];
      //   }
      // }
      //TODO: case-mix average at time, sample?
}

