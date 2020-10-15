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

  real a_cf;                  // cure fraction ~ Beta(a,b)
  real b_cf;
  // vector[n] h_bg;
}

parameters {
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
  lambda_bg = exp(linpred_bg); // background survival with uncertainty
  // lambda_bg = h_bg;           // _known_ point estimate for background survival

  //TODO:
  // lambda_bg = is_fixed ? h_bg : exp(linpred_bg);
}

model {
  beta0 ~ normal(mu_beta, sigma_beta);
  beta_bg ~ normal(mu_bg, sigma_bg);

  curefrac ~ beta(a_cf, b_cf);

  for (i in 1:n) {

    // target += log_mix(curefrac,
    //                   surv_exp_lpdf(t[i] | d[i], lambda_bg[i]),
    //                   surv_exp_lpdf(t[i] | d[i], lambda_bg[i] + lambda0[i]));

    // equivalently
    target += log_sum_exp(log(curefrac) +
                surv_exp_lpdf(t[i] | d[i], lambda_bg[i]),
                log1m(curefrac) +
                surv_exp_lpdf(t[i] | d[i], lambda_bg[i] + lambda0[i]));
  }
}

generated quantities {
  real rate0;
  real rate_bg;
  vector[60] S_bg;
  vector[60] S0;
  vector[60] S_pred;

  # intercept
  rate0 = exp(beta0[1]);
  rate_bg = exp(beta_bg[1]);

  for (i in 1:60) {
    S_bg[i] = exp_Surv(i, rate_bg);
    S0[i] = exp_Surv(i, rate_bg + rate0);
    S_pred[i] = curefrac*S_bg[i] + (1 - curefrac)*S0[i];
  }

  // for (i in 1:n) {    //TODO:
  //   ppv[i] = exp_mix_rng(curefrac, lambda0[i], lambda_bg[i])
  // }

  //// for each individual
  // matrix[60,n] S_bg;
  // matrix[60,n] S0;
  // matrix[60,n] S_pred;
  //
  // for (j in 1:n) {
    //   for (i in 1:60) {
      //     S_bg[i,j] = exp_Surv(i, lambda_bg[j]);
      //     S0[i,j] = exp_Surv(i, lambda_bg[j] + lambda0[j]);
      //     S_pred[i,j] = curefrac*S_bg[i,j] + (1 - curefrac)*S0[i,j];
      //   }
      // }
      //TODO: case-mix average at time, sample?
}

