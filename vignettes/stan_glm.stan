data {
  int<lower=0> N; // number of observations
  int<lower=0> M; // number of mutations

  vector[N] y;    // dependent variable
  matrix[N, M] x; // matrix of predictors (including intercept)
}
parameters {
  vector[M] beta;      // coefficient estimates
  real<lower=0> sigma; // variance
}
transformed parameters {
  vector[N] eta;
  eta = x*beta;
}
model {
  y ~ normal(eta, sigma);
}
