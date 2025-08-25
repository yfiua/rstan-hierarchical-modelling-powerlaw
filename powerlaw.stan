data {
  int<lower=1> N;
  array[N] int<lower=1> tau;  // observed durations
}

parameters {
  real alpha;   // exponent
}

transformed parameters {
  // estimate the zeta function for alpha using Kahan summation
  real zeta_approx;
  {
    real s = 0;
    real c = 0;   // compensation
    for (k in 1:1000) {
      real term = exp(-alpha * log(k));  // numerically stable version of pow(k, -alpha)
      real y_k = term - c;
      real t = s + y_k;
      c = (t - s) - y_k;
      s = t;
    }
    zeta_approx = s;
  }
}

model {
  // target
  target += -alpha * sum(log(to_vector(tau))) - N * log(zeta_approx);

  // alpha ~ normal(2, 1);
}
