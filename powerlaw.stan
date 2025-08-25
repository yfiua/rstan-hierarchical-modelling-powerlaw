data {
  int<lower=1> N;
  array[N] int<lower=1> tau;  // observed durations
}

transformed data {
    // log(1), log(2), ..., log(K)
    int K = 1000;
    vector[K] logk;
    for (k in 1:K) logk[k] = log(k);
}

parameters {
  real alpha;   // exponent
}

transformed parameters {
  // estimate the zeta function for alpha
  real log_z = log_sum_exp(-alpha * logk);
}

model {
  // target
  target += -alpha * sum(log(to_vector(tau))) - N * log_z;

  // alpha ~ normal(2, 1);
}
