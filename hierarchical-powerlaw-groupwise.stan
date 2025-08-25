data {
    int<lower=1> N;         // Number of observations
    int<lower=1> I;         // Total number of unique group IDs
    array[N] int<lower=1, upper=I> id_group; // ID of the group (id_a, id_b), commutative
    array[N] int<lower=1> tau;       // Observed data
}

transformed data {
    // log(1), log(2), ..., log(K)
    int K = 1000;
    vector[K] logk;
    for (k in 1:K) logk[k] = log(k);
}

parameters {
    vector<lower=0>[I] eta; // eta associated with each id_group
    //real<lower=0> sigma;       // Stddev of eta prior
}

transformed parameters {
  // estimate the log zeta function
  vector[I] log_z;

  for (i in 1:I) {
    // Store the log zeta function approximation
    log_z[i] = log_sum_exp(-eta[i] * logk);
  }
}

model {
    // Priors
    // normal dist for eta
    //eta ~ normal(1, sigma);

    for (n in 1:N) {
        target += -eta[id_group[n]] * log(tau[n]) - log_z[id_group[n]];       // power law
    }
}
