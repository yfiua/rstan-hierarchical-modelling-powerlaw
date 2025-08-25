data {
    int<lower=1> N;         // Number of observations
    int<lower=1> I;         // Total number of unique group IDs
    array[N] int<lower=1, upper=I> id_group; // ID of the group (id_a, id_b), commutative
    array[N] int<lower=1> tau;       // Observed data
}

parameters {
    vector<lower=0>[I] eta; // eta associated with each id_group
    //real<lower=0> sigma;       // Stddev of eta prior
}

transformed parameters {
  // estimate the zeta function using Kahan summation
  array[I] real z;

  for (i in 1:I) {
    real s = 0;
    real c = 0;   // compensation
    for (k in 1:1000) {
      real term = exp(-eta[i] * log(k));  // numerically stable version of pow(k, -eta[i])
      real y_k = term - c;
      real t = s + y_k;
      c = (t - s) - y_k;
      s = t;
    }

    // Store the zeta function approximation for each group
    z[i] = s;
  }
}

model {
    // Priors
    // normal dist for eta
    //eta ~ normal(1, sigma);

    for (n in 1:N) {
        target += -eta[id_group[n]] * log(tau[n]) - log(z[id_group[n]]);       // power law
    }
}
