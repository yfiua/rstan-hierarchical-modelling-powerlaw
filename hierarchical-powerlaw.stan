data {
    int<lower=1> N;            // Number of observations
    int<lower=1> I;         // Total number of unique IDs
    int<lower=1> G;         // Total number of unique group IDs
    array[N] int<lower=1, upper=I> id_a; // First ID per observation
    array[N] int<lower=1, upper=I> id_b; // Second ID per observation
    array[N] int<lower=1, upper=G> id_group;      // ID of the group (id_a, id_b), with id_a < id_b
    array[N] int<lower=1> tau;       // Observed data
    //# group assignments (group_id, i, j)
    array[G] int<lower=1, upper=I> id_a_assignment;
    array[G] int<lower=1, upper=I> id_b_assignment;
}

transformed data {
    // log(1), log(2), ..., log(K)
    int K = 1000;
    vector[K] logk;
    for (k in 1:K) logk[k] = log(k);
}

parameters {
    vector<lower=0>[I] eta; // eta associated with each ID
    // real<lower=0> sigma;       // Stddev of eta prior
}

transformed parameters {
    // estimate the zeta function
    vector[G] log_z;
    vector[G] alpha;
  
    for (i in 1:G) {
        alpha[i] = eta[id_a_assignment[i]] + eta[id_b_assignment[i]]; // alpha for each group
  
        // Store the log zeta function approximation
        log_z[i] = log_sum_exp(-alpha[i] * logk);
    }
}

model {
    // Priors
    // normal dist for eta
    // eta ~ normal(1, sigma);

    for (n in 1:N) {
        target += -alpha[id_group[n]] * log(tau[n]) - log_z[id_group[n]];       // power law
    }
}
