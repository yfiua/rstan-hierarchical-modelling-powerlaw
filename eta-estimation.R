library(dplyr)
library(rstan)
library(poweRlaw)

# read params from CLI
args <- commandArgs(trailingOnly = TRUE)
N_per_id <- as.integer(args[1])  # expected number of samples per unique ID

# options for rstan
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# Set random seed for reproducibility
set.seed(42)

# params
I <- 10             # number of unique IDs
mu_true <- 1.0      # true mu
sigma_true <- 0.1   # true sigma
tau_min <- 1        # minimum tau

# calculate total N
N <- I / 2 * N_per_id

# generate latent eta for each ID
eta <- rnorm(I, mean = mu_true, sd = sigma_true)

## save etas to a file
write.csv(eta, file = "sim_pl_eta_reliability.csv")

# read stan model
stan_model_hierarchical_pl <- stan_model("hierarchical-powerlaw.stan", model_name = "hierarchical-powerlaw")

for (iter in 1:100) {
  # generate pairs of IDs
  id_a0 <- sample(1:I, N, replace = TRUE)
  id_b0 <- sapply(id_a0, function(x) {
    sample(setdiff(1:I, x), 1)   # ensure id_a0 != id_b0
  })

  ## ensure id_a < id_b
  id_a <- pmin(id_a0, id_b0)
  id_b <- pmax(id_a0, id_b0)

  # Get group IDs for each pair
  id_group <- paste(id_a, id_b, sep = "_")
  id_group_all <- unique(id_group)
  id_group <- match(id_group, id_group_all)

  # empty tau vector
  tau <- numeric(N)

  for (i in 1:I) {
    for (j in 1:I) {
      # compute alpha for id_a == i & id_b == j
      alpha <- eta[i] + eta[j]

      # sample size
      n <- sum(id_a == i & id_b == j)

      # sample tau from power law distribution
      if (n > 0) {
        tau[id_a == i & id_b == j] <- rpldis(n = n, xmin = tau_min, alpha = alpha)
      }
    }
  }

  # make a dataframe
  df <- data.frame(id_a, id_b, id_group, tau)

  group_assignment <- df %>%
    group_by(id_group) %>%
    slice(1) %>%        # take the first row in each group
    ungroup() %>%
    select(id_group, id_a, id_b)

  # prepare data for Stan
  stan_data <- list(
    N = nrow(df),
    I = max(c(df$id_a, df$id_b)), # total unique IDs
    G = max(df$id_group),     # total unique group IDs
    id_a = df$id_a,
    id_b = df$id_b,
    id_group = df$id_group,
    tau = df$tau,
    # group assignments (id_group, id_a, id_b)
    id_a_assignment = group_assignment$id_a,
    id_b_assignment = group_assignment$id_b
  )

  # fit stan model
  fit_hierarchical_pl <- sampling(
    stan_model_hierarchical_pl,
    data = stan_data,
    iter = 2000,
    chains = 4,
    seed = 123,
    refresh = 0
  )

  cat(summary(fit_hierarchical_pl)$summary[ paste0("eta[", 1:10, "]"), "mean" ], sep = ', ')
  cat(', ')
  cat(summary(fit_hierarchical_pl)$summary[ paste0("eta[", 1:10, "]"), "sd" ], sep = ', ')
  cat('\n')
}
