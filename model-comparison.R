library(dplyr)
library(rstan)
library(readr)
library(bridgesampling)

input_file <- 'sim_pl.csv'
n_thres <- 0

# options for rstan
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# read data
data <- read_csv(input_file)

# filter
data <- data %>% add_count(id_group, name = 'n') %>% filter(n > n_thres)

# create unique IDs for groups
all_id_groups <- unique(data$id_group)
data <- data %>% mutate(id_group = match(id_group, all_id_groups))

# create unique IDs for individuals
all_ids <- unique(c(data$id_a, data$id_b))
data <- data %>%
    mutate(
      id_a = match(id_a, all_ids),
      id_b = match(id_b, all_ids)
    ) %>%
    rowwise() %>%       # ensure id_a < id_b
    mutate(sorted = list(sort(c(id_a, id_b)))) %>%
    mutate(id_a = sorted[1], id_b = sorted[2]) %>%
    select(-sorted) %>%
    ungroup()

group_assignment <- data %>%
    group_by(id_group) %>%
    slice(1) %>%        # take the first row in each group
    ungroup() %>%
    select(id_group, id_a, id_b)

# prepare data for Stan
stan_data <- list(
                  N = nrow(data),
                  I = max(c(data$id_a, data$id_b)), # total unique IDs
                  G = max(data$id_group),     # total unique group IDs
                  id_a = data$id_a,
                  id_b = data$id_b,
                  id_group = data$id_group,
                  tau = data$tau,
                  # group assignments (id_group, id_a, id_b)
                  id_a_assignment = group_assignment$id_a,
                  id_b_assignment = group_assignment$id_b
)

stan_data_groupwise <- list(
                            N = nrow(data),
                            I = max(data$id_group), # total unique group IDs
                            id_group = data$id_group,
                            tau = data$tau
)

# read stan model
stan_model_pl <- stan_model("edges-powerlaw.stan", model_name = "powerlaw")
stan_model_hierarchical_pl <- stan_model("edges-hierarchical-powerlaw.stan", model_name = "hierarchical-powerlaw")
stan_model_hierarchical_pl_2 <- stan_model("edges-hierarchical-powerlaw-2.stan", model_name = "hierarchical-powerlaw-2")
stan_model_hierarchical_pl_groupwise <- stan_model("edges-hierarchical-powerlaw-groupwise.stan", model_name = "hierarchical-powerlaw-groupwise")

# fit the models
fit_pl <- sampling(
  stan_model_pl,
  data = stan_data,
  iter = 2000,
  chains = 4,
  seed = 123
)

fit_hierarchical_pl <- sampling(
  stan_model_hierarchical_pl,
  data = stan_data,
  iter = 2000,
  chains = 4,
  seed = 123
)

fit_hierarchical_pl_2 <- sampling(
  stan_model_hierarchical_pl_2,
  data = stan_data,
  iter = 2000,
  chains = 4,
  seed = 123
)

fit_hierarchical_pl_groupwise <- sampling(
  stan_model_hierarchical_pl_groupwise,
  data = stan_data_groupwise,
  iter = 2000,
  chains = 4,
  seed = 123
)

# bridge samplers
set.seed(456) # for reproducibility

bs_pl <- bridge_sampler(fit_pl, silent = TRUE, repetitions = 10)
bs_hierarchical_pl <- bridge_sampler(fit_hierarchical_pl, silent = TRUE, repetitions = 10)
bs_hierarchical_pl_2 <- bridge_sampler(fit_hierarchical_pl_2, silent = TRUE, repetitions = 10)
bs_hierarchical_pl_groupwise <- bridge_sampler(fit_hierarchical_pl_groupwise, silent = TRUE, repetitions = 10)

# print log marginal likelihoods nicely
sprintf("Log marginal likelihoods for power-law model: %.4f", bs_pl$logml)
sprintf("Log marginal likelihoods for hierarchical power-law model: %.4f", bs_hierarchical_pl$logml)
sprintf("Log marginal likelihoods for hierarchical power-law model 2: %.4f", bs_hierarchical_pl_2$logml)
sprintf("Log marginal likelihoods for hierarchical power-law model (groupwise): %.4f", bs_hierarchical_pl_groupwise$logml)