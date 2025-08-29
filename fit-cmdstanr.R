library(dplyr)
library(cmdstanr)
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

# read stan model to CmdStanR
cmdstan_model_pl <- cmdstan_model("powerlaw.stan")
cmdstan_model_hierarchical_pl <- cmdstan_model("hierarchical-powerlaw.stan")
cmdstan_model_hierarchical_pl_groupwise <- cmdstan_model("hierarchical-powerlaw-groupwise.stan")

# fit the models
fit_pl <- cmdstan_model_pl$sample(
  data = stan_data,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  seed = 123,
  save_warmup = TRUE
)

fit_hierarchical_pl <- cmdstan_model_hierarchical_pl$sample(
  data = stan_data,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  seed = 123,
  save_warmup = TRUE
)

fit_hierarchical_pl_groupwise <- cmdstan_model_hierarchical_pl_groupwise$sample(
  data = stan_data_groupwise,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  seed = 123,
  save_warmup = TRUE
)

# bridge samplers - not working directly with CmdStanR
# set.seed(456) # for reproducibility
# 
# bs_pl <- bridge_sampler(fit_pl, silent = TRUE, repetitions = 10)
# bs_hierarchical_pl <- bridge_sampler(fit_hierarchical_pl, silent = TRUE, repetitions = 10)
# bs_hierarchical_pl_groupwise <- bridge_sampler(fit_hierarchical_pl_groupwise, silent = TRUE, repetitions = 10)
# 
# # print log marginal likelihoods nicely
# print("Log marginal likelihoods for power-law model:")
# print(bs_pl$logml)
# print("Log marginal likelihoods for hierarchical power-law model")
# print(bs_hierarchical_pl$logml)
# print("Log marginal likelihoods for hierarchical power-law model (groupwise)")
# print(bs_hierarchical_pl_groupwise$logml)

# get mean parameter from fit stan model
print(summary(fit_hierarchical_pl)$summary[ paste0("eta[", 1:10, "]"), "mean" ])
