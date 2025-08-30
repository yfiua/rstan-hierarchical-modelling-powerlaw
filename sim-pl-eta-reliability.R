# set random seed for reproducibility
set.seed(42)

# params
I <- 10             # number of unique IDs
mu_true <- 1.0      # true mu
sigma_true <- 0.1   # true sigma

# generate latent eta for each ID
eta <- rnorm(I, mean = mu_true, sd = sigma_true)

# save etas to file
df <- as.data.frame(eta)
write.csv(df, file = "sim-pl-eta-reliability.csv", row.names = FALSE)
