# Load libraries
library(dplyr)
library(poweRlaw)

# Set random seed for reproducibility
set.seed(42)

# Define parameters
N <- 50000          # Number of observations
I <- 10             # Number of unique IDs
mu_true <- 1.0     # True mu
sigma_true <- 0.2  # True sigma
tau_min <- 1       # Minimum tau

# Generate latent eta for each ID
eta <- rnorm(I, mean = mu_true, sd = sigma_true)

# Generate pairs of IDs
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

# save etas to a file
write.csv(eta, file = "sim_pl_eta.csv")


# Pack into a data frame
synthetic_data <- data.frame(
  id_a = id_a,
  id_b = id_b,
  id_group = id_group,
  tau = tau
)

# View first few rows
head(synthetic_data)

# save data
write.csv(synthetic_data, "sim_pl.csv", row.names = FALSE)
