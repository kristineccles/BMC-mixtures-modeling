library(dplyr)
library(tidyr)
library(truncnorm)

set.seed(1332)

# Simulate parameters for 5 individual chemicals
chemical_names <- paste0("Chem", 1:5)
n_doses <- 8
dose_levels <- 10^seq(-5, 2, length.out = n_doses)

# Set one common top value for all chemicals
common_top <- runif(1, 85, 100)

sim_data <- lapply(chemical_names, function(chem) {
  # Use common top value for all
  top <- common_top
  ec50 <- runif(1, 1, 20)       # narrower EC50 range (less spread)
  slope <- runif(1, 0.5, 1.1)   # less variable slope

  # Calculate responses using Hill equation with lower noise
  resp <- top / (1 + (ec50 / dose_levels)^slope) + rnorm(n_doses, 0, 2)
  resp <- pmax(pmin(resp, 100), 0)  # constrain between 0 and 100

  data.frame(
    Chemicalname = chem,
    Dose_uM = dose_levels,
    MaxResp = resp,
    group = "individual"
  )
}) %>% bind_rows()

# Keep mix_ratios defined for potential future use
mix_ratios <- list(
  EM   = c(0.2, 0.2, 0.2, 0.2, 0.2),
  EXP1 = c(0.4, 0.3, 0.1, 0.1, 0.1),
  EXP2 = c(0.1, 0.1, 0.4, 0.2, 0.2),
  EXP3 = c(0.3, 0.1, 0.3, 0.1, 0.2)
)

# Final dataset: only individual chemicals
df <- sim_data
