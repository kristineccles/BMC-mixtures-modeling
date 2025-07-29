################################################
# Mixture Modeling BMCs - individual chems
# Written By: Kristin Eccles
# Date: June 11, 2025
# Note: Needs 01-traditional-individual-chem.R to run
#################################################

#Load libraries
library(tcplfit2)
library(tidyr)

df <- read.csv("PAC_responses_6.csv")
#convert to uM
df$Dose_uM <- (df$Dose_M*1e6)

model_bmd_calc <- list()

for (chem in unique(df$Chemicalname )) {

  # Subset the data for this chemical
  chem_data <- (df[df$Chemicalname == chem, ])
  #chem_data <- subset(data_unlist, id %in% chem)

  # Define BMR setup
  row <- list(
    conc = chem_data$Dose_uM,
    resp = chem_data$MaxResp,
    bmed = 0,        # Baseline centered at 0
    cutoff = 0,   # cut off for hit call
    onesd = 10,       #10% increase = BMR but need to set bmr_scale = 1
    assay = "assay",
    name = chem
  )

  # Run Hill model + BMC calculation
  res <- concRespCore(
    row,
    fitmodels = c("hill"),
    conthits = TRUE,    # This enables BMC calculation
    aicc = FALSE,
    bidirectional = FALSE,
    errfun = "dnorm",
    bmr_scale = 1
  )

  # Store result
  model_bmd_calc[[chem]] <- res[,c("name", "assay", "bmd", "ac10", "bmdl", "bmdu", "tp")]
}

bmd_melt <-reshape2::melt(model_bmd_calc)

# Use dcast to reshape
tidy_bmc<- na.omit(reshape2::dcast(bmd_melt, L1 ~ variable, value.var = "value"))
tidy_bmc$include <- tidy_bmc$L1 %in% MIX_FRACTIONS$Chemical
bmc_individual <- subset(tidy_bmc, include == TRUE)
bmc_mix <- subset(tidy_bmc, include == FALSE)
bmc_individual <- left_join(bmc_individual, MIX_FRACTIONS, by = c("L1" = "Chemical"))


ggplot(tidy_bmc, aes(x = L1, y = log10(bmd), color = L1)) +
  geom_point(size = 3) +  # Dot for BMD
  geom_errorbar(aes(ymin = log10(bmdl), ymax = log10(bmdu)), width = 0.2) +  # Whiskers
  coord_flip() +
  theme_minimal() +
  labs(
    x = "Chemical",
    y = "Log10 BMC",
    color = "Chemical"
  )

##### Bootstrap ####
bootstrap_results <- list()

# Number of bootstrap samples
n_samples <- 1000

# Generate bootstrapped BMDs for each chemical using truncated normal
for (i in 1:nrow(bmc_individual)) {
  chem <- bmc_individual$L1[i]
  mean_bmd <- bmc_individual$bmd[i]
  lower <- bmc_individual$bmdl[i]
  upper <- bmc_individual$bmdu[i]

  # Estimate standard deviation (approximate from CI assuming normality)
  est_sd <- (upper - lower) / (2 * 1.96)

  # Generate truncated normal samples
  bmd_samples <- rtruncnorm(n_samples, a = lower, b = upper, mean = mean_bmd, sd = est_sd)

  # Store with chemical name
  bootstrap_results[[chem]] <- bmd_samples
}

bootstrap_df <- stack(bootstrap_results)
colnames(bootstrap_df) <- c("bmd_boot", "chemical")

