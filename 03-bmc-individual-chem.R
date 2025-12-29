################################################
# Mixture Modeling BMCs - individual chems
# Written By: Kristin Eccles
# Date: June 11, 2025
# Note: Needs 01-traditional-individual-chem.R to run
#################################################
# Load libraries
library(tcplfit2)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)

# Store results
model_bmd_calc <- list()

# Loop over chemicals and BMR levels (10–50)
for (chem in unique(df$Chemicalname)) {
  chem_data <- df[df$Chemicalname == chem, ]

  for (bmr_val in seq(10, 50, 10)) {
    row <- list(
      conc = chem_data$Dose_uM,
      resp = chem_data$MaxResp,
      bmed = 0,
      cutoff = 0,
      onesd = bmr_val,   # vary BMR
      assay = "assay",
      name = chem
    )

    res <- concRespCore(
      row,
      fitmodels = c("hill"),
      conthits = TRUE,
      aicc = FALSE,
      bidirectional = FALSE,
      errfun = "dnorm",
      bmr_scale = 1
    )

    # Add a column to track BMR level
    res$bmr_level <- paste0("BMC", bmr_val)

    model_bmd_calc[[paste(chem, bmr_val, sep = "_")]] <-
      res[, c("name", "assay", "bmd", "bmdl", "bmdu", "tp", "bmr_level")]
  }
}

# Fill missing CI with point estimate
model_bmd_calc <- lapply(model_bmd_calc, function(df) {
  if (anyNA(df$bmdl)) df$bmdl <- ifelse(is.na(df$bmdl), df$bmd, df$bmdl)
  if (anyNA(df$bmdu)) df$bmdu <- ifelse(is.na(df$bmdu), df$bmd, df$bmdu)
  df
})

# Melt and reshape
bmd_melt <- reshape2::melt(model_bmd_calc)
tidy_bmc <- na.omit(
  reshape2::dcast(bmd_melt, name + bmr_level ~ variable, value.var = "value")
)

# Merge with mixing fractions
bmc_individual <- left_join(tidy_bmc, MIX_FRACTIONS, by = c("name" = "Chemical"))

# Plot example for BMC10 only
indiv_bmc <- ggplot(filter(tidy_bmc, bmr_level == "BMC10"),
       aes(x = name, y = log10(bmd), color = name)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = log10(bmdl), ymax = log10(bmdu)), width = 0.2) +
  coord_flip() +
  scale_color_viridis_d(option = "D") +
  theme_minimal() +
  labs(x = "Chemical", y = "Log10 BMC10", color = "Chemical")

##### Bootstrap ####
bootstrap_results <- list()
n_samples <- MCiter

for (i in 1:nrow(bmc_individual)) {
  chem <- bmc_individual$name[i]
  mean_bmd <- bmc_individual$bmd[i]
  lower <- bmc_individual$bmdl[i]
  upper <- bmc_individual$bmdu[i]

  est_sd <- (upper - lower) / (2 * 1.96)

  bmd_samples <- rtruncnorm(n_samples, a = lower, b = upper,
                            mean = mean_bmd, sd = est_sd)

  bootstrap_results[[paste0(chem, "_", bmc_individual$bmr_level[i])]] <- bmd_samples
}

bootstrap_df <- stack(bootstrap_results)
colnames(bootstrap_df) <- c("bmd_boot", "chemical_bmr")

bootstrap_df <- bootstrap_df %>%
  # split "Chem1_BMC10" -> chemical = Chem1, bmr_level = BMC10
  separate(chemical_bmr, into = c("chemical", "bmr_level"), sep = "_(?=[^_]+$)") %>%
  # clean bmr_level: "BMC10" -> "10"
  mutate(
    bmr_level = gsub("BMC", "", bmr_level),
    bmr_col   = paste0("bmd_boot_", bmr_level)
  ) %>%
  group_by(chemical, bmr_col) %>%
  mutate(iter = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c(chemical, iter),
    names_from = bmr_col,
    values_from = bmd_boot
  ) %>%
  arrange(chemical, iter)