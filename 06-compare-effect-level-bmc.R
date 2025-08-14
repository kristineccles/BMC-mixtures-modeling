###################################################################
# Calculate BMCs and compare how effect level affects deviation from additivity
# Written By: Kristin Eccles
# Date: July 29th, 2025
# Note: Needs 04-bmc-mix-model ro run
###################################################################
# Load libraries
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(truncnorm)
library(purrr)
library(ggpmisc)

set.seed(4571)
MCiter <- 1000
BMR_levels <- seq(10, 70, by = 10)  # BMC10 to BMC70

### Load Data ####
# Load individual component estimates used in CA/IA predictions
df_mixture <- tidy_bmc_perc
# Identify mixture proportion columns (e.g., EM_percent, EXP1_percent)
mix_cols <- grep("_percent$", names(df_mixture), value = TRUE)

df <- read.csv("PAC_responses_6.csv")
df$Dose_uM <- df$Dose_M * 1e6
# Isolate only the mixture group
mixture_df <- df %>% filter(group == "mixture")

#####################################################################
# Iterate over the measured mixtures effect levels
# Initialize list to store all results
mixture_bmd_by_bmr <- list()

# Loop over each mixture and BMR level
for (mix in unique(mixture_df$Chemicalname)) {
  mix_data <- mixture_df %>% filter(Chemicalname == mix)

  for (bmr in BMR_levels) {

    # Define input for tcplfit2
    row <- list(
      conc = mix_data$Dose_uM,
      resp = mix_data$MaxResp,
      bmed = 0,
      cutoff = 0,
      onesd = bmr,
      assay = mix,  # Use actual mixture name
      name = paste0(mix, "_", bmr)
    )

    # Run the model
    res <- concRespCore(
      row,
      fitmodels = c("hill"),
      conthits = TRUE,
      aicc = FALSE,
      bidirectional = FALSE,
      errfun = "dnorm",
      bmr_scale = 1
    )

    # Store result
    mixture_bmd_by_bmr[[paste0(mix, "_", bmr)]] <- res[, c("name", "assay", "bmd", "bmdl", "bmdu", "tp")]
  }
}

# Convert to long then wide format
mixture_bmd_long <- reshape2::melt(mixture_bmd_by_bmr)
bmc_mix_all <- reshape2::dcast(mixture_bmd_long, L1 ~ variable, value.var = "value") %>%
  separate(L1, into = c("mixture", "BMR"), sep = "_(?=[0-9]+$)") %>%
  mutate(BMR = as.numeric(BMR))

# Plot: BMC vs BMR for all mixtures
ggplot(bmc_mix_all, aes(x = BMR, y = log10(bmd), color = mixture)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = log10(bmdl), ymax = log10(bmdu)), width = 1) +
  theme_bw() +
  labs(
    x = "Benchmark Response (%)",
    y = "Log10 BMC (µM)",
    title = "Mixture BMC at Increasing Effect Levels"
  )


#####################################################################
### CA and IA Prediction Functions ####

# Load required data
df <- read.csv("PAC_responses_6.csv")
df$Dose_uM <- df$Dose_M * 1e6
MIX_FRACTIONS <- na.omit(read.csv("mixing_fractions.csv"))
df_individual <- df %>% filter(group == "individual")

# Fit model and extract BMDs at multiple BMRs
library(tcplfit2)
model_bmd_list <- list()

for (bmr in BMR_levels) {
  temp_results <- list()

  for (chem in unique(df_individual$Chemicalname)) {
    chem_data <- df_individual[df_individual$Chemicalname == chem, ]
    row <- list(
      conc = chem_data$Dose_uM,
      resp = chem_data$MaxResp,
      bmed = 0,
      cutoff = 0,
      onesd = bmr,
      assay = "assay",
      name = chem
    )
    res <- concRespCore(row, fitmodels = "hill", conthits = TRUE, bmr_scale = 1)
    res$BMR <- bmr
    temp_results[[chem]] <- res[, c("name", "assay", "bmd", "bmdl", "bmdu", "tp", "BMR")]
  }

  model_bmd_list[[as.character(bmr)]] <- do.call(rbind, temp_results)
}

# Combine and prepare
bmd_all <- bind_rows(model_bmd_list)
bmd_melt <- reshape2::melt(split(bmd_all, bmd_all$name), id.vars = c("name", "assay", "bmd", "bmdl", "bmdu", "tp", "BMR"))
bmd_wide <- left_join(bmd_melt, MIX_FRACTIONS, by = c("name" = "Chemical"))
mix_cols <- grep("_percent$", names(bmd_wide), value = TRUE)

# Function for bootstrapped predictions
calc_bmc_mix_boot <- function(bmd_df_individual, mix_cols, MCiter, method = "CA") {
  results <- list()
  for (mix_col in mix_cols) {
    mix_df_individual <- bmd_df_individual %>% filter(!is.na(.data[[mix_col]]))
    bmd_boot <- sapply(1:nrow(mix_df_individual), function(i) {
      rtruncnorm(
        n = MCiter,
        a = mix_df_individual$bmdl[i],
        b = mix_df_individual$bmdu[i],
        mean = mix_df_individual$bmd[i],
        sd = (mix_df_individual$bmdu[i] - mix_df_individual$bmdl[i]) / (2 * 2.2)
      )
    })
    colnames(bmd_boot) <- mix_df_individual$name

    boot_melt <- melt(bmd_boot, varnames = c("iter", "name"), value.name = "bmd_value") %>%
      mutate(mix_fraction = rep(mix_df_individual[[mix_col]], each = MCiter))

    if (method == "CA") {
      mix_vals <- boot_melt %>%
        group_by(iter) %>%
        summarise(bmd = 1 / sum(mix_fraction / bmd_value), .groups = "drop")
    } else {
      mix_vals <- boot_melt %>%
        group_by(iter) %>%
        summarise(bmd = sum(mix_fraction * bmd_value), .groups = "drop")
    }

    results[[mix_col]] <- data.frame(
      bmd = quantile(mix_vals$bmd, 0.5),
      bmdl = quantile(mix_vals$bmd, 0.025),
      bmdu = quantile(mix_vals$bmd, 0.975),
      mixture = sub("_percent$", "", mix_col)
    )
  }

  bind_rows(results)
}

# Loop over BMRs
all_preds <- list()
for (bmr in BMR_levels) {
  subset_df_individual <- subset(bmd_wide, BMR == bmr)
  CA <- calc_bmc_mix_boot(subset_df_individual, mix_cols, MCiter, method = "CA") %>%
    mutate(model = "CA", BMR = bmr)
  IA <- calc_bmc_mix_boot(subset_df_individual, mix_cols, MCiter, method = "IA") %>%
    mutate(model = "IA", BMR = bmr)
  all_preds[[as.character(bmr)]] <- bind_rows(CA, IA)
}

final_preds <- bind_rows(all_preds)

# Save results
write.csv(final_preds, "bmc_predictions_all_BMRs.csv", row.names = FALSE)

# Optional: visualize
ggplot(final_preds, aes(x = BMR, y = log10(bmd), color = model)) +
  geom_point() +
  geom_line() +
  facet_wrap(~mixture) +
  theme_bw() +
  labs(x = "Effect Level (BMR%)", y = "Log10 BMC", title = "Mixture BMC Predictions Across Effect Levels")

# Add missing `model` column to bmc_mix_all
bmc_mix_all$model <- "Measured"

# Drop `tp` column if not needed for binding
bmc_mix_all <- bmc_mix_all[,c(mixture, BMR , bmd, bmdl, bmdu, model)]
final_preds <- final_preds %>% select(mixture, BMR, bmd, bmdl, bmdu, model)

# Now bind
bmc_combined <- bind_rows(bmc_mix_all, final_preds)

# Plot BMC vs BMR with model comparison
p4 <- ggplot() +
  # Measured BMCs
  geom_point(data = bmc_mix_all, aes(y = BMR, x = log10(bmd), color = model), size = 3) +
  geom_errorbarh(data = bmc_mix_all, aes(y = BMR, xmin = log10(bmdl), xmax = log10(bmdu), color = model), height = 1) +

  # Predicted BMCs (CA, IA)
  geom_point(data = final_preds, aes(y = BMR, x = log10(bmd), color = model), size = 3) +
  geom_errorbarh(data = final_preds, aes(y = BMR, xmin = log10(bmdl), xmax = log10(bmdu), color = model), height = 1) +

  facet_wrap(~mixture) +
  theme_bw() +
  labs(
    x = "Log10 BMC (µM)",
    y = "Benchmark Response (%)",
    title = "Mixture BMC at Increasing Effect Levels: Measured vs Model-Predicted"
  ) +
  scale_color_manual(
    name = "Method",
    values = c("CA" = "#FDE725FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"),
    labels = c("CA", "IA", "Measured")
  )
p4

ggsave("bootstrapped_bmc_log10_plot.tiff", plot = p4, width = 8, height = 5, dpi = 300)

##################################################################################
#### Compare Differences ####
# Join measured and predicted by mixture and BMR
comparison_df <- final_preds %>%
  left_join( bmc_mix_all,by = c("BMR", "mixture")) %>%
  mutate(
    log10_diff = log10(bmd.y) - log10(bmd.x),
    abs_diff = bmd.y - bmd.x
  )

# Set number of bootstrap iterations
n_boot <- 1000

# Initialize list to collect results
boot_results <- list()

# Loop over each row to perform bootstrapping
for (i in 1:nrow(comparison_df)) {

  pred_mean <- comparison_df$bmd.x[i]
  pred_lower <- comparison_df$bmdl.x[i]
  pred_upper <- comparison_df$bmdu.x[i]

  meas_mean <- comparison_df$bmd.y[i]
  meas_lower <- comparison_df$bmdl.y[i]
  meas_upper <- comparison_df$bmdu.y[i]

  # Estimate SDs assuming normality
  pred_sd <- (pred_upper - pred_lower) / (2 * 1.96)
  meas_sd <- (meas_upper - meas_lower) / (2 * 1.96)

  # Sample from truncated normal distributions
  pred_samples <- rtruncnorm(n_boot, a = pred_lower, b = pred_upper, mean = pred_mean, sd = pred_sd)
  meas_samples <- rtruncnorm(n_boot, a = meas_lower, b = meas_upper, mean = meas_mean, sd = meas_sd)

  # Compute log10(predicted / measured)
  log10_diff <- log10(pred_samples / meas_samples)

  # Store results
  boot_results[[i]] <- data.frame(
    log10_diff = log10_diff,
    model = comparison_df$model.x[i],
    mixture = comparison_df$mixture[i],
    BMR = comparison_df$BMR[i]
  )
}

# Combine into one dataframe
boot_df <- bind_rows(boot_results)

# Summarize per group
summary_stats <- boot_df %>%
  group_by(model, mixture, BMR) %>%
  summarise(
    median_log10_diff = median(log10_diff),
    lower = quantile(log10_diff, 0.025),
    upper = quantile(log10_diff, 0.975),
    .groups = "drop"
  )

# Plot: log10(Predicted / Measured) with bootstrapped 95% CI
p6 <- ggplot(summary_stats, aes(x = BMR, y = median_log10_diff, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 1) +
  facet_wrap(~mixture) +
  theme_bw() +
  labs(
    x = "Benchmark Response (%)",
    y = "Log10(Predicted BMC / Measured BMC)",
    title = "Model Prediction Error with Uncertainty (CA and IA vs Measured)"
  ) +
  scale_color_manual(
    name = "Method",
    values = c("CA" = "#FDE725FF", "IA" = "#2A788EFF"),
    labels = c("CA", "IA")
  )
p6

ggsave("bootstrapped_bmc_log10diff_plot.png", plot = p6, width = 8, height = 5, dpi = 300)
