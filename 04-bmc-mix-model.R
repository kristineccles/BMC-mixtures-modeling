################################################
# BMC Comparison: Measured, Concentration Addition, IA
# Includes individual chemicals and mixtures
# Author: Kristin Eccles
# Date: June 2025
################################################

library(dplyr)
library(reshape2)
library(ggplot2)
library(MCMCglmm)

set.seed(2353)
MCiter <- 1000

# prepare data
MIX_FRACTIONS <- na.omit(read.csv("mixing_fractions.csv"))
tidy_bmc_perc <- left_join(bootstrap_df, MIX_FRACTIONS, by = c("chemical" = "Chemical"))


### ---- Function: Concentration Addition ----
calculate_ca_bmd <- function(df_mixture, mix_cols, MCiter = 10000) {
  results <- list()

  for (mix_col in mix_cols) {

    # Step 1: Generate bootstrap BMD values for each chemical
    bmd_boot <- sapply(1:nrow(df_mixture), function(i) {
      rtruncnorm(
        n = MCiter,
        a = df_mixture$bmdl[i],
        b = df_mixture$bmdu[i],
        mean = df_mixture$bmd[i],
        sd = (df_mixture$bmdu[i] - df_mixture$bmdl[i]) / (2 * 2.2)
      )
    })

    colnames(bmd_boot) <- df_mixture$L1

    # Step 2: Melt to long format
    boot_melt <- reshape2::melt(bmd_boot, varnames = c("iter", "L1"), value.name = "bmd_value") %>%
      mutate(mix_fraction = rep(df_mixture[[mix_col]], each = MCiter))

    # Step 3: Compute CA BMD for each iteration
    ca_bmd_iter <- boot_melt %>%
      group_by(iter) %>%
      summarise(ca_bmd = 1 / sum(mix_fraction / bmd_value), .groups = "drop")

    # Step 4: Calculate quantiles
    results[[mix_col]] <- data.frame(
      bmd = quantile(ca_bmd_iter$ca_bmd, 0.5),
      bmdl = quantile(ca_bmd_iter$ca_bmd, 0.025),
      bmdu = quantile(ca_bmd_iter$ca_bmd, 0.975),
      L1 = sub("_percent$", "", mix_col)
    )
  }

  # Combine results from all mix columns
  bind_rows(results)
}
### ---- Function: Approximate IA ----
calculate_ia_bmd <- function(df_mixture, mix_cols, MCiter = 10000) {
  results <- list()

  for (mix_col in mix_cols) {

    # Step 1: Generate bootstrap BMD values for each chemical
    bmd_boot <- sapply(1:nrow(df_mixture), function(i) {
      rtruncnorm(
        n = MCiter,
        a = df_mixture$bmdl[i],
        b = df_mixture$bmdu[i],
        mean = df_mixture$bmd[i],
        sd = (df_mixture$bmdu[i] - df_mixture$bmdl[i]) / (2 * 2.2)
      )
    })

    colnames(bmd_boot) <- df_mixture$L1

    # Step 2: Melt to long format
    boot_melt <- reshape2::melt(bmd_boot, varnames = c("iter", "L1"), value.name = "bmd_value") %>%
      mutate(mix_fraction = rep(df_mixture[[mix_col]], each = MCiter))

    # Step 3: Compute CA BMD for each iteration
    ia_bmd_iter <- boot_melt %>%
      group_by(iter) %>%
      summarise(ia_bmd = sum(mix_fraction * bmd_value), .groups = "drop")

    # Step 4: Calculate quantiles
    results[[mix_col]] <- data.frame(
      bmd = quantile(ia_bmd_iter$ia_bmd, 0.5),
      bmdl = quantile(ia_bmd_iter$ia_bmd, 0.025),
      bmdu = quantile(ia_bmd_iter$ia_bmd, 0.975),
      L1 = sub("_percent$", "", mix_col)
    )
  }

  # Combine results from all mix columns
  bind_rows(results)
}

#### Predictions ####

mix_cols <- c("EM_percent", "EXP1_percent", "EXP2_percent", "EXP3_percent")
CA_BMC <- calculate_ca_bmd(df_mixture, mix_cols, MCiter = 1000)%>%
  mutate(method = "BMC CA")
IA_BMC <- calculate_ia_bmd(df_mixture, mix_cols, MCiter = 1000)%>%
  mutate(method = "BMC IA")

measured <- bmc_mix %>%
  dplyr::select(bmd, bmdl, bmdu, L1) %>%
  mutate(method = "BMC Measured")

combined_df <- rbind(CA_BMC, IA_BMC, measured)
# Export
write.csv(combined_df, "combined_bmc_IA_CA_Pred.csv", row.names = FALSE)

### ---- Plot: All Chemicals + Mixtures ----

# Plot

p1 <- ggplot(data = combined_df, aes(y = L1, x = log10(bmd), color = L1)) +
  geom_point(aes(shape = method), size = 3) +
  geom_errorbar(aes(xmin = log10(bmdl), xmax = log10(bmdu)), width = 0.1) +
  theme_bw() +
  labs(
    shape = "Model",
    color = "Mixture",
    x = "Log10 BMC (uM)",
    y = "Mixture"
  )
p1

ggsave("bmc_plot_final_DAIA.tiff", plot = p1, device = "tiff",
       width = 5, height = 7, units = "in", dpi = 300)
