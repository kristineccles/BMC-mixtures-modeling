###################################################################
# Calculate BMCs and compare how effect level affects deviation from additivity
# Written By: Kristin Eccles
# Date: July 29th, 2025
# Note: Needs 04-bmc-mix-model ro run
###################################################################
# --- If bootstrap_df already exists in your workspace, you can skip the read_csv line ---
# library(readr)
# bootstrap_df <- readr::read_csv("/mnt/data/bootstrap_df.csv")

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)

# ----------------------------
# 1) Tidy bootstrap_df to long
# ----------------------------
# expects columns: chemical, iter, bmd_boot_10, bmd_boot_20, ...
boot_long <- bootstrap_df %>%
  pivot_longer(
    cols = starts_with("bmd_boot_"),
    names_to   = "BMR_lab",
    values_to  = "bmd_value"
  ) %>%
  mutate(
    # extract numeric BMR from "bmd_boot_XX"
    BMR = as.numeric(str_extract(BMR_lab, "\\d+"))
  ) %>%
  dplyr::select(chemical, iter, BMR, bmd_value)

# ----------------------------
# 2) Mix fractions to long
# ----------------------------
# MIX_FRACTIONS must have columns: Chemical, <mixture1>_percent, <mixture2>_percent, ...
mix_long <- MIX_FRACTIONS %>%
  pivot_longer(
    cols = ends_with("_percent"),
    names_to = "mixture",
    values_to = "mix_fraction"
  ) %>%
  mutate(mixture = str_remove(mixture, "_percent$")) %>%
  filter(!is.na(mix_fraction) & mix_fraction > 0) %>%
  rename(chemical = Chemical)

# ----------------------------
# 3) Join bootstraps to mix fractions
# ----------------------------
boot_mix <- boot_long %>%
  inner_join(mix_long, by = "chemical")
# now we have: chemical, iter, BMR, bmd_value, mixture, mix_fraction

# ----------------------------
# 4) Compute mixture BMC per iter for CA and IA
# ----------------------------
mix_CA <- boot_mix %>%
  group_by(mixture, BMR, iter) %>%
  summarise(
    bmd = 1 / sum(mix_fraction / bmd_value),
    .groups = "drop"
  ) %>%
  group_by(mixture, BMR) %>%
  summarise(
    bmd_mean  = median(bmd, na.rm = TRUE),
    bmdl = quantile(bmd, 0.025, na.rm = TRUE),
    bmdu = quantile(bmd, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(model = "CA")

mix_IA <- boot_mix %>%
  group_by(mixture, BMR, iter) %>%
  summarise(
    bmd = sum(mix_fraction * bmd_value),
    .groups = "drop"
  ) %>%
  group_by(mixture, BMR) %>%
  summarise(
    bmd_mean  = median(bmd, na.rm = TRUE),
    bmdl = quantile(bmd, 0.025, na.rm = TRUE),
    bmdu = quantile(bmd, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(model = "IA")

final_preds <- bind_rows(mix_CA, mix_IA)
# ==========================================================
# 5) Combine CA and IA predictions (no measured data)
# ==========================================================
# final_preds contains columns: mixture, BMR, bmd, bmdl, bmdu, model
bmc_combined <- final_preds

# ==========================================================
# 6) Plot: CA vs IA only
# ==========================================================
library(ggplot2)

p4 <- ggplot(bmc_combined, aes(y = BMR, x = log10(bmd_mean), color = model)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = log10(bmdl), xmax = log10(bmdu)), height = 1) +
  facet_wrap(~ mixture) +
  theme_bw() +
  labs(
    x = "Log10 BMC (µM)",
    y = "Benchmark Response (%)",
    title = "Mixture BMC at Increasing Effect Levels: CA vs IA"
  ) +
  scale_color_manual(
    name   = "Model",
    values = c("CA" = "#FDE725FF", "IA" = "#2A788EFF"),
    labels = c("CA", "IA")
  )

print(p4)
ggsave("bootstrapped_bmc_log10_CA_IA_plot.tiff",
       plot = p4, width = 8, height = 5, dpi = 300)

# ==========================================================
# 7) Comparison: log10(IA / CA) with bootstrap uncertainty
# ==========================================================
# Rebuild iteration-level predictions for both CA and IA
pred_iter_CA <- boot_mix %>%
  dplyr::group_by(mixture, BMR, iter) %>%
  dplyr::summarise(bmd = 1 / sum(mix_fraction / bmd_value), .groups = "drop") %>%
  dplyr::mutate(model = "CA")

pred_iter_IA <- boot_mix %>%
  dplyr::group_by(mixture, BMR, iter) %>%
  dplyr::summarise(bmd = sum(mix_fraction * bmd_value), .groups = "drop") %>%
  dplyr::mutate(model = "IA")

# Join CA and IA by mixture/BMR/iteration to compute ratio
comparison_iter <- dplyr::inner_join(
  pred_iter_CA,
  pred_iter_IA,
  by = c("mixture", "BMR", "iter"),
  suffix = c("_CA", "_IA")
) %>%
  dplyr::mutate(log10_diff = log10(bmd_IA / bmd_CA))

# Summarize per mixture & BMR
summary_stats <- comparison_iter %>%
  dplyr::group_by(mixture, BMR) %>%
  dplyr::summarise(
    median_log10_diff = median(log10_diff, na.rm = TRUE),
    lower = quantile(log10_diff, 0.025, na.rm = TRUE),
    upper = quantile(log10_diff, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# ==========================================================
# 8) Plot: log10(IA / CA)
# ==========================================================
p6 <- ggplot(summary_stats, aes(x = BMR, y = median_log10_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 3, color = "#2A788EFF") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 1, color = "#2A788EFF") +
  facet_wrap(~ mixture) +
  theme_bw() +
  labs(
    x = "Benchmark Response (%)",
    y = "Log10(IA BMC / CA BMC)",
    title = "Deviation from Additivity: IA vs CA"
  )

print(p6)
ggsave("bootstrapped_log10_IA_CA_diff_plot.png",
       plot = p6, width = 8, height = 5, dpi = 300)
