################################################
# 02 Mixtures Modeling for simulated data
# Written By: Kristin Eccles
# Date: June 11th, 2025
# Updated July 28th, 2025
# Note: Needs 01-individual_chemical_dr.R to run
# must have created folder in project called "mix_pred_boot"
#################################################

# ----------------------------
# Libraries
# ----------------------------
library(MCMCglmm)
library(drc)
library(dplyr)
library(tidyr)
library(broom)
library(reshape2)
library(purrr)
library(truncnorm)   # rtruncnorm()
library(scales)      # rescale()
library(ggplot2)
library(viridis)

# ----------------------------
# Set Up
# ----------------------------
NUM_PTS   <- 100
MIN_LOGX  <- -5
MAX_LOGX  <- 5

MAX_Y <- 1
yVec  <- vector(mode = "numeric", length = 100 - 1)
y     <- data.frame(y = MAX_Y * 1:100 / 100)

MCiter <- 100

# x-grid for forward predictions
x <- 10^seq(MIN_LOGX, MAX_LOGX, length.out = NUM_PTS)

# ----------------------------
# Data filter for refits
# ----------------------------
df$include <- df$Chemicalname %in% individual_coeff_final$curve
df_reduced <- subset(df, include == TRUE)

# ----------------------------
# Helpers (DRY!)
# ----------------------------

# FIXED: robust coef extractor that survives fixed params (missing SEs)
extract_coefs <- function(model, se_terms = c("Slope", "ED50", "Upper Limit")) {
  wci <- broom::tidy(model, conf.int = TRUE)

  # estimates wide
  coef_wide <- wci %>%
    dplyr::select(term, curve, estimate) %>%
    tidyr::pivot_wider(names_from = term, values_from = estimate)

  # std.errors wide (may be missing for fixed params)
  se_wide <- wci %>%
    dplyr::select(term, curve, std.error) %>%
    dplyr::mutate(term = paste0("SE_", gsub(" ", "", term))) %>%
    tidyr::pivot_wider(names_from = term, values_from = std.error)

  # ensure requested SE columns exist
  needed_se_cols <- paste0("SE_", gsub(" ", "", se_terms))
  for (nm in needed_se_cols) if (!nm %in% names(se_wide)) se_wide[[nm]] <- NA_real_

  dplyr::left_join(coef_wide, se_wide, by = "curve")
}

# Truncated-normal bootstrap helper
boot_truncnorm <- function(m, s, lower, upper, n = MCiter) {
  rtruncnorm(n = n, a = lower, b = upper, mean = m, sd = s * sqrt(2))
}

# Summarize to median and 95% CI
summ_ci <- function(v) {
  c(mean  = quantile(v, 0.5,   na.rm = TRUE),
    lower = quantile(v, 0.025,  na.rm = TRUE),
    upper = quantile(v, 0.975,  na.rm = TRUE))
}

# Safe CSV writer
safe_write <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.csv(df, path, row.names = FALSE)
}

# Mix metadata
mix_meta <- tibble::tibble(
  mix = c("EM","EXP1","EXP2","EXP3"),
  col = c("EM_percent","EXP1_percent","EXP2_percent","EXP3_percent")
)

# ----------------------------
# Mixing ratios joiners / tops
# ----------------------------
MIX_FRACTIONS <- data.frame(
  Chemical     = chemical_names,
  EM_percent   = mix_ratios$EM,
  EXP1_percent = mix_ratios$EXP1,
  EXP2_percent = mix_ratios$EXP2,
  EXP3_percent = mix_ratios$EXP3
)

# Keep original 'tops' join (as in script)
tops <- dplyr::left_join(individual_model_coeff, MIX_FRACTIONS, by = c("curve" = "Chemical")) %>% na.omit()

# All four tops come from common_top (as in original)
EXP3_top <- common_top
EM_top   <- common_top
EXP1_top <- common_top
EXP2_top <- common_top
tops_lookup <- c(EM = EM_top, EXP1 = EXP1_top, EXP2 = EXP2_top, EXP3 = EXP3_top)

# ----------------------------
# Refit: Hill model free slope (for IA/CA)
# ----------------------------
individual_model2 <- drm(
  MaxResp ~ Dose_uM, data = df_reduced, curveid = Chemicalname,
  type = "continuous",
  lowerl = c(-Inf, 0),
  upperl = c(0, Inf),
  fct = LL.4(fixed = c(NA, FIXED_C, FIXED_D, NA),
             names = c("Slope", "Lower Limit", "Upper Limit", "ED50"))
)
summary(individual_model2)
plot(individual_model2)

individual_model_coeff2 <- extract_coefs(individual_model2, se_terms = c("Slope", "ED50", "Upper Limit"))

# Join mixing ratios
individual_model_coeff2 <- dplyr::left_join(
  individual_model_coeff2,
  MIX_FRACTIONS,
  by = c("curve" = "Chemical")
)

# ----------------------------
# Bootstrap slope & ED50 (for IA/CA)
# ----------------------------
coeff_by_chem <- split(individual_model_coeff2, individual_model_coeff2$curve)
chem_names    <- names(coeff_by_chem)

slope_boot <- lapply(chem_names, function(ch) {
  boot_truncnorm(coeff_by_chem[[ch]]$Slope,  coeff_by_chem[[ch]]$SE_Slope,  lower = -Inf, upper = 0,   n = MCiter)
}) %>% stats::setNames(chem_names)

ED50_boot <- lapply(chem_names, function(ch) {
  boot_truncnorm(coeff_by_chem[[ch]]$`ED50`, coeff_by_chem[[ch]]$SE_ED50,   lower = 0,    upper = Inf, n = MCiter)
}) %>% stats::setNames(chem_names)

slope_boot_m <- reshape2::melt(slope_boot)
ED50_boot_m  <- reshape2::melt(ED50_boot)

bootmat <- data.frame(
  ED50     = ED50_boot_m[, 1],
  slope    = slope_boot_m[, 1],
  chemical = slope_boot_m[, 2],
  itter    = rep(seq_len(MCiter), times = length(chem_names))
) %>%
  dplyr::left_join(
    individual_model_coeff2 %>% dplyr::select(curve, EM_percent, EXP1_percent, EXP2_percent, EXP3_percent),
    by = c("chemical" = "curve")
  ) %>%
  stats::na.omit()

bootmat_list <- split(bootmat, bootmat$itter)

# Sanity checks (kept)
sum(individual_model_coeff2$EM_percent)
sum(individual_model_coeff2$EXP1_percent)
sum(individual_model_coeff2$EXP2_percent)
sum(individual_model_coeff2$EXP3_percent)

# ----------------------------
# Refit: slope fixed (GCA)
# ----------------------------
individual_model_b1 <- drm(
  MaxResp ~ Dose_uM, data = df_reduced, curveid = Chemicalname,
  type = "continuous",
  lowerl = c(0, 0),
  upperl = c(100, Inf),
  fct = LL.4(fixed = c(FIXED_B, FIXED_C, NA, NA),
             names = c("Slope", "Lower Limit", "Upper Limit", "ED50"))
)
plot(individual_model_b1)

B1_coefs <- extract_coefs(individual_model_b1, se_terms = c("ED50", "Upper Limit"))

B1_by_chem <- split(B1_coefs, B1_coefs$curve)
B1_names   <- names(B1_by_chem)

Top_boot <- lapply(B1_names, function(ch) {
  boot_truncnorm(B1_by_chem[[ch]]$`Upper Limit`, B1_by_chem[[ch]]$SE_UpperLimit, lower = 0, upper = Inf, n = MCiter)
}) %>% stats::setNames(B1_names)

B1_ED50_boot <- lapply(B1_names, function(ch) {
  boot_truncnorm(B1_by_chem[[ch]]$`ED50`, B1_by_chem[[ch]]$SE_ED50, lower = 0, upper = Inf, n = MCiter)
}) %>% stats::setNames(B1_names)

Top_boot_m     <- reshape2::melt(Top_boot)
B1_ED50_boot_m <- reshape2::melt(B1_ED50_boot)

B1bootmat <- data.frame(
  ED50     = B1_ED50_boot_m[, 1],
  Top      = Top_boot_m[, 1],
  chemical = Top_boot_m[, 2],
  itter    = rep(seq_len(MCiter), times = length(B1_names))
) %>%
  dplyr::left_join(
    individual_model_coeff2 %>% dplyr::select(curve, EM_percent, EXP1_percent, EXP2_percent, EXP3_percent),
    by = c("chemical" = "curve")
  )

B1bootmat_list <- split(B1bootmat, B1bootmat$itter)

# ----------------------------
# IA predictions (loop over mixes)
# ----------------------------
# y(x) = 1 - prod_i [ 1 / (1 + exp(-b_i (log(x*w_i) - ED50_i))) ]
ia_for_mix <- function(mix_col) {
  lapply(seq_along(bootmat_list), function(j) {
    apply(matrix(x), 1, function(xv) {
      1 - prod(1 / (1 + exp(-bootmat_list[[j]]$slope * (log(xv * bootmat_list[[j]][[mix_col]]) - bootmat_list[[j]]$ED50))))
    })
  })
}

build_ia_df <- function(mix_name, top_val) {
  mix_col <- mix_meta$col[mix_meta$mix == mix_name]
  L       <- ia_for_mix(mix_col)
  unlist_df <- cbind(x, reshape2::melt(L))
  colnames(unlist_df) <- c("x", "y", "iteration")
  unlist_df$mix    <- mix_name
  unlist_df$method <- "IA"

  unlist_rescaled <- unlist_df %>% dplyr::mutate(y = scales::rescale(y, to = c(0, top_val)))
  safe_write(unlist_rescaled, file.path("mix_pred_boot", paste0("unlist_", mix_name, "_CI.csv")))

  ci <- unlist_rescaled %>%
    dplyr::group_by(x) %>%
    dplyr::summarise(
      mean    = quantile(y, 0.5,   na.rm = TRUE),
      y_lower = quantile(y, 0.025,  na.rm = TRUE),
      y_upper = quantile(y, 0.975,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dplyr::across(c(mean, y_lower, y_upper), ~ scales::rescale(., to = c(0, top_val)))) %>%
    as.data.frame()
  ci$mix_ratio <- mix_name
  ci
}

IA_list <- purrr::pmap(
  list(mix_meta$mix, tops_lookup[mix_meta$mix]),
  build_ia_df
)

IA_df <- dplyr::bind_rows(IA_list) %>% dplyr::mutate(Method = "IA")
safe_write(IA_df, "IA_df.csv")

# Plot IA
p1adj <- ggplot() +
  geom_line(data = IA_df, aes(x = log10(x), y = mean, color = mix_ratio), linewidth = 1) +
  geom_ribbon(data = IA_df, aes(x = log10(x), y = mean, ymin = y_lower, ymax = y_upper, fill = mix_ratio), alpha = 0.2) +
  theme_bw() +
  labs(x = "Log10 Concentration (uM)", y = "% Max MeBio Response",
       color = "Mixing Ratio \nIndependent Action",
       fill  = "Mixing Ratio \nIndependent Action") +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme(legend.title = element_text(face = "bold", size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title   = element_text(size = 16, hjust = 0.5))
p1adj

# ----------------------------
# CA predictions (loop over mixes)
# ----------------------------
# x(y) = sum_i w_i * ( ED50_i * exp( b_i * log((1 - y) / y) ) )
ca_for_mix <- function(mix_col) {
  lapply(seq_along(bootmat_list), function(j) {
    b  <- bootmat_list[[j]]$slope      # can be negative; that's fine
    w  <- bootmat_list[[j]][[mix_col]]
    ED <- bootmat_list[[j]]$ED50       # must be > 0

    apply(as.matrix(y$y), 1, function(yv) {
      sum( w * ( ED * exp( (1 / b) * log((1 - yv) / yv) ) ) )
    })
  })
}

build_ca_df <- function(mix_name, top_val) {
  mix_col <- mix_meta$col[mix_meta$mix == mix_name]
  L       <- ca_for_mix(mix_col)
  unlist_df <- cbind(y, reshape2::melt(L)) %>%
    dplyr::mutate(y = scales::rescale(y, to = c(0, top_val)))
  colnames(unlist_df) <- c("y", "x", "iteration")
  unlist_df$mix    <- mix_name
  unlist_df$method <- "CA"
  safe_write(unlist_df, file.path("mix_pred_boot", paste0("unlist_CA_", mix_name, "_CI.csv")))

  ci <- unlist_df %>%
    dplyr::group_by(y) %>%
    dplyr::summarise(
      mean   = quantile(x, 0.5,   na.rm = TRUE),
      x_lower= quantile(x, 0.025,  na.rm = TRUE),
      x_upper= quantile(x, 0.975,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(y = scales::rescale(y, to = c(0, top_val))) %>%
    as.data.frame()
  ci$mix_ratio <- mix_name
  ci
}

CA_list <- purrr::pmap(
  list(mix_meta$mix, tops_lookup[mix_meta$mix]),
  build_ca_df
)

CA_df <- dplyr::bind_rows(CA_list) %>% dplyr::mutate(Method = "CA")
safe_write(CA_df, "CA_df.csv")

# Plot CA
p2adj <- ggplot() +
  geom_line(data = CA_df, aes(x = log10(mean), y = y, color = mix_ratio), linewidth = 1) +
  geom_ribbon(data = CA_df, aes(x = log10(mean), y = y, xmin = log10(x_lower), xmax = log10(x_upper), fill = mix_ratio), alpha = 0.2) +
  theme_bw() +
  labs(x = "Log10 Concentration (uM)", y = "% Max MeBio Response",
       color = "Mixing Ratio \nConcentration Addition",
       fill  = "Mixing Ratio \nConcentration Addition") +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme(legend.title = element_text(face = "bold", size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title   = element_text(size = 16, hjust = 0.5))
p2adj

# ----------------------------
# GCA predictions (loop over mixes)
# ----------------------------
# y(x) = sum_i Top_i * ((x*w_i)/ED50_i) / [1 + sum_i ((x*w_i)/ED50_i)]
gca_for_mix <- function(mix_col) {
  lapply(seq_along(B1bootmat_list), function(j) {
    apply(matrix(x), 1, function(xv) {
      num <- sum(B1bootmat_list[[j]]$Top * ((xv * B1bootmat_list[[j]][[mix_col]]) / B1bootmat_list[[j]]$ED50))
      den <- 1 + sum((xv * B1bootmat_list[[j]][[mix_col]]) / B1bootmat_list[[j]]$ED50)
      num / den
    })
  })
}

build_gca_df <- function(mix_name) {
  mix_col <- mix_meta$col[mix_meta$mix == mix_name]
  L       <- gca_for_mix(mix_col)
  unlist_df <- cbind(x, reshape2::melt(L))
  colnames(unlist_df) <- c("x", "y", "iteration")
  unlist_df$mix    <- mix_name
  unlist_df$method <- "GCA"
  safe_write(unlist_df, file.path("mix_pred_boot", paste0("unlist_GCA_", mix_name, "_CI.csv")))

  ci <- unlist_df %>%
    dplyr::group_by(x) %>%
    dplyr::summarise(
      mean  = quantile(y, 0.5,   na.rm = TRUE),
      lower = quantile(y, 0.025,  na.rm = TRUE),
      upper = quantile(y, 0.975,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame()
  ci$mix_ratio <- mix_name
  ci
}

GCA_list <- purrr::map(mix_meta$mix, build_gca_df)

GCA_df <- dplyr::bind_rows(GCA_list) %>%
  dplyr::mutate(
    Method = "GCA",
    y       = mean,
    y_lower = lower,
    y_upper = upper
  )

safe_write(GCA_df, "GCA_df.csv")

# Plot GCA
p3adj <- ggplot() +
  geom_line(data = GCA_df, aes(x = log10(x), y = y, color = mix_ratio), linewidth = 1) +
  geom_ribbon(data = GCA_df, aes(x = log10(x), y = y, ymin = y_lower, ymax = y_upper, fill = mix_ratio), alpha = 0.2) +
  theme_bw() +
  labs(x = "Log10 Concentration (uM)", y = "% Max MeBio Response",
       color = "Mixing Ratio \nGCA",
       fill  = "Mixing Ratio \nGCA") +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme(legend.title = element_text(face = "bold", size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title   = element_text(size = 16, hjust = 0.5))
p3adj


