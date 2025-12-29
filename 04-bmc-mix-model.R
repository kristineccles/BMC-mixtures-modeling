################################################
# BMC Comparison: Concentration Addition, IA
# Includes individual chemicals and mixtures
# Author: Kristin Eccles
# Date: June 27th, 2025
# Note: Needs 03-bmc-individual-chem.R to run
################################################

library(dplyr)
library(reshape2)
library(ggplot2)
library(rlang)
library(purrr)
library(tcplfit2)

MCiter <- 100  # number of bootstrap iterations

# prepare data
df_mixture <- left_join(bootstrap_df, MIX_FRACTIONS, by = c("chemical" = "Chemical"))%>%
  as.data.frame()

#### Function: Concentration Addition ####
calculate_ca_bmd <- function(df_mixture, mix_cols, bmr_col = "bmd_boot_10") {
  results <- list()

  for (mix_col in mix_cols) {
    ca_iter <- df_mixture %>%
      group_by(iter) %>%
      rename(mix_fraction = !!sym(mix_col), bmd = !!sym(bmr_col)) %>%

      summarise(
        ca_bmd = 1 / sum(mix_fraction / bmd),
        .groups = "drop"
      )

    results[[mix_col]] <- data.frame(
      bmd  = quantile(ca_iter$ca_bmd, 0.5, na.rm = TRUE),
      bmdl = quantile(ca_iter$ca_bmd, 0.025, na.rm = TRUE),
      bmdu = quantile(ca_iter$ca_bmd, 0.975, na.rm = TRUE),
      L1   = sub("_percent$", "", mix_col),
      method = "BMC CA"
    )
  }
  bind_rows(results)
}

#### Function: Approximate IA ####
# Check Concentration Alignmnet #
check_dose_alignment <- function(df) {
  # Extract unique doses per chemical
  dose_lists <- df %>%
    group_by(Chemicalname) %>%
    summarise(Doses = list(sort(unique(Dose_uM))), .groups = "drop")

  # Use first chemical as reference
  ref_doses <- dose_lists$Doses[[1]]
  ref_chem <- dose_lists$Chemicalname[1]

  # Compare all chemicals to reference
  dose_check <- dose_lists %>%
    mutate(
      matches_ref = map_lgl(Doses, ~ identical(.x, ref_doses)),
      missing_in_ref = map(Doses, ~ setdiff(ref_doses, .x)),
      extra_in_chem = map(Doses, ~ setdiff(.x, ref_doses))
    )

  # Report summary
  if (all(dose_check$matches_ref)) {
    message("TRUE: All chemicals have identical concentration levels.")
  } else {
    message("FALSE: Concentrations do NOT match across all chemicals. Details below:")
    print(dose_check)
  }

  invisible(dose_check)
}
# If FALSE then a concentration grid will need to be used to estimate the response at unmeasured concentrations
dose_check_results <- check_dose_alignment(df)

#### Calculate IA ####
# ---- Step 1: scale responses per chemical ----
scaled_df <- df %>%
  group_by(Chemicalname) %>%
  mutate(
    min_resp = min(MaxResp, na.rm = TRUE),
    max_resp = max(MaxResp, na.rm = TRUE),
    ScaledResp = ifelse(
      max_resp > min_resp,
      (MaxResp - min_resp) / (max_resp - min_resp),
      0
    )
  ) %>%
  ungroup()

# ---- Step 2: compute Independent Action mixture responses ----
calc_mixture <- function(mix_name) {
  # Extract numeric vector of weights for this mixture
  weights <- mix_ratios[[mix_name]]

  # Assign names based on chemical order in scaled_df
  names(weights) <- unique(scaled_df$Chemicalname)

  scaled_df %>%
    group_by(Dose_uM) %>%
    summarise(
      IA_Response = 1 - prod(1 - weights[Chemicalname] * ScaledResp),
      .groups = "drop"
    ) %>%
    mutate(Mixture = mix_name)
}

# Run for all mixtures
mixture_results <- bind_rows(
  calc_mixture("EM"),
  calc_mixture("EXP1"),
  calc_mixture("EXP2"),
  calc_mixture("EXP3")
)

mixture_results

# Back Scale
global_min <- min(df$MaxResp, na.rm = TRUE)
global_max <- max(df$MaxResp, na.rm = TRUE)

IA_results <- mixture_results %>%
  mutate(
    IA_back = IA_Response * (global_max - global_min) + global_min
  )

# Cacluate IA BMC
IA_bmd_calc <- list()

for (mix in unique(IA_results$Mixture)) {

  # Subset data for this mixture
  mix_data <- IA_results %>% filter(Mixture == mix)

  # Prepare input row for concRespCore
  row <- list(
    conc = mix_data$Dose_uM,
    resp = mix_data$IA_back,  # use back-calculated response
    bmed = 0,                 # baseline
    cutoff = 0,
    onesd = 10,               # 10% increase = BMR
    assay = "BMC",
    name = mix
  )

  # Fit Hill model & compute BMC
  res <- concRespCore(
    row,
    fitmodels = c("hill"),
    conthits = TRUE,
    aicc = FALSE,
    bidirectional = FALSE,
    errfun = "dnorm",
    bmr_scale = 1
  )

  # Store selected results
  IA_bmd_calc[[mix]] <- res[, c("bmd", "bmdl", "bmdu", "name")]
}


#### Run Analysis ####
mix_cols <- c("EM_percent", "EXP1_percent", "EXP2_percent", "EXP3_percent")

CA_BMC <- calculate_ca_bmd(df_mixture, mix_cols) %>%
  mutate(method = "BMC CA") %>%
  rename(name = L1)

IA_BMC <- bind_rows(IA_bmd_calc)
IA_BMC$method <- "BMC IA"

combined_df <- rbind(CA_BMC, IA_BMC)

# Export results
write.csv(combined_df, "combined_bmc_IA_CA_Pred.csv", row.names = FALSE)

#### Plot ####
p1 <- ggplot(data = combined_df, aes(y = name, x = log10(bmd), color = name)) +
  geom_point(aes(shape = method), size = 3) +
  geom_errorbar(aes(xmin = log10(bmdl), xmax = log10(bmdu)), width = 0.1) +
  theme_bw() +
  labs(shape = "Model", color = "Mixture", x = "Log10 BMC (uM)", y = "Mixture")

p1
ggsave("bmc_plot_final_DAIA.tiff", plot = p1, device = "tiff",
       width = 5, height = 7, units = "in", dpi = 300)
