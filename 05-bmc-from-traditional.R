################################################
# Calculate BMCs of models using trad methods
# Written By: Kristin Eccles
# Date: July 29th, 2025
# Note: needs 01-traditional-individual-chem.R
# and 02-traditional-mixture-model.R to run
#################################################

# Load libraries
library(drc)
library(viridis)
library(broom)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(ggh4x)
library(ggpubr)
library(ggpmisc)
library(tidyverse)
library(reshape2)
library(tcplfit2)
library(tibble)

FIXED_C = 0

set.seed (4540)
#################################################
#### Import Data ####

file_paths <- list.files(path = "mix_pred_boot", pattern = "*.csv", full.names = TRUE)
data_list <- map(file_paths, read.csv)
data_unlist <- na.omit(as.data.frame(bind_rows(data_list)))
data_unlist[data_unlist == Inf | data_unlist == -Inf] <- NA
data_unlist <- na.omit(data_unlist)

# Create the `id` columns
data_unlist$id <- paste(data_unlist$iteration, data_unlist$group, data_unlist$mix, data_unlist$method)
data_unlist$id2 <- paste(data_unlist$mix, data_unlist$method)

unique(data_unlist$id2)
model_bmd_calc <- list()

for (chem in unique(data_unlist$id)) {

  # Subset the data for this chemical
  chem_data <- (data_unlist[data_unlist$id == chem, ])
  #chem_data <- subset(data_unlist, id %in% chem)

  # Define BMR setup
  row <- list(
    conc = chem_data$x,
    resp = chem_data$y,
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

model_bmd_calc

# Extract BMD summary info from each chemical's result
model_melt <- reshape2::melt(model_bmd_calc)

# Use dcast to reshape
tidy_model<- reshape2::dcast(model_melt, L1 ~ variable, value.var = "value")
tidy_model <- tidy_model %>%
  separate(L1, into = c("itter", "group", "mix", "method"), sep = " ")
tidy_model$id2 <- paste(tidy_model$group, tidy_model$mix, tidy_model$method)

model_summary <- na.omit(tidy_model) %>%
  group_by(id2) %>%
  summarise(
    bmdl  = quantile(bmd, 0.025),
    bmd_mean = quantile(bmd, 0.50),
    bmdu = quantile(bmd, 0.975),
    ac10_p5  = quantile(ac10, 0.025),
    ac10_p50 = quantile(ac10, 0.50),
    ac10_p95 = quantile(ac10, 0.975)
  )

model_bmd_coeff_edit<- model_summary %>%
  separate(id2, into = c("group", "mixture","model"), sep = " ", remove = FALSE)%>%
  as.data.frame()

model_bmd_coeff_edit$id <- paste(model_bmd_coeff_edit$mixture, model_bmd_coeff_edit$model)

# Create the plot
ggplot(tidy_model, aes(x = ac10, y = bmd, color = tp)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "left", label.y.npc = 0.95) +
  labs(title = "AC10 vs BMD10", x = "AC10",y = "BMD") +
  theme_minimal()

###########################################################################################################
#### Plot ####
model_bmd_coeff_edit$id <- "Traditional"
model_bmd_coeff_edit$L1 <- model_bmd_coeff_edit$mixture

combined_model <- combined_model %>%
  separate(method, into = c("method", "model"), sep = " ")
combined_model$id <- "BMC"

measured$id <- "Measured"


p2 <- ggplot(data = model_bmd_coeff_edit, aes(y = mixture, x = log10(bmd_mean), color = mixture)) +
  geom_point(aes(shape = model), size = 3) +
  geom_errorbar(aes(xmin = log10(bmdl), xmax = log10(bmdu)), width = 0.1) +
  theme_bw() +
  labs(shape = "Model",color = "Mixture", x = "Log10 BMC (uM)", y = "Mixture")
p2

# combined
p3 <- ggplot() +
  #TRADITIONAL
  geom_point(data = model_bmd_coeff_edit, aes(y = id, x = log10(bmd_mean), color = model), size = 2) +
  geom_errorbar(data = model_bmd_coeff_edit, aes(
    y = id,
    xmin = log10(bmdl), xmax = log10(bmdu),  color = model), width = 0.1) +

  #BMC
  geom_point(data = combined_model, aes(y = id, x = log10(bmd), color = model), size = 2) +
  geom_errorbar(data = combined_model, aes(
    y = id,
    xmin = log10(bmdl), xmax = log10(bmdu), color = model), width = 0.1) +

  #MEASURED
  geom_point(data = measured, aes(y = id, x = log10(bmd), color = id), size = 2) +
  geom_errorbar(data = measured, aes(
    y = id, xmin = log10(bmdl), xmax = log10(bmdu), color = id), width = 0.1) +

  theme_bw() +
  facet_grid(.~L1, drop = TRUE)+
  scale_color_manual(name = "Method",
                     values =  c("CA" = "#FDE725FF", "GCA" = "#7AD151FF", "IA" = "#2A788EFF", "Measured" = "#5A5A5A"),
                     labels = c("CA", "GCA", "IA", "Measured"))+
  labs(x = "Log10 BMC (uM)", y = "Mixture")
p3

ggsave("plot_compare_all.tiff", plot = p3, device = "tiff",
       width = 10, height = 3, units = "in", dpi = 300)
