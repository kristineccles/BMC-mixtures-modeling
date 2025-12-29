################################################
# 01 Mixture Modeling BMCs - individual chems
# Written By: Kristin Eccles
# Date: June 11, 2025
# Updated July 28th, 2025
#################################################
# load libraries
library(drc)
library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(sjPlot)
library(cowplot)
library(data.table)
library(reshape2)
library(ggpubr)
library(broom)
library(viridis)
library(truncnorm)

# Set Up
#set up the data frame
NUM_PTS <- 100
MIN_LOGX <- (-3)
MAX_LOGX <- 15

xVec <- 1:NUM_PTS
x <- as.data.frame(10 ^ (MIN_LOGX + ((MAX_LOGX-MIN_LOGX)*xVec/NUM_PTS)))
colnames(x) <- "x"

#fixed parameters
FIXED_C = 0 #lower limit
FIXED_B = -1 # Slope - increasing slopes are negative in drm
FIXED_D <- 100

#----------------------------------------------------------------------
#---------- Curve Fitting for Individual Chemicals  ---------
#----------------------------------------------------------------------

#### Fit Parameters with Selected Model Hill Model ####
# fit for mixtures modeling
individual_model<- drm(MaxResp~Dose_uM, data=df, curveid = Chemicalname,
                       type = "continuous",
                       lowerl = c(-Inf, 0, 0),
                       upperl = c(0, 100, Inf),
                       fct=LL.4(fixed=c(NA, FIXED_C , NA, NA),
                                names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(individual_model)
# Quick Plot of curve fits
plot(individual_model, col = TRUE)

ED(individual_model, c(10))
ED(individual_model, c(50))

# get coefficients
individual_model_wCI <- tidy(individual_model, conf.int = TRUE)
individual_model_wCI

#reorganize
individual_model_coeff <- individual_model_wCI%>%
  dplyr::select(term, curve, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate)%>%
  as.data.frame
individual_model_coeff

# add SE
slope_SE <- as.data.frame(subset(individual_model_wCI, term == "Slope")$std.error)
colnames(slope_SE) <- "SE_slope"
ED50_SE <- as.data.frame(subset(individual_model_wCI, term == "ED50")$std.error)
colnames(ED50_SE) <- "SE_ED50"
UpperLimit_SE <- as.data.frame(subset(individual_model_wCI, term == "Upper Limit")$std.error)
colnames(UpperLimit_SE) <- "SE_UpperLimit"
individual_model_coeff <- cbind(individual_model_coeff,slope_SE, ED50_SE, UpperLimit_SE)
#write.csv(individual_model_coeff, "individual_model_coeff.csv", row.names = FALSE)

#----------------------------------------------------------------------
#---------- Curve Selection  ---------
#----------------------------------------------------------------------

#### Chemicals to include in Mixtures Modeling ####
# select chemicals for inclusion in mixture modeling
# remove chemicals where p-value cannot be estimated
individual_model_select <- subset(individual_model_wCI, p.value > 0)
slope_select <- subset(individual_model_select, term == "Slope")

# the 95% CI for slope included 0 (i.e., the direction of response could not be determined).
# check to see if signs are the same for the CI
slope_select$slope_ci <- sign(slope_select$conf.high)==sign(slope_select$conf.low)
criteria_chems <- subset(slope_select, slope_ci == "TRUE")$curve
criteria_chems

# chems to include
unique_chem <- unique(individual_model_coeff$curve)
match_chems <- as.data.frame(unique_chem)
match_chems$match <- unique_chem %in% criteria_chems

#Keep only matching
match_chems <- subset(match_chems, match ==TRUE)

# get reduced chem list from all 3 criteria
individual_model_coeff$include <- individual_model_coeff$curve %in%  match_chems$unique_chem
individual_coeff_final<- subset(individual_model_coeff, include == TRUE)
individual_coeff_final

ggplot(individual_coeff_final, aes(x = curve, y = log10(ED50), color = curve)) +
  geom_point(size = 3) +  # Dot for BMD
  geom_errorbar(aes(ymin = log10(ED50-SE_ED50*1.96), ymax = log10(ED50+SE_ED50*1.96)), width = 0.2) +  # Whiskers
  coord_flip() +
  theme_minimal() +
  labs(x = "Chemical",
    y = "Log10 EC50",
    color = "Chemical")

#----------------------------------------------------------------------
#---------- Plotting Curves  ---------
#----------------------------------------------------------------------

# ggplot figure
# set up x var for individual curve predictions
unique_chem <- unique(individual_model_coeff$curve)

# set up output
chem_out_df<- data.frame(x)

for(i in unique_chem) {

  #model selection
  chem_df <- subset(individual_model_coeff, curve == i)
  out <- common_top/(1+ (exp(chem_df$Slope * (log(x) - (chem_df$ED50)))))
  chem_out_df[i] <- out[1]

}
chem_out_df

chem_df_melt <-reshape2::melt(data= chem_out_df, id= "x")

individual_plot <- ggplot(data= chem_df_melt, aes(x = log10(x), y = value, color = variable))+
  geom_line(linewidth = 1)+
  theme_bw()+
  scale_color_viridis_d() +
  labs(x = "Log10 Concentration", y = "% Response", color = "Chemical")
individual_plot

ggsave(individual_plot, file="individual_curves.jpg", height = 5, width = 5)
