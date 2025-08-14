################################################
# Mixtures Modeling for invitro PACs
# Written By: Kristin Eccles
# CAte: September 22nd, 2022
# Updated: November 8th, 2024
# Note: Needs 01-individual_chemical_dr.R to run
# must have created folder in project called "mix_pred_boot"
#################################################

# load libraries
library(MCMCglmm)

# Set Up
NUM_PTS <- 1000
MIN_LOGX <- (-5)
MAX_LOGX <- 5

MAX_Y <- 1
yVec <- vector(mode="numeric", length = 1000 - 1)
y <- as.data.frame(MAX_Y * 1:1000 / 1000)
colnames(y) <-"y"

#set up parameters for CI generation
MCiter <- 1000

#set seed for reproducibility
set.seed(8789)


#----------------------------------------------------------------------
#---------- Refit Curve for CA and IA  ---------
#----------------------------------------------------------------------
#need to have the same top and bottom
#only fit for reduced list
df$include <- df$Chemicalname %in%  individual_coeff_final$curve
df_reduced <- subset(df, include == TRUE)

#### Fit Parameters with Selected Model Hill Model ####
# fit for mixtures modeling
individual_model2<- drm(MaxResp~Dose_uM, data=df_reduced, curveid = Chemicalname,
                        type = "continuous",
                        lowerl = c(-Inf, 0),
                        upperl = c(0, Inf),
                        fct=LL.4(fixed=c(NA, FIXED_C , FIXED_D, NA),
                                 names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(individual_model2)
# Quick Plot of curve fits
plot(individual_model2)

# Extract coefficients
individual_model_wCI2 <- tidy(individual_model2, conf.int = TRUE)
individual_model_coeff2 <- individual_model_wCI2 %>%
  dplyr::select(term, curve, estimate) %>%
  tidyr::pivot_wider(names_from = term, values_from = estimate)

# add SE
slope_SE2 <- as.data.frame(subset(individual_model_wCI2, term == "Slope")$std.error)
colnames(slope_SE2) <- "SE_slope"
EXP2_SE2 <- as.data.frame(subset(individual_model_wCI2, term == "ED50")$std.error)
colnames(EXP2_SE2) <- "SE_ED50"
UpperLimit_SE2 <- as.data.frame(subset(individual_model_wCI2, term == "Upper Limit")$std.error)
colnames(UpperLimit_SE2) <- "SE_UpperLimit"
individual_model_coeff2 <- cbind(individual_model_coeff2,slope_SE2, EXP2_SE2)
#write.csv(individual_model_coeff, "individual_model_coeff.csv", row.names = FALSE)

# Bootstrap for CI
coeff_final_list2 <- split(individual_model_coeff2, f = individual_model_coeff2$curve)
slope_boot <- sapply(names(coeff_final_list2), function(x) rtruncnorm(n = MCiter, a = -Inf, b = 0, mean = coeff_final_list2[[x]]$Slope, sd = coeff_final_list2[[x]]$SE_slope*sqrt(2)))
EXP2_boot <- sapply(names(coeff_final_list2), function(x) rtruncnorm(n = MCiter, a = 0, b = Inf, mean = coeff_final_list2[[x]]$ED50, sd = coeff_final_list2[[x]]$SE_ED50*sqrt(2)))

# Add Mixing Ratios
MIX_FRACTIONS <- na.omit(read.csv("mixing_fractions.csv"))
individual_model_coeff2 <- left_join(individual_model_coeff2, MIX_FRACTIONS, by = c("curve" = "Chemical"))
individual_model_coeff2 <- na.omit(individual_model_coeff2)

tops <- left_join(individual_model_coeff, MIX_FRACTIONS, by = c("curve" = "Chemical"))
tops <- na.omit(tops)

#### Bootstrap for CI ####

coeff_final_list2 <- split(individual_model_coeff2, f=individual_model_coeff2$curve)
n <- names(coeff_final_list2)

slope_boot <- sapply(setNames(n, n), FUN = function(x) {
  rtruncnorm(n = MCiter, a = -Inf, b = 0, mean = coeff_final_list2[[x]]$Slope, sd = coeff_final_list2[[x]]$SE_slope*sqrt(2))},
  simplify = FALSE,USE.NAMES = TRUE)
slope_boot_melt <- reshape2::melt(slope_boot)

EXP2_boot <- sapply(setNames(n, n), FUN = function(x) {
  rtruncnorm(n = MCiter, a = 0, b = Inf, mean = coeff_final_list2[[x]]$ED50, sd = coeff_final_list2[[x]]$SE_ED50*sqrt(2))},
  simplify = FALSE,USE.NAMES = TRUE)

EXP2_boot_melt <- reshape2::melt(EXP2_boot)

bootmat <- as.data.frame(cbind(EXP2_boot_melt[,1], slope_boot_melt))
colnames(bootmat) <- cbind("ED50","slope", "chemical")
#add iteration number for each chemical
bootmat$itter <- 1:MCiter

#Mixing Ratios
MIX_FRACTIONS <- na.omit(read.csv("mixing_fractions.csv"))

#check
sum(individual_model_coeff2$EM_percent)
sum(individual_model_coeff2$EXP1_percent)
sum(individual_model_coeff2$EXP2_percent)
sum(individual_model_coeff2$EXP3_percent)

EXP3_top <- weighted.mean(tops$`Upper Limit`, tops$EXP3_percent, na.rm = TRUE)
EM_top <- weighted.mean(tops$`Upper Limit`, tops$EM_percent, na.rm = TRUE)
EXP1_top <- weighted.mean(tops$`Upper Limit`, tops$EXP1_percent, na.rm = TRUE)
EXP2_top <- weighted.mean(tops$`Upper Limit`, tops$EXP2_percent, na.rm = TRUE)

bootmat <- left_join(bootmat, individual_model_coeff2[,c("curve","EM_percent","EXP1_percent",
                                                         "EXP2_percent", "EXP3_percent")],
                     by= c("chemical" = "curve"), keep=FALSE )
#make new dataframe for CI generation
bootmat_list <- split(bootmat, f=bootmat$itter)

# Refit with slope of -1 for GCA
individual_model_b1<- drm(MaxResp~Dose_uM, data=df_reduced, curveid= Chemicalname,
                          type = "continuous",
                          lowerl = c(0, 0),
                          upperl = c(100, Inf),
                          fct=LL.4(fixed=c(FIXED_B, FIXED_C , NA, NA),
                                   names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

# Quick Plot of curve fits
plot(individual_model_b1)

# get coefficients
B1individual_model_wCI <- tidy(individual_model_b1, conf.int = TRUE)
B1individual_model_wCI

#reorganize
B1individual_model_coeff <- B1individual_model_wCI%>%
  dplyr::select(term, curve, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate)%>%
  as.data.frame
B1individual_model_coeff

# add std. error
B1EXP2_SE <- as.data.frame(subset(B1individual_model_wCI, term == "ED50")$std.error)
colnames(B1EXP2_SE) <- "SE_ED50"

B1UpperLimit_SE <- as.data.frame(subset(B1individual_model_wCI, term == "Upper Limit")$std.error)
colnames(B1UpperLimit_SE) <- "SE_UpperLimit"

B1individual_model_coeff <- cbind(B1individual_model_coeff,B1EXP2_SE, B1UpperLimit_SE)

#Bootstrap for CIs
B1coeff_final_list <- split(B1individual_model_coeff, f=B1individual_model_coeff$curve)
n <- names(B1coeff_final_list)

B1Top_boot <- sapply(setNames(n, n), FUN = function(x) {
  rtnorm(n = MCiter, lower = 0, upper = Inf, mean = B1coeff_final_list[[x]]$`Upper Limit`, sd = B1coeff_final_list[[x]]$SE_UpperLimit*sqrt(2))},
  simplify = FALSE,USE.NAMES = TRUE)
B1Top_boot_melt <- reshape2::melt(B1Top_boot)

B1EXP2_boot <- sapply(setNames(n, n), FUN = function(x) {
  rtnorm(n = MCiter, lower = 0, upper = Inf, mean =  B1coeff_final_list[[x]]$ED50, sd = B1coeff_final_list[[x]]$SE_ED50*sqrt(2))},
  simplify = FALSE,USE.NAMES = TRUE)

B1EXP2_boot_melt <- reshape2::melt(B1EXP2_boot)

B1bootmat <- as.data.frame(cbind(B1EXP2_boot_melt[,1], B1Top_boot_melt))
colnames(B1bootmat) <- cbind("ED50","Top", "chemical")
#add iteration number for each chemical
B1bootmat$itter <- 1:MCiter

#add mixing ratios
B1bootmat <- na.omit(left_join(B1bootmat, individual_model_coeff2[,c("curve","EM_percent","EXP1_percent",
                                                             "EXP2_percent", "EXP3_percent")],
                       by= c("chemical" = "curve"), keep=FALSE ))
#make new dataframe for CI generation
B1bootmat_list <- split(B1bootmat, f=B1bootmat$itter)

#----------------------------------------------------------------------
#---------- Predict Mixture Responses  ---------
#----------------------------------------------------------------------

#### Independent Action ####
IA_EM_CI <-
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      1-prod((1/(1+ (exp(-bootmat_list[[j]]$slope * (log(x*bootmat_list[[j]]$EM_percent) - (bootmat_list[[j]]$ED50)))))))
    })})

unlist_EM_CI <- cbind(x,  reshape2::melt(IA_EM_CI))
colnames(unlist_EM_CI) <- c("x", "y", "iteration")
unlist_EM_CI$mix <- "EM"
unlist_EM_CI$method <- "IA"

unlist_EM_CI<- unlist_EM_CI%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EM_top))))

write.csv(unlist_EM_CI, "mix_pred_boot/unlist_EM_CI.csv", row.names = FALSE)

IA_CI_EM_final <- unlist_EM_CI%>%
  group_by(x)%>%
  summarise(mean = quantile(y, 0.5, na.rm =TRUE),
            y_lower = quantile(y,0.025, na.rm =TRUE),
            y_upper = quantile(y, 0.975, na.rm =TRUE))%>%
  mutate(across(c(mean, y_lower, y_upper), ~ scales::rescale(., to = c(0, EM_top))))%>%
  as.data.frame()
IA_CI_EM_final$mix_ratio <- "EM"

IA_EXP1_CI <-
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      1-prod((1/(1+ (exp(-bootmat_list[[j]]$slope * (log(x*bootmat_list[[j]]$EXP1_percent) - (bootmat_list[[j]]$ED50)))))))
    })})

unlist_EXP1_CI <- cbind(x,  reshape2::melt(IA_EXP1_CI))
colnames(unlist_EXP1_CI) <- c("x", "y", "iteration")
unlist_EXP1_CI$mix <- "EXP1"
unlist_EXP1_CI$method <- "IA"

unlist_EXP1_CI<- unlist_EXP1_CI%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EXP1_top))))

write.csv(unlist_EXP1_CI, "mix_pred_boot/unlist_EXP1_CI.csv", row.names = FALSE)

IA_CI_EXP1_final <- unlist_EXP1_CI%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            y_lower = quantile(y,0.025, na.rm =TRUE),
            y_upper = quantile(y, 0.975, na.rm =TRUE))%>%
  mutate(across(c(mean, y_lower, y_upper), ~ scales::rescale(., to = c(0, EXP1_top))))%>%
  as.data.frame()
IA_CI_EXP1_final$mix_ratio <- "EXP1"

IA_EXP2_CI <-
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      1-prod((1/(1+ (exp(-bootmat_list[[j]]$slope * (log(x*bootmat_list[[j]]$EXP2_percent) - (bootmat_list[[j]]$ED50)))))))
    })})

unlist_EXP2_CI <- cbind(x,  reshape2::melt(IA_EXP2_CI))
colnames(unlist_EXP2_CI) <- c("x", "y", "iteration")
unlist_EXP2_CI$mix <- "EXP2"
unlist_EXP2_CI$method <- "IA"

unlist_EXP2_CI<- unlist_EXP2_CI%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EXP2_top))))
write.csv(unlist_EXP2_CI, "mix_pred_boot/unlist_EXP2_CI.csv", row.names = FALSE)

IA_CI_EXP2_final <- unlist_EXP2_CI%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            y_lower = quantile(y,0.025, na.rm =TRUE),
            y_upper = quantile(y, 0.975, na.rm =TRUE))%>%
  mutate(across(c(mean, y_lower, y_upper), ~ scales::rescale(., to = c(0, EXP2_top))))%>%
  as.data.frame()
IA_CI_EXP2_final$mix_ratio <- "EXP2"

IA_EXP3_CI <-
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      1-prod((1/(1+ (exp(-bootmat_list[[j]]$slope * (log(x*bootmat_list[[j]]$EXP3_percent) - (bootmat_list[[j]]$ED50)))))))
    })})

unlist_EXP3_CI <- cbind(x,  reshape2::melt(IA_EXP3_CI))
colnames(unlist_EXP3_CI) <- c("x", "y", "iteration")
unlist_EXP3_CI$mix <- "EXP3"
unlist_EXP3_CI$method <- "IA"

unlist_EXP3_CI<- unlist_EXP3_CI%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EXP3_top))))
write.csv(unlist_EXP3_CI, "mix_pred_boot/unlist_EXP3_CI.csv", row.names = FALSE)

IA_CI_EXP3_final <- unlist_EXP3_CI%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            y_lower = quantile(y,0.025, na.rm =TRUE),
            y_upper = quantile(y, 0.975, na.rm =TRUE))%>%
  mutate(across(c(mean, y_lower, y_upper), ~ scales::rescale(., to = c(0, EXP3_top))))%>%
  as.data.frame()
IA_CI_EXP3_final$mix_ratio <- "EXP3"

# Combined data frames
IA_df <- rbind(IA_CI_EM_final, IA_CI_EXP1_final, IA_CI_EXP2_final, IA_CI_EXP3_final)
IA_df$Method <- "IA"
write.csv(IA_df, "IA_df.csv", row.names = FALSE)

#plot
p1adj <- ggplot()+
  geom_line(data = IA_df, aes(x= log10(x), y = mean, color = mix_ratio), linewidth = 1)+
  geom_ribbon(data=IA_df, aes(x=log10(x), y= mean, ymin=y_lower, ymax=y_upper, fill = mix_ratio), alpha=0.2) +
  theme_bw()+
  labs(x="Log10 Concentration (uM)", y="% Max MeBio Response", color = "Mixing Ratio \nIndependent Action",
       fill = "Mixing Ratio \nIndependent Action") +
  scale_fill_viridis(discrete= TRUE)+
  scale_color_viridis(discrete= TRUE)+
  theme_bw() +
  theme(legend.title = element_text(face="bold", size=11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size = 16, hjust=0.5))
p1adj

#### Concentration Addition ####
# x = (e^(e *b)* (d/y - 1))^(1/b)

CA_EM_CI <-
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(y, MARGIN = 1, FUN = function(y) {
      sum(bootmat_list[[j]]$EM_percent*(bootmat_list[[j]]$ED50 * exp((bootmat_list[[j]]$slope)*log((1-y)/y))))
    })})

unlist_CA_EM_CI <- cbind(y,  reshape2::melt(CA_EM_CI))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EM_top))))

colnames(unlist_CA_EM_CI) <- c("y", "x", "iteration")
unlist_CA_EM_CI$mix <- "EM"
unlist_CA_EM_CI$method <- "CA"
write.csv(unlist_CA_EM_CI, "mix_pred_boot/unlist_CA_EM_CI.csv", row.names = FALSE)

CA_CI_EM_final <- unlist_CA_EM_CI%>%
  group_by(y)%>%
  summarise(mean = quantile(x,0.5, na.rm =TRUE),
            x_lower = quantile(x,0.025, na.rm =TRUE),
            x_upper = quantile(x, 0.975, na.rm =TRUE))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EM_top))))%>%
  as.data.frame()
CA_CI_EM_final$mix_ratio <- "EM"

CA_EXP1_CI <-
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(y, MARGIN = 1, FUN = function(y) {
      sum(bootmat_list[[j]]$EXP1_percent*(bootmat_list[[j]]$ED50 * exp((bootmat_list[[j]]$slope)*log((1-y)/y))))
    })})

unlist_CA_EXP1_CI <- cbind(y,  reshape2::melt(CA_EXP1_CI))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EXP1_top))))

colnames(unlist_CA_EXP1_CI) <- c("y", "x", "iteration")
unlist_CA_EXP1_CI$mix <- "EXP1"
unlist_CA_EXP1_CI$method <- "CA"
write.csv(unlist_CA_EXP1_CI, "mix_pred_boot/unlist_CA_EXP1_CI.csv", row.names = FALSE)

CA_CI_EXP1_final <- unlist_CA_EXP1_CI%>%
  group_by(y)%>%
  summarise(mean = quantile(x,0.5, na.rm =TRUE),
            x_lower = quantile(x,0.025, na.rm =TRUE),
            x_upper = quantile(x, 0.975, na.rm =TRUE))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EXP1_top))))%>%
  as.data.frame()
CA_CI_EXP1_final$mix_ratio <- "EXP1"

CA_EXP2_CI <-
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(y, MARGIN = 1, FUN = function(y) {
      sum(bootmat_list[[j]]$EXP2_percent*(bootmat_list[[j]]$ED50 * exp((bootmat_list[[j]]$slope )*log((1-y)/y))))
    })})

unlist_CA_EXP2_CI <- cbind(y,  reshape2::melt(CA_EXP2_CI))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EXP2_top))))

colnames(unlist_CA_EXP2_CI) <- c("y", "x", "iteration")
unlist_CA_EXP2_CI$mix <- "EXP2"
unlist_CA_EXP2_CI$method <- "CA"
write.csv(unlist_CA_EXP2_CI, "mix_pred_boot/unlist_CA_EXP2_CI.csv", row.names = FALSE)

CA_CI_EXP2_final <- unlist_CA_EXP2_CI%>%
  group_by(y)%>%
  summarise(mean = quantile(x,0.5, na.rm =TRUE),
            x_lower = quantile(x,0.025, na.rm =TRUE),
            x_upper = quantile(x, 0.975, na.rm =TRUE))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EXP2_top))))%>%
  as.data.frame()
CA_CI_EXP2_final$mix_ratio <- "EXP2"

CA_EXP3_CI <-
  lapply(1:length(bootmat_list), FUN =  function(j) {
    apply(y, MARGIN = 1, FUN = function(y) {
      sum(bootmat_list[[j]]$EXP3_percent*(bootmat_list[[j]]$ED50 * exp((bootmat_list[[j]]$slope )*log((1-y)/y))))
    })})

unlist_CA_EXP3_CI <- cbind(y,  reshape2::melt(CA_EXP3_CI))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EXP3_top))))

colnames(unlist_CA_EXP3_CI) <- c("y", "x", "iteration")
unlist_CA_EXP3_CI$mix <- "EXP3"
unlist_CA_EXP3_CI$method <- "CA"
write.csv(unlist_CA_EXP3_CI, "mix_pred_boot/unlist_CA_EXP3_CI.csv", row.names = FALSE)

CA_CI_EXP3_final <- unlist_CA_EXP3_CI%>%
  group_by(y)%>%
  summarise(mean = quantile(x,0.5, na.rm =TRUE),
            x_lower = quantile(x,0.025, na.rm =TRUE),
            x_upper = quantile(x, 0.975, na.rm =TRUE))%>%
  mutate(across(c(y), ~ scales::rescale(., to = c(0, EXP3_top))))%>%
  as.data.frame()
CA_CI_EXP3_final$mix_ratio <- "EXP3"

# Combined dataframes
CA_df <- rbind(CA_CI_EM_final, CA_CI_EXP1_final, CA_CI_EXP2_final, CA_CI_EXP3_final)
CA_df$Method <- "CA"
write.csv(CA_df, "CA_df.csv", row.names = FALSE)

#plot
p2adj <- ggplot()+
  geom_line(data = CA_df, aes(x= log10(mean), y = y, color = mix_ratio), linewidth = 1)+
  geom_ribbon(data = CA_df, aes(x=log10(mean), y = y, xmin=log10(x_lower), xmax=log10(x_upper), fill = mix_ratio), alpha=0.2) +
  theme_bw()+
  labs(x="Log10 Concentration (uM)", y="% Max MeBio Response", color = "Mixing Ratio \nConcentration Addition",
       fill = "Mixing Ratio \nConcentration Addition") +
  theme_bw() +
  scale_fill_viridis(discrete= TRUE)+
  scale_color_viridis(discrete= TRUE)+
  theme(legend.title = element_text(face="bold", size=11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size = 16, hjust=0.5))
p2adj

#### Generalized Concentration Addition ####
# GCA Predictions

GCA_EM_CI <-
  lapply(1:length(B1bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      sum((B1bootmat_list[[j]]$Top)*((x*B1bootmat_list[[j]]$EM_percent)/B1bootmat_list[[j]]$ED50))/
        (1+sum((x*B1bootmat_list[[j]]$EM_percent)/B1bootmat_list[[j]]$ED50))
    })})

unlist_GCA_EM_CI <- cbind(x,  reshape2::melt(GCA_EM_CI))
colnames(unlist_GCA_EM_CI) <- c("x", "y", "iteration")
unlist_GCA_EM_CI$mix <- "EM"
unlist_GCA_EM_CI$method <- "GCA"
write.csv(unlist_GCA_EM_CI, "mix_pred_boot/unlist_GCA_EM_CI.csv", row.names = FALSE)


GCA_CI_EM_final <- unlist_GCA_EM_CI%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            lower = quantile(y,0.025, na.rm =TRUE),
            upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()
GCA_CI_EM_final$mix_ratio <- "EM"

GCA_EXP1_CI <-
  lapply(1:length(B1bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      sum((B1bootmat_list[[j]]$Top)*((x*B1bootmat_list[[j]]$EXP1_percent)/B1bootmat_list[[j]]$ED50))/
        (1+sum((x*B1bootmat_list[[j]]$EXP1_percent)/B1bootmat_list[[j]]$ED50))
    })})

unlist_GCA_EXP1_CI <- cbind(x,  reshape2::melt(GCA_EXP1_CI))
colnames(unlist_GCA_EXP1_CI) <- c("x", "y", "iteration")
unlist_GCA_EXP1_CI$mix <- "EXP1"
unlist_GCA_EXP1_CI$method <- "GCA"
write.csv(unlist_GCA_EXP1_CI, "mix_pred_boot/unlist_GCA_EXP1_CI.csv", row.names = FALSE)

GCA_CI_EXP1_final <- unlist_GCA_EXP1_CI%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            lower = quantile(y,0.025, na.rm =TRUE),
            upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()
GCA_CI_EXP1_final$mix_ratio <- "EXP1"

GCA_EXP2_CI <-
  lapply(1:length(B1bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      sum((B1bootmat_list[[j]]$Top)*((x*B1bootmat_list[[j]]$EXP2_percent)/B1bootmat_list[[j]]$ED50))/
        (1+sum((x*B1bootmat_list[[j]]$EXP2_percent)/B1bootmat_list[[j]]$ED50))
    })})

unlist_GCA_EXP2_CI <- cbind(x,  reshape2::melt(GCA_EXP2_CI))
colnames(unlist_GCA_EXP2_CI) <- c("x", "y", "iteration")
unlist_GCA_EXP2_CI$mix <- "EXP2"
unlist_GCA_EXP2_CI$method <- "GCA"
write.csv(unlist_GCA_EXP2_CI, "mix_pred_boot/unlist_GCA_EXP2_CI.csv", row.names = FALSE)

GCA_CI_EXP2_final <- unlist_GCA_EXP2_CI%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            lower = quantile(y,0.025, na.rm =TRUE),
            upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()
GCA_CI_EXP2_final$mix_ratio <- "EXP2"

GCA_EXP3_CI <-
  lapply(1:length(B1bootmat_list), FUN =  function(j) {
    apply(x, MARGIN = 1, FUN = function(x) {
      sum((B1bootmat_list[[j]]$Top)*((x*B1bootmat_list[[j]]$EXP3_percent)/B1bootmat_list[[j]]$ED50))/
        (1+sum((x*B1bootmat_list[[j]]$EXP3_percent)/B1bootmat_list[[j]]$ED50))
    })})

unlist_GCA_EXP3_CI <- cbind(x,  reshape2::melt(GCA_EXP3_CI))
colnames(unlist_GCA_EXP3_CI) <- c("x", "y", "iteration")
unlist_GCA_EXP3_CI$mix <- "EXP3"
unlist_GCA_EXP3_CI$method <- "GCA"
write.csv(unlist_GCA_EXP3_CI, "mix_pred_boot/unlist_GCA_EXP3_CI.csv", row.names = FALSE)

GCA_CI_EXP3_final <- unlist_GCA_EXP3_CI%>%
  group_by(x)%>%
  summarise(mean = quantile(y,0.5, na.rm =TRUE),
            lower = quantile(y,0.025, na.rm =TRUE),
            upper = quantile(y, 0.975, na.rm =TRUE))%>%
  as.data.frame()
GCA_CI_EXP3_final$mix_ratio <- "EXP3"

# Combined dataframes
GCA_df <- rbind(GCA_CI_EM_final, GCA_CI_EXP1_final, GCA_CI_EXP2_final, GCA_CI_EXP3_final)
GCA_df$Method <- "GCA"

GCA_df$y <- GCA_df$mean
GCA_df$y_lower <- GCA_df$lower
GCA_df$y_upper <- GCA_df$upper

write.csv(GCA_df, "GCA_df.csv", row.names = FALSE)

#plot
p3adj <- ggplot()+
  geom_line(data = GCA_df, aes(x= log10(x), y = y, color = mix_ratio), linewidth = 1)+
  geom_ribbon(data=GCA_df, aes(x=log10(x), y=y, ymin=y_lower, ymax=y_upper, fill = mix_ratio), alpha=0.2) +
  theme_bw()+
  labs(x="Log10 Concentration (uM)", y="% Max MeBio Response", color = "Mixing Ratio \nGCA",
       fill = "Mixing Ratio \nGCA") +
  scale_fill_viridis(discrete= TRUE)+
  scale_color_viridis(discrete= TRUE)+
  theme_bw() +
  theme(legend.title = element_text(face="bold", size=11),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size = 16, hjust=0.5))
p3adj

