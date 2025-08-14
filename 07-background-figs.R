################################################
# Plots for Mixing Ratios
# Written By: Kristin Eccles
# Date: August 14th, 2025
#################################################

# load libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)

# load data
MIX_FRACTIONS <- na.omit(read.csv("mixing_fractions.csv"))

#################################################
#### Pie Plot ####

# Reshape data to long format
mixture_long <- pivot_longer(
  MIX_FRACTIONS,
  cols = ends_with("_percent"),
  names_to = "Mixture",
  values_to = "Percent"
)

# Clean up Mixture names
mixture_long$Mixture <- gsub("_percent", "", mixture_long$Mixture)

# Create pie chart using ggplot
pie_p1 <- ggplot(mixture_long, aes(x = "", y = Percent, fill = Chemical)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~Mixture) +
  scale_fill_viridis_d(option = "D") +
  theme_void() +
  labs(
    fill = "Chemical"
  )
pie_p1
ggsave("pie_chart.tiff", pie_p1, dpi= 600)
