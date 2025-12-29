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
library(scales)

#################################################
#### Contribution Plot ####

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
pie_p1 <- ggplot(mixture_long, aes(x = Mixture, y = Percent, fill = Chemical)) +
  geom_bar(stat = "identity", width = .7) +  # Use your Percent column
  scale_fill_viridis_d(option = "D") +
  scale_y_continuous(labels = percent_format(scale = 100)) +  # Show % on y-axis
  theme_minimal() +
  labs(
    y = "Percent Contribution",
    fill = "Chemical"
  )
pie_p1
ggsave("contribution_plot.tiff", pie_p1, dpi= 600)

# Combined Plot:

combined_plot_all <- ggarrange(individual_plot, # from 01
                               indiv_bmc, # from 03
                               pie_p1, #above
                               ncol = 3,
                               vjust =3,
                               labels = "AUTO",
                               common.legend = TRUE,
                               legend = "bottom")
combined_plot_all
ggsave("Background_fig.jpg", combined_plot_all,  height =5, width =10)
