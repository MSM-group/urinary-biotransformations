#Clear environment
rm(list = ls())
#Load required packages
library(pacman)
pacman::p_load("tidyverse", "ggpubr", "ggpattern", "RColorBrewer")
#Load biotransformation binary probability matrix
dat_in <- readr::read_csv("data/ec_biotransformation_probability_matrix.csv")
#Load EC class metadata
ec <- readxl::read_excel("data/metadata/ec_representatives_metadata.xlsx") %>%
  dplyr::select(ec_class, descsription) %>%
  dplyr::distinct() %>%
  mutate(ec_class = as.factor(ec_class))
#Create data frame for plotting
dat <- dat_in %>%
  tidyr::pivot_longer(-x2nd_lvl, names_to = "compound", values_to = "activity") %>%
  dplyr::mutate(x2nd_lvl = forcats::as_factor(x2nd_lvl) %>%
                  forcats::fct_rev() %>%
                  forcats::fct_relevel("2.4" , after = 5),
                activity = forcats::as_factor(activity) %>%
                  forcats::fct_relevel("envipath", "drugbug", "consensus", "0")) %>%
  dplyr::select(x2nd_lvl, compound, activity)
#Make plot
plot <- ggplot2::ggplot(dat, aes(x = compound, y = x2nd_lvl, fill = activity)) +
  ggplot2::geom_tile(show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = c(RColorBrewer::brewer.pal(3, "Set2")[1:2], RColorBrewer::brewer.pal(3, "Set2")[1], "grey"), labels = c("enviPath", "DrugBug", "consensus", "no prediction")) +
  ggpattern::geom_tile_pattern(aes(x = compound, y = x2nd_lvl, fill = activity),
                               dplyr::filter(dat, activity == "consensus"),
                               color = NA,
                               linewidth = 0.05,
                               pattern_fill = RColorBrewer::brewer.pal(3, "Set2")[2],
                               pattern_linetype = "blank",
                               pattern_density = 0.5,
                               height = 1,
                               width = 1) +
  ggpubr::theme_pubr() +
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
                 axis.text.y = element_text(size = 12),
                 axis.title.y = element_text(size = 12),
                 axis.title.x = element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size = 6),
                 legend.title = element_blank()) +
  labs(x = "Compound", y = "EC class")
#Save plot
ggsave("plots/fig3intermediate.png",
       plot,
       units = "cm",
       device = "png",
       width = 9.4,
       height = 9.4)
