#Clear environment
rm(list = ls())
#Load required packages
library(pacman)
pacman::p_load("tidyverse", "ggpubr", "janitor", "ggnewscale", "pBrackets")
#Create function to import blast results
import_blast <- function(file, path){
  read_tsv(paste0(path, file), col_names = c("qseqid", "qlen", "saccver", "ssciname", "pident", "qcovs", "length", "evalue", "bitscore"), col_types = "ciccdiidd")
}
#Create tidy eval functions to parse data
#Function to evaluate whether or not enzyme is present in an organism according to blast results
test_presence = function(n){
  dplyr::case_when(n == 0 ~ FALSE,
                   TRUE ~ TRUE)
}
#Function to add text description of first-level EC classes
classify_enyzme <- function(class){
  dplyr::case_when(
    stringr::str_sub(class, start = 1, end = 1) == "1" ~ "Oxidoreductase",
    stringr::str_sub(class, start = 1, end = 1) == "2" ~ "Transferase",
    stringr::str_sub(class, start = 1, end = 1) == "3" ~ "Hydrolase"
  )
}
#Function to color plot by first-level EC class
color_plot <- function(presence, enzyme_class){
  dplyr::case_when(
    presence == FALSE ~ "absent",
    TRUE ~ enzyme_class
  )
}
#Import assembly metadata
metadat <- read_tsv("data/metadata/PRJNA316969_AssemblyDetails.txt", skip = 1) %>%
  janitor::clean_names()
#Import blast results
dir <- "data/thomas_white_blast/"
dat_in <- tibble(file = list.files(dir)) %>%
  dplyr::mutate(accession = word(file, sep = "_", start = 1, end = 2),
         blast = map(file, import_blast, path = dir)) %>%
  dplyr::select(-file) %>%
  tidyr::unnest(blast)
#Retain blast hits with a bitscore >50
dat <- dat_in %>%
  filter(bitscore > 50) %>%
  dplyr::mutate(qseqid = stringr::str_replace(qseqid, pattern = "1.1_NADP-dependent_oxidoreductase", "1.3_NADP-dependent_oxidoreductase")) #Correct EC class to match uniprot annotation
#Process blast data and metadata for plotting
presence_absence <- dat %>%
  dplyr::group_by(accession, qseqid) %>%
  dplyr::summarise(count = n(), bitscore = max(bitscore)) %>% #count blast hits per strain for each enzyme and retain max bitscore
  dplyr::ungroup() %>%
  tidyr::complete(accession, qseqid) %>% #add strain-enzyme combinations for which there were no blast hits
  tidyr::replace_na(list(count = 0)) %>%
  dplyr::mutate(present = test_presence(count), ec_class = stringr::word(qseqid, sep = "_", start = 1, end = 1), top_level_class = classify_enyzme(ec_class)) %>% #convert blast hit counts into binary presence-absence, read EC classes from sequence names
  dplyr::arrange(ec_class) %>%
  dplyr::mutate(top_level_class = forcats::fct_inorder(top_level_class)) %>%
  dplyr::left_join(metadat, by = c("accession" = "number_assembly")) %>% #add strain metadata including taxonomy
  dplyr::mutate(lab = stringr::str_c(taxonomy, strain, sep = " ")) %>%
  dplyr::arrange(taxonomy) %>%
  dplyr::mutate(lab = forcats::fct_inorder(lab)) %>%
  dplyr::mutate(qseqid = stringr::str_replace_all(qseqid, "_", " ")) %>%
  dplyr::arrange(ec_class) %>%
  dplyr::mutate(qseqid = forcats::fct_inorder(qseqid))
#Make plot
presence_absence_plot_bitscore_multicolor <- ggplot2::ggplot() +
  geom_tile(aes(x = strain, y = qseqid, fill = bitscore), filter(presence_absence, top_level_class == "Oxidoreductase")) +
  scale_fill_gradient(low = "firebrick1", high = "firebrick4", na.value = "white", limits = c(min(na.omit(presence_absence$bitscore)), max(na.omit(presence_absence$bitscore)))) +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = strain, y = qseqid, fill = bitscore), filter(presence_absence, top_level_class == "Transferase")) +
  scale_fill_gradient(low = "gold1", high = "gold4", na.value = "white", limits = c(min(na.omit(presence_absence$bitscore)), max(na.omit(presence_absence$bitscore)))) +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = strain, y = qseqid, fill = bitscore), filter(presence_absence, top_level_class == "Hydrolase")) +
  scale_fill_gradient(low = "dodgerblue1", high = "dodgerblue4", na.value = "white", limits = c(min(na.omit(presence_absence$bitscore)), max(na.omit(presence_absence$bitscore)))) +
  facet_grid(cols = vars(taxonomy), rows = vars(top_level_class), scales = "free", space = "free", switch = "both") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic", size = 6),
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic", size = 6),
        strip.text.y = element_text(size = 6),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        panel.spacing.x = unit(0.05, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.1, "cm"),
        legend.text = element_text(size = 6, angle = 90, hjust = 1),
        legend.title = element_text(size = 6, vjust = 1))
#Save plot and add brackets
jpeg(filename = "plots/fig4intermdiate.jpg",
     units = "cm",
     width = 18.8,
     height = 12,
     res = 300)
presence_absence_plot_bitscore_multicolor
grid.brackets(1050, 1020, 980, 1020, h = 0.018, ticks = 0.2)
grid.brackets(1220, 1020, 1057, 1020, h = 0.018, ticks = 0.1)
grid.brackets(1410, 1020, 1300, 1020, h = 0.018, ticks = 0.1)
grid.brackets(1500, 1020, 1460, 1020, h = 0.018, ticks = 0.4)
grid.brackets(1680, 1020, 1610, 1020, h = 0.018, ticks = 0.2)
grid.brackets(1960, 1020, 1890, 1020, h = 0.018, ticks = 0.2)
grid.brackets(2010, 1020, 1970, 1020, h = 0.018, ticks = 0.4)
grid.brackets(2135, 1020, 2095, 1020, h = 0.018, ticks = 0.4)
grid.brackets(2190, 1020, 2150, 1020, h = 0.018, ticks = 0.4)
dev.off()
