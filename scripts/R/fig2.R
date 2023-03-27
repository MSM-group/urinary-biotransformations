#clear environment
rm(list = ls())
#load required packages
library(pacman)
pacman::p_load("tidyverse", "phyloseq", "readxl", "ggpubr", "janitor", "ggalluvial", "RColorBrewer")
#import abundance and taxonomy data
dat_in <- phyloseq::import_biom("data/biehl/mihop.biom")
#consider NA values as such (I hate using base R but tidyverse doesn't work on phyloseq objects)
tax_table(dat_in)[tax_table(dat_in) == "NA"] <- NA
#correct taxonomy table ranks
colnames(tax_table(dat_in)) <- c("domain", "phylum", "class", "order", "family", "genus", "species")
#import sample metadata
sample_dat <- read_excel("data/biehl/sample_data.xlsx") %>%
  janitor::clean_names() %>%
  sample_data()
#integrate sample metadata into phyloseq object
sample_names(sample_dat) <-  sample_dat$number_sample_id
dat <- phyloseq::merge_phyloseq(dat_in, sample_dat)
sample_names(dat) <- dat@sam_data$sid
#remove midstream urine samples
dat <- phyloseq::subset_samples(dat, sample_type != "MU")
#normalize data using rarefaction
dat_norm <- phyloseq::rarefy_even_depth(dat, rngseed = 420)
#Agglomerate at the genus level
dat_norm_genera <- phyloseq::tax_glom(dat_norm, taxrank = "genus", NArm = FALSE)
dat_norm_genera_sites <- phyloseq::merge_samples(dat_norm_genera, group = "sample_type", fun = mean)
#Caulculate relative abundance of genera
dat_norm_genera_sites_rel <- phyloseq::transform_sample_counts(dat_norm_genera_sites, function(x) x/sum(x))
#Bin low abundance phyla as other
dat_norm_genera_sites_rel_low_abundance <- phyloseq::filter_taxa(dat_norm_genera_sites_rel, function(x) sum(x) < 0.01, TRUE)
merge_taxa <- phyloseq::taxa_names(dat_norm_genera_sites_rel_low_abundance)
dat_norm_genera_sites_rel <- phyloseq::merge_taxa(dat_norm_genera_sites_rel, merge_taxa, archetype = 1)
#Merge Eubacteria into one genus
eubacterium <- phyloseq::tax_table(dat_norm_genera_sites_rel) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "asv") %>%
  dplyr::filter(stringr::str_detect(genus, "Eubacterium")) %>%
  pull(asv)
dat_norm_genera_sites_rel <- phyloseq::merge_taxa(dat_norm_genera_sites_rel, eubacterium, archetype = 1)
#Convert phyloseq format to tibble for plot
#Create function to handle merged taxa
unknown_vs_other <- function(a, b){
  dplyr::case_when(is.na(b) ~ dplyr::case_when(a == merge_taxa[1] ~ "unknown",
                                               a == eubacterium[1] ~ "Eubacterium",
                                               TRUE ~ "other"),
                   TRUE ~ b)
}
#convert phyloseq taxonomic data to tibble
tax_tab <- phyloseq::tax_table(dat_norm_genera_sites_rel) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "asv") %>%
  dplyr::select(asv, genus) %>%
  dplyr::mutate(genus = unknown_vs_other(asv, genus))
#convert phyloseq abundance data to tibble
asv_tab <- phyloseq::otu_table(dat_norm_genera_sites_rel) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "asv")
#merge taxonomic and abundance data
tab <- left_join(tax_tab, asv_tab, by = "asv") %>%
  select(-asv) %>%
  dplyr::rename(urine = CU, vagina = vswab, gut = stool) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(total = sum(urine, gut, vagina)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(genus = forcats::fct_reorder(genus, total, .desc = TRUE)) %>%
  dplyr::mutate(genus = forcats::fct_relevel(as.factor(genus), "other", after = Inf)) %>%
  dplyr::mutate(genus = forcats::fct_relevel(as.factor(genus), "unknown", after = Inf)) %>%
  dplyr::select(-total)
#convert data to appropriate format for alluvial plot
alluvial_data <- tab %>%
  pivot_longer(-genus, names_to = "site", values_to = "count")
#make alluvial plot
n_colors <- length(levels(alluvial_data$genus))-2
get_palette <- colorRampPalette(brewer.pal(8, "Set1"))
set.seed(121)
pal <- c(sample(get_palette(n_colors), n_colors, replace = FALSE), "gray70", "gray30")
names(pal) <- levels(alluvial_data$genus)
alluvial <- ggplot2::ggplot(alluvial_data, ggplot2::aes(y = count, axis1 = site, axis2 = genus, fill = genus)) +
  ggalluvial::geom_alluvium(width = 1/3) +
  ggalluvial::geom_stratum(width = 1/3, fill = "gray80", color = "gray50", size = 0.1) +
  ggplot2::scale_x_discrete(limits = c("genus", "site")) +
  ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, min.y = 0.05, fontface = "italic") +
  ggplot2::scale_fill_manual(values = pal) +
  ggpubr::theme_pubr() +
  ggplot2::theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.position = "none")
#save plot
ggplot2::ggsave("plots/fig2.jpg",
                alluvial,
                device = "jpeg",
                dpi = 300,
                units = "cm",
                width = 18.8,
                height = 9.4)
