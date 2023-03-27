#clear environment
rm(list = ls())
#load required packages
library(pacman)
pacman::p_load("tidyverse", "readxl", "janitor", "pdftools")
#Retrieve strain names from supplemental table in paper
metadat <- readxl::read_excel("data/metadata/Copy_of_41467_2018_3968_MOESM3_ESM.xlsx") %>%
  janitor::clean_names()  %>%
  dplyr::mutate(genera = word(species, sep = " ", start = 1, end = 1))
#Load NCBI bioproject info table for PRJNA316969
ncbi_acc <- readr::read_tsv("data/metadata/PRJNA316969_AssemblyDetails.txt", skip = 1) %>%
  janitor::clean_names() %>%
  dplyr::mutate(strain = stringr::str_to_upper(strain))
#Select only strains from asymptomatic patients
metadat_controls <- metadat %>%
  dplyr::filter(symptom_status == "Asymptomatic")
controls <- ncbi_acc %>%
  filter(strain %in% metadat_controls$strain_collection_id)
#Save list of accessions
out_list <- controls %>%
  dplyr::select(number_assembly) %>%
  dplyr::distinct()
readr::write_csv(out_list, "data/metadata/PRJNA316969_accession.csv")

