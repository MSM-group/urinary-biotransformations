#!/bin/bash
accession=( $(tail -n +2 data/file_lists/PRJNA316969_accession.csv | cut -d ',' -f1) )
cd data/thomas_white_assemblies_2
for i in "${accession[@]}"; do
	url="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$i/download?include_annotation_type=GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA&filename=$i.zip"
	curl -OJX GET ${url} -H "Accept: application/zip"
done