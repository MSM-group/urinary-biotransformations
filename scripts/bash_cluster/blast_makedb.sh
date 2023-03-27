#!/bin/bash
#BSUB -J "thomas_white_aa_blast_db"
#BSUB -R "rusage[mem=64]"
#BSUB -n 1
#BSUB -W 4:00
#BSUB -o log_files/thomas_white_aa_blast_db.out
#BSUB -e error_files/thomas_white_aa_blast_db_error.out

#read csv file with file names
accession=( $(tail -n +2 file_lists/PRJNA316969_accession.csv | cut -d ',' -f1) )

module load blast-plus/2.9.0
for i in "${accession[@]}"; do
	makeblastdb -in data/thomas_white_assemblies/ncbi_dataset/data/$i/protein.faa -dbtype prot -out data/reverse_blast_database/$i
done