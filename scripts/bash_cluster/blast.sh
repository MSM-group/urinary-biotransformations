#!/bin/bash
#BSUB -J "thomas_white_blast_2"
#BSUB -R "rusage[mem=1024]"
#BSUB -n 1
#BSUB -W 4:00
#BSUB -o log_files/thomas_white_blast_2.out
#BSUB -e error_files/thomas_white_blast_error_2.out

accession=( $(tail -n +2 file_lists/PRJNA316969_accession.csv | cut -d ',' -f1) )

module load blast-plus/2.9.0
for i in "${accession[@]}"; do
	blastp -query data/21_ec_reps_namfix_including_drugbug.fasta -db data/reverse_blast_database/$i -out output/thomas_white_blast_3/${i}_blast_3.out -outfmt "6 qseqid qlen saccver ssciname pident qcovs length evalue bitscore" -evalue 0.1 -qcov_hsp_perc 20
done