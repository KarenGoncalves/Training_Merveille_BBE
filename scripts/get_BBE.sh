#!/bin/sh

DOSSIER=$SCRATCH/Merveille/Training_Merveille_BBE/
mySequence=$DOSSIER/baits_BBE.fasta
database=$DOSSIER/Amaryllidoideae

########################
##### Run the code #####
########################

module load StdEnv/2023 blast+/2.14.1

cd $DOSSIER

######################################################
##### If you have a protein sequence, use blastp #####
######################################################


srun --account=def-desgagne --time=0:05:00 --mem=1G\
 blastp -query ${mySequence}\
 -db ${database}\
 -evalue 0.0005\
 -outfmt '6'\
 > $DOSSIER/blastp_result_280525.out
