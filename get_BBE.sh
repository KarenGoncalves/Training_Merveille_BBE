#!/bin/sh

DOSSIER=$SCRATCH/Merveille/findCandidates/
mySequence=$DOSSIER/baits_BBE.fasta
database=$DOSSIER/Amaryllidoideae

########################
##### Run the code #####
########################

module load StdEnv/2023 blast+

cd $DOSSIER

######################################################
##### If you have a protein sequence, use blastp #####
######################################################

srun --account=def-desgagne --time=0:05:00 --mem=1G\
 blastp -query ${mySequence}\
 -db ${database}\
 -max_target_seqs 10\
 -outfmt '7'\
 > $DOSSIER/blastp_result_270525.out
 
