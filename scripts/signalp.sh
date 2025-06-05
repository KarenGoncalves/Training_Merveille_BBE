'#!/bin/sh
#SBATCH --account=def-desgagne
#SBATCH --job-name=signalP_BBE
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --output=/scratch/karencgs/Merveille/Training_Merveille_BBE/signalP-%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G

module restore annotation_modules
THREADS=8
        


cd /scratch/karencgs/Merveille/Training_Merveille_BBE
ml StdEnv/2023 cuda scipy-stack/2023b
source ~/ENV2023/bin/activate

signalp6 --fastafile Sequences_for_tree.fasta\
 --model_dir ~/ENV2023/lib/python3.11/site-packages/signalp/model_weights\
 --organism other\
 --output_dir ./\
 --format none\
 --write_procs ${THREADS}\
 --mode fast

grep -v "#" prediction_results.txt |\
 awk -F '\t' '{if ($2 == "SP") {pos=$9} else {pos="1\t0"};\
 print $1, pos, pr}' |\
 sed -E 's/CS pos: [0-9]+-([0-9]+). Pr: ([0-9\.]+)/\1\t\2/'\
 > signalP_prediction.txt
