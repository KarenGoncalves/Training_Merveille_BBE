#cd $SCRATCH/Merveille/findCandidates
#
#scp cedar:/scratch/karencgs/Amaryllidoideae_proteome.fasta ./
#
#srun --account=def-desgagne --time=0:10:00 --mem=15G makeblastdb -in Amaryllidoideae_proteome.fasta  -dbtype prot -out Amaryllidoideae
#
#sh get_BBE.sh; 
#
#grep -v "#" blastp_result_270525.out | cut -f 2 | sort | uniq\
# > Amaryllidoideae_BBE_id.txt
#
#grep -f Amaryllidoideae_BBE_id.txt\
# Amaryllidoideae_proteome.fasta --no-group-separator\
# -A1 > Amaryllidoideae_BBE.fasta
#
#awk 'BEGIN {FS="\n"; OFS="\t"; RS=">"; ORS="\n"}\
# {print $1, length($2)}' Amaryllidoideae_BBE.fasta |\
# sed -E "s/^([A-Za-z0-9_\.\-]+) .+\t([0-9]+)$/\1\t\2/g" |\
# sort -k 2 > Protein_lengths_amaryllidoideae.txt

awk '$2 > 450' Protein_lengths_amaryllidoideae.txt
