#cd $SCRATCH/Merveille/Training_Merveille_BBE
#
#scp cedar:/scratch/karencgs/Amaryllidoideae_proteome.fasta ./
# module load StdEnv/2023 blast+
#srun --account=def-desgagne --time=0:10:00 --mem=15G makeblastdb -in Amaryllidoideae_proteome.fasta  -dbtype prot -out Amaryllidoideae
#
#sh scripts/get_BBE.sh; 

grep -v "#" blastp_result_280525.out | cut -f 2 | sort | uniq\
 > Amaryllidoideae_BBE_id.txt

grep -f Amaryllidoideae_BBE_id.txt\
 Amaryllidoideae_proteome.fasta --no-group-separator\
 -A1 > Amaryllidoideae_BBE.fasta

awk 'BEGIN {FS="\n"; OFS="\t"; RS=">"; ORS="\n"}\
 {if ($1 ~ /internal/) {type="internal"}\
 else if ($1 ~ /3prime_partial/) {type="3prime_partial"}\
 else if ($1 ~ /5prime_partial/) {type="5prime_partial"}\
 else if ($1 ~ /complete/) {type="complete"};\
 print $1, length($2), type}'\
 Amaryllidoideae_BBE.fasta |\
 sed -E "s/^([A-Za-z0-9_\.\-]+) .+\t([0-9]+)\t([a-z0-9_]+)$/\1\t\2\t\3/g" |\
 sort -k 2 > Protein_lengths_amaryllidoideae.txt

awk '$2 > 450 {print $1}' Protein_lengths_amaryllidoideae.txt > Amaryllidoideae_BBE450_id.txt

grep -f Amaryllidoideae_BBE450_id.txt blastp_result_280525.out\
 | awk '$3 > 35 && $4 > 400 { print $2}' |\
 sort | uniq |\
 sed -E 's/^(.+)(_i[0-9]+|_seq[0-9]+)*(\.p[0-9]+)/\1\t\1\2\3/' > Gene_to_protein.txt

sed 's/*//' Amaryllidoideae_BBE.fasta > Amaryllidoideae_BBE_noStop.fasta
