#### Telecharger et ouvrir les librairies ####
devtools::source_gist("https://gist.github.com/KarenGoncalves/0db105bceff4ff69547ee25460dda978")

install_from_dif_sources(
  cran_packages = c("tidyverse", "tinytex", "patchwork"),
  bioconductor_packages = c("Biostrings", "msa", "treeio", "ggtree", "ape", "seqinr", "phangorn"),
  github_packages = "YuLab-SMU/ggmsa"
)

##### Criteria ######
minimum_length = 450
min_alnLength = 400
min_pct_id = 35

##### Open files #####
protein_info = read_delim("Protein_lengths_amaryllidoideae.txt", skip=1, 
                             col_names = c("ProtID", "Length", "Type")) %>% 
  mutate(GeneID = gsub("(_i[0-9]+|_seq[0-9]+)*\\.p[0-9]+",
                       "", ProtID))

Amaryllidaceae_hits = 
  "Amaryllidoideae_BBE_noStop.fasta" %>%
  readAAStringSet(format = "fasta")
names(Amaryllidaceae_hits) = gsub(" .+", "", names(Amaryllidaceae_hits))

baits =  
  "baits_BBE.fasta" %>%
  readAAStringSet(format = "fasta")

blastp_result = 
  read_delim("blastp_result_270525.out",
             col_names = c("Bait", "Hit", "Pct_identity", "Aln_length",
                           "Mismatches", "gap", "Bait_start", "Bait_end",
                           "Hit_start", "Hit_end", "evalue", "bitscore"))
###### Filter data #####
protein_minlength =
  protein_info %>% 
  filter(Length > minimum_length)

blastp_best_alignments =
  blastp_result %>% 
  filter(Pct_identity > min_pct_id,
         Aln_length > min_alnLength)

# Find protein with best matches and characteristics in blastp result
grouped_blastp_result =
  left_join(protein_minlength, blastp_best_alignments,
            by = join_by("ProtID" == "Hit")) %>% 
  group_by(Bait, GeneID) %>% 
  arrange(Pct_identity, Aln_length, Length, ProtID) %>% 
  slice_head(n = 1)

selected_prots = grouped_blastp_result$ProtID %>% unique


###### Filter fasta
Amaryllidaceae_candidates = 
  Amaryllidaceae_hits[selected_prots]

# Clean sequence names and export data
Amaryllidaceae_candidates_cleanNames = 
  names(Amaryllidaceae_candidates) %>% 
  gsub("_TRINITY", "", x = .) %>% 
  gsub("-[A-Za-z].+(\\.p[0-9]+)", "\\1", x = .)

names(Amaryllidaceae_candidates) = Amaryllidaceae_candidates_cleanNames

clean_baits = names(baits) %>% 
  gsub("^sp\\|[0-9A-Z\\.]+\\|", "", x = .) %>% 
  gsub(" ", "", x = .)

names(baits) = clean_baits

c(Amaryllidaceae_candidates, baits) %>%
  writeXStringSet(filepath = "Sequences_for_tree.fasta", 
                  width=100000)
