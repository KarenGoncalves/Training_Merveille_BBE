library(tidyverse)
library(Biostrings)

#### Load fasta and prepare signalP prediction table
og_seqs = readAAStringSet("Sequences_for_tree.fasta", "fasta")
# For signalP, the names of the columns are in the second row
signalP = 
  read_delim('prediction_results.txt', skip=2, 
             col_names=c("ID", "Prediction", "Other", "SP_SEC_SPI", 
                         "Lipo_SEC_SPII", "TAT_Tat_SPI", "TATLIPO_Tat_SPII",
                         "PILIN_Sec_SPIII", "CS_position")) %>% 
  mutate(Cleavage_site = 
           ifelse(Prediction == "OTHER", 1, # If no SP is detected, start the protein at position 1
                  gsub("CS pos: \\d+-(\\d+)\\. .+", "\\1", 
                       CS_position) %>% as.numeric))

#  Trim the sequences, starting them according to the cleavage site
mature_sequences =
  sapply(names(og_seqs), \(peptide_id) {
    # Find start position
    Start_mature_protein = 
      (signalP %>% filter(ID == peptide_id))$Cleavage_site
    
    # Use Biostrings::subseq to indicate where the mature protein starts
    subseq(og_seqs[peptide_id], start = Start_mature_protein) %>% 
      as.character %>% unname
  }) %>% AAStringSet

mature_sequences %>% 
  writeXStringSet("mature_sequences_for_tree.fasta",
                  format = "fasta", width = 100000)

