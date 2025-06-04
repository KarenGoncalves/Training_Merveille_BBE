## Telecharger et ouvrir les librairies ##
devtools::source_gist("https://gist.github.com/KarenGoncalves/0db105bceff4ff69547ee25460dda978")

install_from_dif_sources(
  cran_packages = c("tidyverse", "tinytex", "patchwork"),
  bioconductor_packages = c("Biostrings", "msa", "treeio", "ggtree", "ape", "seqinr", "phangorn"),
  github_packages = "YuLab-SMU/ggmsa"
)
fasta_name = "Sequences_for_tree.fasta"
## Ouvrir le fichier fasta ##
# To open a fasta file with multiple sequences, use the function readDNAStringSet() or readAAStringSet()
fasta_for_alignment = fasta_name %>%
  readAAStringSet(format = "fasta")

## Multiple sequence alignment ##
# use ?msa to know what you need to put in the function and what the default values are
myFirstAlignment <- msa(fasta_for_alignment, 
                        method = "ClustalOmega",
                        verbose = T, maxiters = 10
)
### clustalo --R -o tempClustalOmega.aln --outfmt=clustal --seqtype=protein --force --gapopen=6.000000 --gapext=1.000000 --cluster-size=100 --iter=0 --output-order=tree-order

## Tree with bootstrap ##
## We need to transform the alignment into another type of R file
myFirstAlignment2 <- 
    myFirstAlignment %>% as.phyDat()


# Then, we create a distance matrix
bionj_tree <- 
  dist.ml(x = myFirstAlignment2, model = "JTT") %>% 
  bionj()

fit_ml <- 
  pml(bionj_tree, data = myFirstAlignment2, model = "JTT")

fitJC <- optim.pml(fit_ml, model = "JTT", rearrangement="stochastic")

anova(fit_ml, fitJC)
AIC(fit_ml)
AIC(fitJC)

bootML_tree <- 
  bootstrap.pml(fitJC, bs=300, 
                optNni=TRUE,
                control = pml.control(trace = 0))
cnet <- consensusNet(bootML_tree, p = 0.5)
plot(cnet, show.edge.label = T)

write.tree(bootML_tree, 
           paste0("Maximum_likelihood_wholeProts_", Sys.Date(), ".nwk"))

# bs_fastme_tree <- bootstrap.phyDat(
#   myFirstAlignment2, 
#   bs = 300, \(x){
#     dist.ml(x = x, 
#             model = "Blosum62") |> 
#       fastme.bal()
#   })
# 
# bs_plot <- plotBS(fastme_tree, bs_fastme_tree)
