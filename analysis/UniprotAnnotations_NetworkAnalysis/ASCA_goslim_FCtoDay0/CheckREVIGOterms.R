#load libraries
library(dplyr)

ncols <- max(count.fields("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/ASCA_goslim_FCtoDay0/ASCA_GOslim.sif", sep = "\t"))
MEDsif <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/ASCA_goslim_FCtoDay0/ASCA_GOslim.sif", header = FALSE , sep = "\t",col.names = paste0("V", seq_len(ncols)), stringsAsFactors = FALSE)

#subset out GO slim terms from REVIGO .xgmml file
MED_REVIGO_term_list_unique <- unique(as.character(unlist(MEDsif[1:77,c(1,3)])))
#subset out GO slim terms from OntologyX GO slim analysis
MED_uniprot_term_list_unique <- unique(unlist(MEDsif[78:700, 3]))

terms_only_in_REVIGO <- MED_REVIGO_term_list_unique[which(!(MED_REVIGO_term_list_unique %in% MED_uniprot_term_list_unique))]
terms_only_in_uniprot_list <- MED_uniprot_term_list_unique[which(!(MED_uniprot_term_list_unique %in% MED_REVIGO_term_list_unique))]

length(terms_only_in_REVIGO)
#[1] 17

length(terms_only_in_uniprot_list)
#[1] 26
