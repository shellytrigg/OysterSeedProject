#load libraries
library(plyr)
library(tidyr)
#install.packages("ontologyIndex")
#install.packages("ontologySimilarity")
library(ontologyIndex)
library(ontologySimilarity)
library(GSEABase)
library(reshape2)


#####WITH SR LAB GO SLIMS#####
#load Robert's lab go slim terms
srlabGOSLIM <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/background/GO-GOslim.sorted.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(srlabGOSLIM)<- c("GOid", "GOterm", "GOslim", "category")
#load protein/GOid data
STACKED_sig0.1_pro_GOid_term <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1_ASCA/intermediate_files/STACKED_sig0.1_pro_GOid_term.csv", stringsAsFactors = FALSE)
colnames(STACKED_sig0.1_pro_GOid_term)[2] <- "GOid"
STACKED_sig0.1_pro_GOid_term_srlab <- merge(STACKED_sig0.1_pro_GOid_term,srlabGOSLIM[,c(1,3)], by = "GOid",all.x = TRUE)
##############################

uniqsrlabGOSLIMs <- unique(STACKED_sig0.1_pro_GOid_term_srlab$GOslim)
length(uniqsrlabGOSLIMs)
#37
length(unique(srlabGOSLIM$GOslim))
#39; only 38 unique GO slim terms exist in the Roberts lab GO slim list and these include CC, MP and BP; they aren't the same terms as in the .obo file

#semantic sim with SR lab GO slim terms

beach2 <- list()
for(i in 1:length(uniqsrlabGOSLIMs)){
  temp_beach2 <- try(go$name[[uniqsrlabGOSLIMs[i]]], TRUE)
  if(isTRUE(class(temp_beach2)=="try-error")) {next} else {beach2[[i]] = temp_beach2}
}

test <- get_term_info_content(go, term_sets = uniqsrlabGOSLIMs)
####
##CAN'T MAP SR LAB GO SLIM TERMS TO TERMS IN goslim_generic.OBO BECAUSE THE TERMS DON'T MATCH!!!!!
###WOULD NEED THE GO SLIM IDs THAT MATCH THE SR LAB GO SLIM TERMS IN ORDER TO GET THE SEMANTIC SIMILARITIES TO MAKE THE GO SLIM TERM-TERM RELATIONSHIPS

