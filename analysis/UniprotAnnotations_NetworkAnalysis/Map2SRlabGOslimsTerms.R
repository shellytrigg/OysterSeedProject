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
srlabGOSLIM <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/background/GOSlim_terms.csv", stringsAsFactors = FALSE)
colnames(srlabGOSLIM)[1]<- "GO"
#load protein/GOid data
STACKED_sig0.1_pro_GOid_term <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1_ASCA/intermediate_files/STACKED_sig0.1_pro_GOid_term.csv", stringsAsFactors = FALSE)
STACKED_sig0.1_pro_GOid_term_srlab <- merge(STACKED_sig0.1_pro_GOid_term,srlabGOSLIM[,c(1,3)], by = "GO",all.x = TRUE)
##############################

uniqsrlabGOSLIMs <- unique(STACKED_sig0.1_pro_GOid_term_srlab$GOSlim_bin)
length(uniqsrlabGOSLIMs)
#37
length(unique(srlabGOSLIM$GOSlim_bin))
#38; only 38 unique GO slim terms exist in the Roberts lab GO slim list and these include CC, MP and BP
#I can't map srlab GOslim terms to GOIDs; they aren't the same terms as in the .obo file

#need to get GO IDs for Roberts Lab slim terms
#Use GSEA to generate list of all GO Slim BP, MP, and CC
#DF will have "term", "GOid", "GOcategory"

#BP first
#goslims with GSEA
myCollection <- GOCollection(sig0.1_sig_GOids)
#I downloaded goslim_generic.obo from http://geneontology.org/docs/go-subset-guide/
#then i moved it to the R library for GSEABase in the extdata folder
fl <- system.file("extdata", "goslim_generic.obo", package="GSEABase")
slim <- getOBOCollection(fl)
slims <- data.frame(goSlim(myCollection, slim, "BP"))
slims$GOid <- rownames(slims)
slims$Term <- as.character(slims$Term)
rownames(slims) <- NULL

GSEA_bp <- slims[,c("Term", "GOid")]
GSEA_bp$GOcategory <- "BP"

#MP next
slims <- data.frame(goSlim(myCollection, slim, "MF"))
slims$GOid <- rownames(slims)
slims$Term <- as.character(slims$Term)
rownames(slims) <- NULL

GSEA_MF <- slims[,c("Term", "GOid")]
GSEA_MF$GOcategory <- "MF"

#CC next
slims <- data.frame(goSlim(myCollection, slim, "CC"))
slims$GOid <- rownames(slims)
slims$Term <- as.character(slims$Term)
rownames(slims) <- NULL

GSEA_CC <- slims[,c("Term", "GOid")]
GSEA_CC$GOcategory <- "CC"

GSEA_BP_MF_CC <- rbind(GSEA_bp, GSEA_MF, GSEA_CC)
colnames(GSEA_BP_MF_CC) <- c("Term", "GOslim", "GOcategory")

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

