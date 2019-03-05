#load libraries
library(plyr)
library(tidyr)
#install.packages("ontologyIndex")
library(ontologyIndex)

#merge Uniprot mapping data with foldchange data
FC_logFC_pval <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/TotNumSpecRatio_FC_Pval/sumNUMSPECSTOT_plus1_ratioFC_logFC_pval.csv", stringsAsFactors = FALSE)
colnames(FC_logFC_pval)[1] <- "protein_ID"
###read in uniprot table with GO information
uniprot <- read.csv("/Volumes/web/metacarcinus/Cgigas/all_giga-uniprot-blastP-out.nopipe.annotations.tab", sep ="\t", header = FALSE, stringsAsFactors = FALSE)
#name columns with meaningful names
colnames(uniprot) <- c("protein_ID","Entry", "Entry_name", "perc_ident_match", "align_len", "num_mismatch", "num_gaps","querStart", "querEnd", "subjStart", "subjEnd", "evalue", "bitscore","Entry.1","Entry_name.1", "Protein_names", "Gene_names", "Organism", "Protein_length","Pathway", "GO_bp", "GO","GO_IDs", "Protein_fams")
#combine uniprot and fold change/pvalue data
FC_logFC_pval_uniprot<- merge(FC_logFC_pval, uniprot, by = "protein_ID", all.x = TRUE)
#separate out unmapped proteins 
FC_logFC_pval_uniprot_unmapped <- FC_logFC_pval_uniprot[which(FC_logFC_pval_uniprot$Entry == "unmapped"),]
nrow(FC_logFC_pval_uniprot_unmapped)
#[1] 270
#separate out mapped proteins
FC_logFC_pval_uniprot_mapped <- FC_logFC_pval_uniprot[-grep("unmapped", FC_logFC_pval_uniprot$Entry),]

#first attempt at playing with data, used ASCA temperature affected proteins only from Feb 5.
#ASCA_temp <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_all_proteins_avgADJNSAF/ASCA_TempAffectedProteins.csv", stringsAsFactors = FALSE)
#ASCA_temp_t <- data.frame(t.data.frame(ASCA_temp[,-1]))
#colnames(ASCA_temp_t) <- ASCA_temp$X
#ASCA_temp_t$protein_ID <- rownames(ASCA_temp_t)
#ASCA_temp_annos <- merge(ASCA_temp_t, uniprot, by = "protein_ID", all.x = TRUE)

###Analysis with all ASCA proteins (time, temp, and time x temp affected)
ASCA_all <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_all_proteins_avgADJNSAF/ASCA_TimeTempTimexTemp_proteins_loadings.csv")
colnames(ASCA_all)[1] <- "protein_ID"
#make a data frame with ASCA selected proteins and their fold change, pvalue, and uniprot mapping data
ASCA_all_FC_pval <- merge(ASCA_all, FC_logFC_pval_uniprot_mapped, by = "protein_ID", all.x = TRUE)

#create a data frame of ASCA proteins and their GO IDs
ASCA_pro_GO <- ASCA_all_FC_pval[,c("protein_ID","GO_IDs")]
ASCA_pro_GOid_term <- data.frame()
for (i in 1:nrow(ASCA_pro_GO)){
  ASCA_pro_GOid_term_row <- data.frame(t(data.frame(strsplit(as.character(ASCA_pro_GO$GO_IDs[i]),'; ', fixed = TRUE))))
  ASCA_pro_GOid_term <- rbind.fill(ASCA_pro_GOid_term,ASCA_pro_GOid_term_row)
}

#add protein IDs back to GO terms
ASCA_pro_GOid_term <- cbind(ASCA_all_FC_pval[,"protein_ID"], ASCA_pro_GOid_term)

#reshape data
STACKED_ASCA_pro_GOid_term <- tidyr::gather(ASCA_pro_GOid_term,"protein_ID","GO", 2:ncol(ASCA_pro_GOid_term))
STACKED_ASCA_pro_GOid_term <- STACKED_ASCA_pro_GOid_term[,c(1,3)]
STACKED_ASCA_pro_GOid_term <- unique(STACKED_ASCA_pro_GOid_term)
colnames(STACKED_ASCA_pro_GOid_term)[1] <- "protein_ID"
STACKED_ASCA_pro_GOid_term <- STACKED_ASCA_pro_GOid_term[which(!is.na(STACKED_ASCA_pro_GOid_term$GO)),]

write.csv(data.frame(unique(STACKED_ASCA_pro_GOid_term$GO)), "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations/ASCA_goslim_FCtoDay0/ASCA_all_GOterms.csv", row.names = FALSE, quote = FALSE)

###
#get ancestors
ASCA_sig_GOids <- unique(STACKED_ASCA_pro_GOid_term$GO)
#goslims wiht GSEA
myCollection <- GOCollection(ASCA_sig_GOids)
#I downloaded goslim_generic.obo from http://geneontology.org/docs/go-subset-guide/
#then i moved it to the R library for GSEABase in the extdata folder
fl <- system.file("extdata", "goslim_generic.obo", package="GSEABase")
slim <- getOBOCollection(fl)
slims <- data.frame(goSlim(myCollection, slim, "BP"))
slims$GOid <- rownames(slims)
rownames(slims) <- NULL

###this creates a list of ancestor GO terms for each GO term in ASCA list
##if the GO term is not in the "go" object, the entry is "NULL"
##for example, GO:0103046 is not in the "go" object so in the term_prop list
#> term_prop[[457]]
#NULL
#> ASCA_sig_GOids[[457]]
#[1] "GO:0103046"
#https://jonlefcheck.net/2013/05/20/continuing-a-for-loop-in-r-after-an-error/

term_prop <- list()
for(i in 1:length(ASCA_sig_GOids)){
  temp_term_prop <- try(go$ancestors[[ASCA_sig_GOids[i]]], TRUE)
  if(isTRUE(class(temp_term_prop)=="try-error")) {next} else {term_prop[[i]] = temp_term_prop}
}

ancestors <- data.frame(matrix(0,1190,2))
colnames(ancestors) <- c("orig.GO","GOan")

ances_test <- data.frame()
for(i in 1:length(term_prop)){
  ancestors$orig.GO[i] <- ASCA_sig_GOids[i]
  ancestors$GOan[i] <- paste(term_prop[[i]], collapse = "_")
  ancestors_row <- data.frame(t(data.frame(strsplit(as.character(ancestors$GOan[i]),'_', fixed = TRUE))))
  ances_test <- rbind.fill(ances_test,ancestors_row)
}

ances_test <- cbind(ASCA_sig_GOids, ances_test)

###GO terms not in 'go' object (don't have ancestors)
ances_test[which(is.na(ances_test$X1)),1]
#[1] GO:0103046 GO:0103025 GO:0062023 GO:0101031 GO:0061844 GO:0106037

STACKED_ASCA_anc_term <- tidyr::gather(ances_test,"Ancestor","Ancterm", 2:ncol(ances_test))
STACKED_ASCA_anc_term <- STACKED_ASCA_anc_term[,c(1,3)]
STACKED_ASCA_anc_term <- unique(STACKED_ASCA_anc_term)
STACKED_ASCA_anc_term <- STACKED_ASCA_anc_term[which(!is.na(STACKED_ASCA_anc_term$Ancterm)),]

STACKED_ASCA_anc_term_slim <- STACKED_ASCA_anc_term[which(STACKED_ASCA_anc_term$Ancterm %in% slims$GOid),]

#compare GO slim frequencies with ancestors
#make frequency table from ancestors
anc_slim <- data.frame(table(STACKED_ASCA_anc_term_slim$Ancterm))
colnames(anc_slim)[1] <- "GOid"

View(merge(slims[,c("GOid", "Count")], anc_slim, all = TRUE))


###try adding children terms
term_propC <- list()
for(i in 1:length(ASCA_sig_GOids)){
  temp_term_propC <- try(go$children[[ASCA_sig_GOids[i]]], TRUE)
  if(isTRUE(class(temp_term_propC)=="try-error")) {next} else {term_propC[[i]] = temp_term_propC}
}

children <- data.frame(matrix(0,1190,2))
colnames(children) <- c("orig.GO","GOch")

child_test <- data.frame()
for(i in 1:length(term_propC)){
  children$orig.GO[i] <- ASCA_sig_GOids[i]
  children$GOch[i] <- paste(term_propC[[i]], collapse = "_")
  children_row <- data.frame(t(data.frame(strsplit(as.character(children$GOch[i]),'_', fixed = TRUE))))
  child_test <- rbind.fill(child_test,children_row)
}

child_test <- cbind(ASCA_sig_GOids, child_test)

###GO terms not in 'go' object (don't have ancestors)
length(child_test[which(is.na(child_test$X1)),1])
#627

STACKED_ASCA_ch_term <- tidyr::gather(child_test,"Child","Chterm", 2:ncol(child_test))
STACKED_ASCA_ch_term <- STACKED_ASCA_ch_term[,c(1,3)]
STACKED_ASCA_ch_term <- unique(STACKED_ASCA_ch_term)
STACKED_ASCA_ch_term <- STACKED_ASCA_ch_term[which(!is.na(STACKED_ASCA_ch_term$Chterm)),]


###try adding parent terms
term_propP <- list()
for(i in 1:length(ASCA_sig_GOids)){
  temp_term_propP <- try(go$parents[[ASCA_sig_GOids[i]]], TRUE)
  if(isTRUE(class(temp_term_propP)=="try-error")) {next} else {term_propP[[i]] = temp_term_propP}
}

parents <- data.frame(matrix(0,1190,2))
colnames(parents) <- c("orig.GO","GOpar")

par_test <- data.frame()
for(i in 1:length(term_propP)){
  parents$orig.GO[i] <- ASCA_sig_GOids[i]
  parents$GOpar[i] <- paste(term_propP[[i]], collapse = "_")
  parents_row <- data.frame(t(data.frame(strsplit(as.character(parents$GOpar[i]),'_', fixed = TRUE))))
  par_test <- rbind.fill(par_test,parents_row)
}

par_test <- cbind(ASCA_sig_GOids, par_test)

###GO terms not in 'go' object (don't have ancestors)
length(par_test[which(is.na(par_test$X1)),1])
#681

STACKED_ASCA_par_term <- tidyr::gather(par_test,"Parent","Parterm", 2:ncol(par_test))
STACKED_ASCA_par_term <- STACKED_ASCA_par_term[,c(1,3)]
STACKED_ASCA_par_term <- unique(STACKED_ASCA_par_term)
STACKED_ASCA_par_term <- STACKED_ASCA_par_term[which(!is.na(STACKED_ASCA_par_term$Parterm)),]

#change colnames to match
colnames(STACKED_ASCA_par_term)[2] <- "term"
colnames(STACKED_ASCA_ch_term)[2] <- "term"
colnames(STACKED_ASCA_anc_term)[2] <- "term"

#create an additional column to describe the term
STACKED_ASCA_par_term$term_type <- "parent"
STACKED_ASCA_ch_term$term_type <- "child"
STACKED_ASCA_anc_term$term_type <- "ancestor"
ASCA_orig_term <- cbind(as.data.frame(ASCA_sig_GOids), as.data.frame(ASCA_sig_GOids))
colnames(ASCA_orig_term)[2] <- "term"
ASCA_orig_term$term_type <- "orig"


par_ch_anc <- rbind(STACKED_ASCA_anc_term, STACKED_ASCA_ch_term, STACKED_ASCA_par_term, ASCA_orig_term)
nrow(par_ch_anc)
#[1] 21771

par_ch_anc_slimBP <- par_ch_anc[which(par_ch_anc$term %in% slims$GOid),]
#how many proteins have "biological process" term
nrow(par_ch_anc_slimBP[grep("GO:0008150", par_ch_anc_slimBP$term),])
#[1] 690

#unique terms to upload to Revigo
write.csv(unique(par_ch_anc_slim[-grep("GO:0008150",par_ch_anc_slim$term),2]),"~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations/ASCA_goslim_FCtoDay0/ASCA_all_GOSLIMterms.csv")

#checking how GSEA and ontologyX GOid frequncies compare 
par_ch_anc_slim_freq <- data.frame(table(par_ch_anc_slimBP$term))
colnames(par_ch_anc_slim_freq)[1] <- "GOid"

View(merge(slims[,c("GOid", "Count")], par_ch_anc_slim_freq, by = "GOid", all = TRUE))
###There is a difference with how GSEAbase calculates counts and the GOid frequencies from ontologyX
#this provides a little explanation as to how GSEA is calculating it: https://support.bioconductor.org/p/100403/
#for now, I'm going to go with ontologyX counts

#make a list of protein-slim terms
colnames(par_ch_anc_slimBP) <- c("GO", "GOslim")
goslim_protein <- merge(STACKED_ASCA_pro_GOid_term,par_ch_anc_slimBP[-grep("GO:0008150",par_ch_anc_slimBP$GOslim),], by = "GO")
length(unique(goslim_protein$protein_ID))
#[1] 186
###There are 221 proteins that are excluded. They seem to be excluded because they don't have GO BP terms, only cell component or molec. function terms

#add slim terms to table
colnames(goslim_protein)[3] <- "GOid"
goslim_protein <- merge(goslim_protein,slims[,3:4], by = "GOid")
colnames(goslim_protein)[1] <- "GOslim"
colnames(goslim_protein)[4] <- "GOslimTerm"

###merge GO slim terms with uniprot info
#goslim_protein_uniprot <- merge(goslim_protein,ASCA_all_FC_pval[,-grep("GO", colnames(ASCA_all_FC_pval))], by = "protein_ID")
#write.csv(goslim_protein_uniprot,"~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations/ASCA_goslim_FCtoDay0/ASCA_all_proteins_GOSLIMterms_uniprot.csv", row.names = FALSE, quote = FALSE)

goslim_protein_uniprot <- ASCA_all_FC_pval[which(ASCA_all_FC_pval$protein_ID %in% unique(goslim_protein$protein_ID)),-grep("GO|_FC|_Chi", colnames(ASCA_all_FC_pval))]

#export table for nodes attributes

write.csv(goslim_protein_uniprot[,-c(29,31:38,40:48)], "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations/ASCA_goslim_FCtoDay0/ASCA_all_UNIQUE_proteins_wGSLIM_uniprot.csv", row.names = FALSE, quote = FALSE)


################################
################################
#Everything below is from Feb 5 preliminary analysis
################################
################################


###find all unique proteins which have adj Chisq pval < 0.05
adjChiSqpvalColumns <- colnames(FC_logFC_pval_uniprot_mapped)[grep("adj.ChiSq.pval",colnames(FC_logFC_pval_uniprot_mapped))]

all_sig_pro <- data.frame()
for (i in 1:length(adjChiSqpvalColumns)){
  column <- adjChiSqpvalColumns[i]
  sig_pro <- data.frame(FC_logFC_pval_uniprot_mapped[which(FC_logFC_pval_uniprot_mapped[,column] <0.05),1])
  all_sig_pro <- rbind(all_sig_pro, sig_pro)
}

nrow(all_sig_pro)

all_sig_pro <- unique(all_sig_pro)
nrow(all_sig_pro)

sig_FC_logFC_pval_uniprot_mapped <- FC_logFC_pval_uniprot_mapped[which(FC_logFC_pval_uniprot_mapped$protein_ID %in% all_sig_pro[,1]),]

write.table(sig_FC_logFC_pval_uniprot_mapped, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations/sigProteins_FC_logFC_pval_uniprotAnnos.tsv", sep = "\t",quote = FALSE, row.names = FALSE)


write.table(sig_FC_logFC_pval_uniprot_mapped[,grep("adj|Entry$")], "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations/sigProteins_FC_logFC_pval_uniprotAnnos.tsv", sep = "\t",quote = FALSE, row.names = FALSE)
write.csv(sig_FC_logFC_pval_uniprot_mapped[,c(50,grep("adj", colnames(sig_FC_logFC_pval_uniprot_mapped)))], "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations/sigProteins_FC_logFC_pval_uniprotAccessions.csv", quote = FALSE, row.names = FALSE)

######
##make a heatmap of log fold change values for all sig proteins 
####

logFC_sig_uniprotmapped <- sig_FC_logFC_pval_uniprot_mapped[,grep("logFC", colnames(sig_FC_logFC_pval_uniprot_mapped))]
rownames(logFC_sig_uniprotmapped) <- sig_FC_logFC_pval_uniprot_mapped[,1]
logFC_sig_uniprotmapped_m <- as.matrix(logFC_sig_uniprotmapped)
heatmap3(logFC_sig_uniprotmapped_m, cexCol = 0.7)

D3 <- sig_FC_logFC_pval_uniprot_mapped[which(sig_FC_logFC_pval_uniprot_mapped$D_3_T_23_adj.ChiSq.pval < 0.05 | sig_FC_logFC_pval_uniprot_mapped$D_3_T_29_adj.ChiSq.pval < 0.05),grep("logFC", colnames(sig_FC_logFC_pval_uniprot_mapped))]
rownames(D3) <- sig_FC_logFC_pval_uniprot_mapped[which(sig_FC_logFC_pval_uniprot_mapped$D_3_T_23_adj.ChiSq.pval < 0.05 | sig_FC_logFC_pval_uniprot_mapped$D_3_T_29_adj.ChiSq.pval < 0.05),1]
D3_m <- as.matrix(D3[,grep("D_3", colnames(D3))])
heatmap3(D3_m, cexRow = 0.5,cexCol = 0.7, scale = "none")

###ANNOTATIONS ARE NOT FILTERED BY BLASTp e-values !!!!!!!!!
###make list of proteins and each GO term on separate line

pro_GO <- sig_FC_logFC_pval_uniprot_mapped[,c("protein_ID","GO")]
pro_GOid_term <- data.frame()
for (i in 1:nrow(pro_GO)){
  pro_GOid_term_row <- data.frame(t(data.frame(strsplit(as.character(pro_GO$GO[i]),'; ', fixed = TRUE))))
  pro_GOid_term <- rbind.fill(pro_GOid_term,pro_GOid_term_row)
}

pro_GOid_term <- cbind(sig_FC_logFC_pval_uniprot_mapped[,"protein_ID"], pro_GOid_term)
STACKED_pro_GOid_term <- tidyr::gather(pro_GOid_term,"protein_ID","GO", 2:ncol(pro_GOid_term))
STACKED_pro_GOid_term <- STACKED_pro_GOid_term[,c(1,3)]
STACKED_pro_GOid_term <- unique(STACKED_pro_GOid_term)
colnames(STACKED_pro_GOid_term)[1] <- "protein_ID"
STACKED_pro_GOid_term <- STACKED_pro_GOid_term[which(!is.na(STACKED_pro_GOid_term$GO)),]


STACKED_pro_GOid_term <- STACKED_pro_GOid_term %>% separate(GO, c("GO_term", "GO_ID"), " \\[GO:")
STACKED_pro_GOid_term[,"GO_ID"] <- gsub("]","", STACKED_pro_GOid_term[,"GO_ID"])
STACKED_pro_GOid_term[,"GO_ID"] <- paste("GO:",STACKED_pro_GOid_term$GO_ID, sep = "")

sig_uniq_GO_IDs <- data.frame(paste("GO:",unique(STACKED_pro_GOid_term$GO_ID), sep = ""))


write.csv(sig_uniq_GO_IDs, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations/sigProt_unique_GOids.csv",row.names = FALSE, quote = FALSE)

#####my own enrichment

D3T23_sig <- sig_FC_logFC_pval_uniprot_mapped[which(sig_FC_logFC_pval_uniprot_mapped$D_3_T_23_adj.ChiSq.pval < 0.05),1:2]
D3T23_sig_GO <- STACKED_pro_GOid_term[which(STACKED_pro_GOid_term$protein_ID %in% D3T23_sig$protein_ID),]

D3T23_sig_GO_freq <- data.frame(table(D3T23_sig_GO$GO_ID))
#convert factors to characters
D3T23_sig_GO_freq$Var1 <- as.character(D3T23_sig_GO_freq$Var1)

All_sig_GO_freq <- data.frame(table(STACKED_pro_GOid_term$GO_ID))
#convert factors to characters
All_sig_GO_freq$Var1 <- as.character(All_sig_GO_freq$Var1)


background_terms <- nrow(STACKED_pro_GOid_term)
D3T23_terms <- nrow(D3T23_sig_GO)

GO_enrich <- data.frame()
for (i in 1:nrow(D3T23_sig_GO_freq)){
  sample_GO <- D3T23_sig_GO_freq[i,1]
  sample_GO_freq <- D3T23_sig_GO_freq[i,2]
  all_GO_freq <- All_sig_GO_freq[which(All_sig_GO_freq$Var1 == sample_GO),2]
  x <- data.frame(sample = c(sample_GO_freq, D3T23_terms), all = c(all_GO_freq, background_terms))
  y <- fisher.test(x)
  p.val <- y$p.value
  GO_enrich_row <- cbind.data.frame(sample_GO, sample_GO_freq, all_GO_freq, p.val)
  GO_enrich <- rbind.data.frame(GO_enrich, GO_enrich_row)
}
GO_enrich$adj.p.val <- p.adjust(GO_enrich$p.val)
GO_enrich$sample_GO <- as.character(GO_enrich$sample_GO)

write.csv(GO_enrich, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations/GOenrich_D3T23.csv",row.names = FALSE, quote = FALSE)


#####FOR 29 DEGREES

D3T29_sig <- sig_FC_logFC_pval_uniprot_mapped[which(sig_FC_logFC_pval_uniprot_mapped$D_3_T_29_adj.ChiSq.pval < 0.05),1:2]
D3T29_sig_GO <- STACKED_pro_GOid_term[which(STACKED_pro_GOid_term$protein_ID %in% D3T29_sig$protein_ID),]

D3T29_sig_GO_freq <- data.frame(table(D3T29_sig_GO$GO_ID))
#convert factors to characters
D3T29_sig_GO_freq$Var1 <- as.character(D3T29_sig_GO_freq$Var1)

D3T29_terms <- nrow(D3T29_sig_GO)

GO_enrich <- data.frame()
for (i in 1:nrow(D3T29_sig_GO_freq)){
  sample_GO <- D3T29_sig_GO_freq[i,1]
  sample_GO_freq <- D3T29_sig_GO_freq[i,2]
  all_GO_freq <- All_sig_GO_freq[which(All_sig_GO_freq$Var1 == sample_GO),2]
  x <- data.frame(sample = c(sample_GO_freq, D3T29_terms), all = c(all_GO_freq, background_terms))
  y <- fisher.test(x)
  p.val <- y$p.value
  GO_enrich_row <- cbind.data.frame(sample_GO, sample_GO_freq, all_GO_freq, p.val)
  GO_enrich <- rbind.data.frame(GO_enrich, GO_enrich_row)
}
GO_enrich$adj.p.val <- p.adjust(GO_enrich$p.val)
GO_enrich$sample_GO <- as.character(GO_enrich$sample_GO)

write.csv(GO_enrich, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations/GOenrich_D3T29.csv",row.names = FALSE, quote = FALSE)





######
##Make nodes and edge table to add to cytoscape network
#####

D3T23_sig_prot_nodes <- sig_FC_logFC_pval_uniprot_mapped[which(sig_FC_logFC_pval_uniprot_mapped$D_3_T_23_adj.ChiSq.pval < 0.05),c(1:2,5)]
colnames(D3T23_sig_prot_nodes) <- c("protein_ID", "log10p.adj", "logFC")
D3T23_sig_prot_nodes$log10p.adj <- log(D3T23_sig_prot_nodes$log10p.adj, 10)
write.csv(D3T23_sig_prot_nodes, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations/D3T23_sig_prot_nodes.csv",row.names = FALSE, quote = FALSE)


D3T29_sig_prot_nodes <- sig_FC_logFC_pval_uniprot_mapped[which(sig_FC_logFC_pval_uniprot_mapped$D_3_T_29_adj.ChiSq.pval < 0.05),c(1,6,9)]
colnames(D3T29_sig_prot_nodes) <- c("protein_ID", "log10p.adj", "logFC")
D3T29_sig_prot_nodes$log10p.adj <- log(D3T29_sig_prot_nodes$log10p.adj, 10)

#combine D3 protein node info
D3_sig_prot_nodes <- rbind(D3T29_sig_prot_nodes, D3T29_sig_prot_nodes)

#make edge attribute table

write.csv(D3T23_sig_GO[,1:2], "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations/D3T23_sig_prot_edges.csv",row.names = FALSE, quote = FALSE)




#######
##Alternatively, convert GO terms to GO slim and compare across samples without enrichment
#########



GO_assoc <- data.frame()


https://www.slideshare.net/dgrapov/proteomics-workshop-2014-lab-dmitry-grapov


ftp://ftp.geneontology.org/go/www/GO.downloads.files.shtml

#heatmap based on pfams

library(heatmap3)

ASCA_temp_annos_pfams <- ASCA_temp_annos[,c(2:14,37)]


test <- heatmap3(as.matrix(cut_data_ord_t), Colv = NA, cexRow = 0.5)

dendro_indx <- data.frame(test$rowInd)
dendro_indx$dir <- c(rep("down",64),rep("up",74))
dendro_indx$dendro_list_num <- as.numeric(rownames(dendro_indx))
colnames(dendro_indx)[1] <- "levels_num"
dendro_indx <- dendro_indx[order(dendro_indx$levels_num),]
dendro_indx$protein_ID <- rownames(cut_data_ord_t)
colnames(dendro_indx)[4] <- "protein_ID"

View(merge(dendro_indx,ASCA_temp_annos[,c("protein_ID", "Protein_fams")], by = "protein_ID"))

cut_data_ord_t <- data.frame(t.data.frame(cut_data_ord))
#make row names which contain the protein IDs a column
cut_data_ord_t <- cbind(rownames(cut_data_ord_t), cut_data_ord_t)
#add the column name
colnames(cut_data_ord_t)[1] <- "Protein.ID"
#convert the Protein ID column to character
cut_data_ord_t[,1] <- as.character(cut_data_ord_t[,1])
cut_data_ord_t <- cut_data_ord_t[,-1]
colnames(cut_data_ord_t) <- rownames(cut_data_ord)


