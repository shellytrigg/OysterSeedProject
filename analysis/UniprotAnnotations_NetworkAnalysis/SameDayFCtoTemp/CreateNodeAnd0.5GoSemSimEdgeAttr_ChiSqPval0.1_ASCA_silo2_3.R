#load libraries
library(plyr)
library(tidyr)
#install.packages("ontologyIndex")
#install.packages("ontologySimilarity")
library(ontologyIndex)
library(ontologySimilarity)
library(GSEABase)
library(reshape2)


#Cytoscape needs two files to build a network: 1. Node Attribute file and 2. Edge attribute file. 
#Node attribute file can be a list of proteins with information about each protein like alternative names, fold change and pvalue information, etc. 
#All of this data can be used to change the appearace of the nodes in the network. 
#For example, you can color nodes by their foldchange, node sizes can be based on their p-value, etc.

##########################################################
###Getting node attributes containing same day comparisons
###########################################################

#read in uniprot mapping
#uniprot <- read.csv("/Volumes/web/metacarcinus/Cgigas/all_giga-uniprot-blastP-out.nopipe.annotations.tab", sep ="\t", header = FALSE, stringsAsFactors = FALSE)
#select only some uniprot columns
#colnames(uniprot) <- c("protein_ID","Entry", "Entry_name", "perc_ident_match", "align_len", "num_mismatch", "num_gaps","querStart", "querEnd", "subjStart", "subjEnd", "evalue", "bitscore","Entry.1","Entry_name.1", "Protein_names", "Gene_names", "Organism", "Protein_length","Pathway", "GO_bp", "GO","GO_IDs", "Protein_fams")

#read in same day log FC and pval data
sameday_logFC_pval <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/TotNumSpecRatio_FC_Pval/sumNUMSPECSTOT_plus1_AllSilos_ratioFC_logFC_pval_DAYSCOMPARED.csv", stringsAsFactors = FALSE)
colnames(sameday_logFC_pval)[1] <- "protein_ID"
sameday_logFC_pval$protein_ID <- gsub("\\|","\\.",sameday_logFC_pval$protein_ID)

#combine uniprot and foldchange data
sameday_logFC_pval_uniprot<- merge(sameday_logFC_pval, uniprot, by = "protein_ID", all.x = TRUE)
#exclude proteins that didn't map to uniprot DB
sameday_logFC_pval_uniprot_mapped <- sameday_logFC_pval_uniprot[-grep("unmapped", sameday_logFC_pval_uniprot$Entry),]


#read in lists of proteins that were selected by ASCA and by 0.1 cutoff
all_sig0.1_ASCA_pro_23 <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/allS2S3_sig0.1_ASCA_pro.csv", stringsAsFactors = FALSE)
all_sig0.1_ASCA_pro_23 <- all_sig0.1_ASCA_pro_23[,c("protein_ID", "method")]
colnames(all_sig0.1_ASCA_pro_23)[2] <- "method_2v3"

#all_sig0.1_ASCA_pro_93 <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/all_sig0.1_ASCA_clust_pro.csv", stringsAsFactors = FALSE)
#all_sig0.1_ASCA_pro_93 <- all_sig0.1_ASCA_pro_93[,c("protein_ID", "method")]
#colnames(all_sig0.1_ASCA_pro_93)[2] <- "method_9v3"

#combine 2v3 and 9v3 protein selects
#all_sig0.1_ASCA_pro <- merge(all_sig0.1_ASCA_pro_23,all_sig0.1_ASCA_pro_93,by = "protein_ID", all = TRUE)
#nrow(all_sig0.1_ASCA_pro)
#231

all_sig0.1_ASCA_pro <- all_sig0.1_ASCA_pro_23
nrow(all_sig0.1_ASCA_pro)
#172

#merge with uniprot and FC info
sameday_logFC_pval_uniprot_mapped$evalue <- as.numeric(sameday_logFC_pval_uniprot_mapped$evalue)
all_sig0.1_pro_logFC_pval <- merge(sameday_logFC_pval_uniprot_mapped[which(sameday_logFC_pval_uniprot_mapped$evalue <= 10^-10),],all_sig0.1_ASCA_pro, by = "protein_ID")
nrow(all_sig0.1_pro_logFC_pval)
#146

#make column with unique shortened protien names base on gene names
all_sig0.1_pro_logFC_pval$gene_root <- gsub("_.*","", all_sig0.1_pro_logFC_pval$Entry_name)

gene_freq <- data.frame(table(gsub("_.*","", all_sig0.1_pro_logFC_pval$Entry_name)))
gene_freq <- gene_freq[which(gene_freq$Freq > 1),]

#remove choyp_gene part of protein_ID
vers <- gsub("^.*?\\.","", all_sig0.1_pro_logFC_pval$protein_ID)

#remove everything after ".m."
#vers <- gsub("\\.m\\..*$","", vers)

#add version number to gene name
all_sig0.1_pro_logFC_pval$gene_vers <- paste(all_sig0.1_pro_logFC_pval$gene, vers, sep = ".")

#only append version numbers to genes that are in the list more than once
for(i in 1:nrow(all_sig0.1_pro_logFC_pval)){
  if(all_sig0.1_pro_logFC_pval$gene_root[i] %in% gene_freq$Var1){
    all_sig0.1_pro_logFC_pval$gene[i] <- paste(all_sig0.1_pro_logFC_pval$gene_root[i], vers[i], sep = ".")
  }
  else{
    all_sig0.1_pro_logFC_pval$gene[i] <- all_sig0.1_pro_logFC_pval$gene_root[i]
  }
}

#count rows 
nrow(all_sig0.1_pro_logFC_pval)
#146

#remove extra columns (e.g. bit score, map length, GO terms)
#all_sig0.1_ASCA_pro_abbrv <- all_sig0.1_pro_logFC_pval[,-c(56:77,81,82)]
all_sig0.1_ASCA_pro_abbrv <- all_sig0.1_pro_logFC_pval[,-c(20:81)]
all_sig0.1_ASCA_pro_abbrv$type <- "protein"

#save node attribute file
write.table(all_sig0.1_ASCA_pro_abbrv,"~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1_ASCA/all_SILOS23only_sig0.1_pro_logFC_pval_abbrv_Evalcutoff_NodeAttb.txt", sep = "\t",quote = FALSE, row.names = FALSE)



############################################################

##########################################################
###Getting edge attributes; these are the GO term relationships
###########################################################
####I want to use GO slim terms instead of GO terms so the network is less busy
###So I need to map my GO terms to GO slims while keeping the protein info.
###That way, I am able to incorporate the protein data into the networks

#Get GO IDs from all_sig0.1_pro_logFC_pval
sig0.1_pro_GO <- all_sig0.1_pro_logFC_pval[,c("gene","GO_IDs")]
sig0.1_pro_GOid_term <- data.frame()
for (i in 1:nrow(sig0.1_pro_GO)){
  sig0.1_pro_GOid_term_row <- data.frame(t(data.frame(strsplit(as.character(sig0.1_pro_GO$GO_IDs[i]),'; ', fixed = TRUE))))
  sig0.1_pro_GOid_term <- rbind.fill(sig0.1_pro_GOid_term,sig0.1_pro_GOid_term_row)
}

#add protein IDs back to GO IDs
sig0.1_pro_GOid_term <- cbind(all_sig0.1_pro_logFC_pval[,"gene"], sig0.1_pro_GOid_term)
#this results in a table with protein ID listed in one column and each next column contains a GO ID that was listed with the protein in Uniprot DB.
#convert all factors to  characters
sig0.1_pro_GOid_term <- data.frame(lapply(sig0.1_pro_GOid_term, as.character), stringsAsFactors = FALSE)
#there are two proteins that don't have GO IDs

#reshape data so that all GO ID columns are gathered in one column called "GO" 
STACKED_sig0.1_pro_GOid_term <- tidyr::gather(sig0.1_pro_GOid_term,"gene","GO", 2:ncol(sig0.1_pro_GOid_term))
#exlude middle column which just contains the string "gene" in each row
STACKED_sig0.1_pro_GOid_term <- STACKED_sig0.1_pro_GOid_term[,c(1,3)]
#remove duplicate rows
STACKED_sig0.1_pro_GOid_term <- unique(STACKED_sig0.1_pro_GOid_term)
colnames(STACKED_sig0.1_pro_GOid_term)[1] <- "gene"
#remove any rows where GO column has NA value. 
STACKED_sig0.1_pro_GOid_term <- STACKED_sig0.1_pro_GOid_term[which(!is.na(STACKED_sig0.1_pro_GOid_term$GO)),]
#this resulting data frame has two columns "gene" and "GO"
write.csv(STACKED_sig0.1_pro_GOid_term, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1_ASCA/intermediate_files/STACKED_SILOS23ONLY_sig0.1_pro_GOid_term.csv", quote = FALSE, row.names = FALSE)
###Next map all GO IDs to GO slim terms


#make list of unique GO terms without a protein ID column
sig0.1_sig_GOids <- unique(STACKED_sig0.1_pro_GOid_term$GO)
#of unique proteins with GO ids
length(unique(STACKED_sig0.1_pro_GOid_term$gene))
#144

#of unique GOids
length(sig0.1_sig_GOids) #all silos: 1043
#867
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

GSEA_BP_MF <- rbind(GSEA_bp, GSEA_MF)
colnames(GSEA_BP_MF) <- c("Term", "GOslim", "GOcategory")

#next make a list of all GO ids and any other GO id that is related to it
# first make a list of all ancestor terms that correspond to GO IDs in my list

#load GO data from OntologyX package
data(go)

###this creates a list of ancestor GO IDs for each GO ID in my list
##if the GO ID is not in the "go" object from the package, the entry is "NULL"
#https://jonlefcheck.net/2013/05/20/continuing-a-for-loop-in-r-after-an-error/

term_prop <- list() # create an empty list that will get filled in by loop
for(i in 1:length(sig0.1_sig_GOids)){ # for each line in my GO IDs list
  temp_term_prop <- try(go$ancestors[[sig0.1_sig_GOids[i]]], TRUE) #make a list of all ancester GO IDs for each GO ID in my list
  if(isTRUE(class(temp_term_prop)=="try-error")) {next} else {term_prop[[i]] = temp_term_prop} # if the "go" data doesn't contain my GO ID, go on to the next GO ID. 
}

#create an empty data frame the length the list of ancestor list made above and with two columns
ancestors <- data.frame(matrix(0,length(term_prop),2))
#names the two columns
colnames(ancestors) <- c("orig.GO","GOan")

#make an empty data frame to get filled in by the loop
ances_test <- data.frame()
for(i in 1:length(term_prop)){ #for each GO ID in the ancestors list (which is the same length and order as the sig0.1_sig_GOids list)
  ancestors$orig.GO[i] <- sig0.1_sig_GOids[i] # fill in the orig. GO ID column 
  ancestors$GOan[i] <- paste(term_prop[[i]], collapse = "_") # fill in the GO ancestor column with all GO ancestor IDs corresponding to original GO term separated by an underscore
  ancestors_row <- data.frame(t(data.frame(strsplit(as.character(ancestors$GOan[i]),'_', fixed = TRUE)))) #spread ancestor IDs out across multiple columns
  ances_test <- rbind.fill(ances_test,ancestors_row) #add each row to the data frame
}
#convert factors to characters
ances_test <- data.frame(lapply(ances_test, as.character), stringsAsFactors = FALSE)
#add original GO IDs back to ancestor GO IDs
ances_test <- cbind(data.frame(sig0.1_sig_GOids, stringsAsFactors = FALSE), ances_test)

###list GO IDs not in 'go' object (don't have ancestors)
length(ances_test[which(is.na(ances_test$X1)),1]) #for all silos #9
#5

#reshape data so that all ancestor GO ID columns are gathered in one column called "Ancterm" 
STACKED_sig0.1_anc_term <- tidyr::gather(ances_test,"Ancestor","Ancterm", 2:ncol(ances_test))
STACKED_sig0.1_anc_term <- STACKED_sig0.1_anc_term[,c(1,3)]
STACKED_sig0.1_anc_term <- unique(STACKED_sig0.1_anc_term)
STACKED_sig0.1_anc_term <- STACKED_sig0.1_anc_term[which(!is.na(STACKED_sig0.1_anc_term$Ancterm)),]

# repeat the above code but this time to make a list of all child terms that correspond to GO IDs in my list
term_propC <- list()
for(i in 1:length(sig0.1_sig_GOids)){
  temp_term_propC <- try(go$children[[sig0.1_sig_GOids[i]]], TRUE)
  if(isTRUE(class(temp_term_propC)=="try-error")) {next} else {term_propC[[i]] = temp_term_propC}
}

children <- data.frame(matrix(0,length(term_propC),2))
colnames(children) <- c("orig.GO","GOch")

child_test <- data.frame()
for(i in 1:length(term_propC)){
  children$orig.GO[i] <- sig0.1_sig_GOids[i]
  children$GOch[i] <- paste(term_propC[[i]], collapse = "_")
  children_row <- data.frame(t(data.frame(strsplit(as.character(children$GOch[i]),'_', fixed = TRUE))))
  child_test <- rbind.fill(child_test,children_row)
}
#convert factors to characters
child_test <- data.frame(lapply(child_test, as.character), stringsAsFactors = FALSE)
#add original GO IDs back to child GO IDs
child_test <- cbind(data.frame(sig0.1_sig_GOids, stringsAsFactors = FALSE), child_test)

###count GO IDs not in 'go' object (these GO IDs don't have child terms)
length(child_test[which(is.na(child_test$X1)),1]) # for all silos #576
#480

#reshape data so that all child GO ID columns are gathered in one column called "Chterm" 
STACKED_sig0.1_ch_term <- tidyr::gather(child_test,"Child","Chterm", 2:ncol(child_test))
STACKED_sig0.1_ch_term <- STACKED_sig0.1_ch_term[,c(1,3)]
STACKED_sig0.1_ch_term <- unique(STACKED_sig0.1_ch_term)
STACKED_sig0.1_ch_term <- STACKED_sig0.1_ch_term[which(!is.na(STACKED_sig0.1_ch_term$Chterm)),]


# repeat the above code but this time to make a list of all parent terms that correspond to GO IDs in my list
term_propP <- list()
for(i in 1:length(sig0.1_sig_GOids)){
  temp_term_propP <- try(go$parents[[sig0.1_sig_GOids[i]]], TRUE)
  if(isTRUE(class(temp_term_propP)=="try-error")) {next} else {term_propP[[i]] = temp_term_propP}
}

parents <- data.frame(matrix(0,length(term_propP),2))
colnames(parents) <- c("orig.GO","GOpar")

par_test <- data.frame()
for(i in 1:length(term_propP)){
  parents$orig.GO[i] <- sig0.1_sig_GOids[i]
  parents$GOpar[i] <- paste(term_propP[[i]], collapse = "_")
  parents_row <- data.frame(t(data.frame(strsplit(as.character(parents$GOpar[i]),'_', fixed = TRUE))))
  par_test <- rbind.fill(par_test,parents_row)
}

#convert factors to characters
par_test <- data.frame(lapply(par_test, as.character), stringsAsFactors = FALSE)
#add original GO IDs back to ancestor GO IDs
par_test <- cbind(data.frame(sig0.1_sig_GOids, stringsAsFactors = FALSE), par_test)

###count GO IDs not in 'go' object (these GO IDs don't have parent terms)
length(par_test[which(is.na(par_test$X1)),1]) #for all silos #599
#500

#reshape data so that all child GO ID columns are gathered in one column called "Parterm" 
STACKED_sig0.1_par_term <- tidyr::gather(par_test,"Parent","Parterm", 2:ncol(par_test))
STACKED_sig0.1_par_term <- STACKED_sig0.1_par_term[,c(1,3)]
STACKED_sig0.1_par_term <- unique(STACKED_sig0.1_par_term)
STACKED_sig0.1_par_term <- STACKED_sig0.1_par_term[which(!is.na(STACKED_sig0.1_par_term$Parterm)),]

#change colnames to match before we bind these data frames together
colnames(STACKED_sig0.1_par_term)[2] <- "term"
colnames(STACKED_sig0.1_ch_term)[2] <- "term"
colnames(STACKED_sig0.1_anc_term)[2] <- "term"

#create an additional column to describe the term
STACKED_sig0.1_par_term$term_type <- "parent"
STACKED_sig0.1_ch_term$term_type <- "child"
STACKED_sig0.1_anc_term$term_type <- "ancestor"

#create another data frame with the original GO IDs in case some of these are in fact already GO slim IDs
sig0.1_orig_term <- cbind(as.data.frame(sig0.1_sig_GOids), as.data.frame(sig0.1_sig_GOids))
colnames(sig0.1_orig_term)[2] <- "term"
sig0.1_orig_term$term_type <- "orig"

#combine all ancestor, child, parent, and original GO IDs into one data frame
par_ch_anc <- rbind(STACKED_sig0.1_anc_term, STACKED_sig0.1_ch_term, STACKED_sig0.1_par_term, sig0.1_orig_term)
#this creates a data frame with three columns:
# sig0.1_sig_GOIDs, which are the original GO IDs from my list
# term, which are the Ancestor, parent, child, or original GO IDs
# type, which specifies "Ancestor", "parent", "child" or, "orig" for each GO ID

nrow(par_ch_anc) #for all silos [1] 18822
#15547

# ###################################
# ###Attempt to improve GO ID mapping to GSEA GO slim IDs by including the GO slim alt.id. This would help in the event that an ancestor, parent, child, or the GO id itself matches a GO slim alt.id rather than a GO slim ID. 
#this code is commented out because it didn't help improve the mapping; no child, parent, ancestor or original ID mapped to an alt. id.
# ##################################
# 
# #Get data
# #I downloaded goslim_generic.obo from http://www.geneontology.org/GO_slims/goslim_generic.obo
# 
# #build an index where one column is the GO slim ID and the second column is the GO slim alt. id; The GO slim ID can be listed in the first column multiple times if it has multiple GO slim alt. ids.
# #made an alt id file in command line by head -2637 ~/Downloads/goslim_generic.obo | grep "id\:" > 
# ~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1_ASCA/intermediate_files/GOslimsIDs_alts.csv; then opening in excel and adding the slim term to alt terms in another column
# ##I probably could have made a loop in R to do this, but in the interest of time I did not
# 
# 
# #read in data
# GOslims_alts <- read.csv("~/Desktop/GOslimsIDs_alts.csv", stringsAsFactors = FALSE)
# colnames(GOslims_alts)[1] <- "term"
# par_ch_anc_GOslims_alts <- merge(unique(par_ch_anc[,1:2]), GOslims_alts, by = "term")
# par_ch_anc_GOslims_uniq <- unique(par_ch_anc_GOslims_alts[,2:3])
# 
# #create a frequency table of GO slim IDs in par_ch_anc_slimBP data frame
# par_ch_anc_GOslims_uniq_freq <- data.frame(table(par_ch_anc_GOslims_uniq$GOslim)) 
# colnames(par_ch_anc_GOslims_uniq_freq)[1] <- "GOid" #rename the first column
# par_ch_anc_GOslims_uniq_freq$GOid <- as.character(par_ch_anc_GOslims_uniq_freq$GOid) #change the class of the first column from factor to character
# 
# #view the par_ch_anc_slim_freq table with the slims table to compare how many original GO ids map to GO Slim IDs
# #View(merge(slims[,c("GOid", "Count")], par_ch_anc_GOslims_uniq_freq, by = "GOid", all.x = TRUE))
# 
# #par_ch_anc_GOslims_uniq_cat <- merge(par_ch_anc_GOslims_uniq, GSEA_BP_MF, by = "GOslim")
# 
# ###I don't get more mapping from including the alternate IDs; maybe this is because they are obsolete IDs
# ##################################################


#extract only unqiue "all terms" (e.g. if an original term and an ancestor term are the same, don't list twice)
par_ch_anc_uniq <- unique(par_ch_anc[,c("sig0.1_sig_GOids","term")])
nrow(par_ch_anc_uniq)
#13338

#rename term column so merge will work (REMEMBER this term column is a list of all terms ever for each sig0.1_sig_GOids)
colnames(par_ch_anc_uniq)[2] <- "GOslim"

par_ch_anc_uniq_BP_MF <- merge(par_ch_anc_uniq, GSEA_BP_MF, by = "GOslim")
#count unique GO ids remaining
length(unique(par_ch_anc_uniq_BP_MF$sig0.1_sig_GOids)) #for all silos #814
#663

colnames(par_ch_anc_uniq_BP_MF) <- c("GOslim", "GO", "GOslimTerm", "GOcategory")

#merge par_ch_anc_slimBP and STACKED_sig0.1_pro_GOid_term and remove proteins with 'Biological Process' GO slim term since these are uninformative
goslim_protein <- merge(STACKED_sig0.1_pro_GOid_term,par_ch_anc_uniq_BP_MF, by = "GO")
#count the number of unique proteins that have GO Slim terms
length(unique(goslim_protein$gene))# started with 231 proteins, down to 197 proteins that mapped to uniprot; down to 194 that had GO IDs, down to now 189 that have GO slim
#141

#how many proteins have "biological process" term (which is a really vague, non-descriptive term)
nrow(goslim_protein[grep("GO:0008150|GO:0003674|GO:0005575", goslim_protein$GOslim),]) #for all silos #[1] 1591
#1197

#remove vague terms
goslim_protein <- goslim_protein[-grep("GO:0008150|GO:0003674|GO:0005575", goslim_protein$GOslim),]
length(unique(goslim_protein$gene))# started with 231 proteins, down to 197 proteins that mapped to uniprot; down to 194 that had GO IDs, down to now 189 that have GO slim; to 170 that have specific GO terms
#126

########
##Using GO semantic similarity to define relationships among GO slim terms. This is similar methodology to REVIGO
#####

#create a list of just unique GO Slim IDs excluding the biological process ID GO:0008150
sig0.1_GOslims <- unique(goslim_protein$GOslim)
length(sig0.1_GOslims) #FOR ALL SILOS #91
#88

#output this list so that we can upload it to REVIGO to see how these terms can be further slimmed/categorized
#write.csv(sig0.1_GOslims, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1proteins_GOsemsim_edges/sig0.1_GOslimIDs.csv", quote = FALSE, row.names = FALSE)

#make a list of GO slim IDs in a format can be used in the OntologySimilarity function
OX_sig0.1_GOslims <- list()
for(i in 1:length(sig0.1_GOslims)){
  temp_OX_sig0.1_GOslims <- try(go$id[[sig0.1_GOslims[i]]], TRUE)
  if(isTRUE(class(temp_OX_sig0.1_GOslims)=="try-error")) {next} else {OX_sig0.1_GOslims[[i]] = temp_OX_sig0.1_GOslims}
}

#creating a GO_IC formatted file for the OntologySimilarity function
GO_IC_data <- get_term_info_content(go, term_sets = OX_sig0.1_GOslims)

# create a GO semantic similarity matrix from our GO slim IDs using get_sim_grid function in the OntologySimilarity package 
sim_matrix <- get_sim_grid(
  ontology=go, 
  information_content=GO_IC_data,
  term_sets=OX_sig0.1_GOslims)
#add column and row names to the semantic similarity matrix
rownames(sim_matrix) <- sig0.1_GOslims
colnames(sim_matrix) <- sig0.1_GOslims
#convert lower triangle of matrix to NA val including the diagonal
sim_matrix[lower.tri(sim_matrix, diag = TRUE)] <- NA
#reshape the data so that each GO ID combination and semantic similarity value is listed on a different row 
term_term <- melt(sim_matrix)
nrow(term_term) #for all silos: 8281; 91 x 91 makes sense
#7744

#make a new data frame with GO ID combinations with semantic similarity values greater than 0.5. This also excludes NAs.
term_term_0.5 <- term_term[which(term_term$value > 0.5),]
nrow(term_term_0.5) #for all silos: 187
#176

#####convert GO IDs to terms####

#first create an empty data frame the length of the term_term_0.5 data frame and with 3 columns
longterm_term <- data.frame(matrix(0,nrow(term_term_0.5),3))
#name the columns
colnames(longterm_term) <- c("Var1","Var2", "value")

#Loop through each row of the term_term_0.5 data frame
for(i in 1:nrow(term_term_0.5)){
  longterm_term$Var1[i] <- GSEA_BP_MF[which(GSEA_BP_MF$GOslim == term_term_0.5$Var1[i]),"Term"] #fill in new data frame with column 1 GO ID's GO term
  longterm_term$Var2[i] <- GSEA_BP_MF[which(GSEA_BP_MF$GOslim == term_term_0.5$Var2[i]),"Term"]#fill in new data frame with column 2 GO ID's GO term
  longterm_term$value[i] <- term_term_0.5$value[i]#fill in new data frame with the semantic similarity value for the GO combination
}
#create column with interaction type information for edge attribute file
longterm_term$type <- "term-term"
#rename column
colnames(longterm_term)[3] <- "semsimValue"
nrow(longterm_term) #for all silos 187
#176

#create data frame containing protein-term data to merge with term-term data for edge attribute file
prot_term <- unique(goslim_protein[,c("gene","GOslim","GOslimTerm")])
prot_term$type <- "protein-term"
prot_term <- prot_term[,-2]
colnames(prot_term) <- c("Var1", "Var2", "type")
nrow(prot_term) #for all silos 902
#646

#merge term-term and protein-term data frames
edge_attb <- rbind.fill(longterm_term,prot_term)
nrow(edge_attb) #for all silos: 1089; makes sense because 902+187 = 1089
#822

#write out edge attribute file
#write.table(edge_attb, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1_ASCA/edge_attb_semsim0.5_ALLSILOSsig0.1ASCA_EValcutoff.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(edge_attb, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1_ASCA/edge_attb_semsim0.5_SILO23ONLYsig0.1ASCA_EValcutoff.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#####Calculate magnitude FC#############
#make list of unique GO Slims associated with these proteins
sig0.1_GOslims_terms_uniq <- unique(edge_attb[grep("protein", edge_attb$type),"Var2"])
days <- c(3,5,7,9,11,13)

go_term_magFC <- data.frame()
for (i in 1:length(sig0.1_GOslims_terms_uniq)){
  mag_df <- data.frame()
  pattern <- paste("^",sig0.1_GOslims_terms_uniq[i],"$", sep = "")
  protein_list <- prot_term[grep(pattern,prot_term$Var2),"Var1"]
  numprots <- length(protein_list)
  for(j in days){
    mag_sum_23 <- sum(abs(all_sig0.1_pro_logFC_pval[which(all_sig0.1_pro_logFC_pval$gene %in% protein_list),paste("D_",j,"_logFC_23",sep = "")]))
    #mag_sum_93 <- sum(abs(all_sig0.1_pro_logFC_pval[which(all_sig0.1_pro_logFC_pval$gene %in% protein_list),paste("D_",j,"_logFC_93",sep = "")]))
    #mag_sum_92 <- sum(abs(all_sig0.1_pro_logFC_pval[which(all_sig0.1_pro_logFC_pval$gene %in% protein_list),paste("D_",j,"_logFC_92",sep = "")]))
    mag_df[1,paste("D",j,"_magFC_23", sep = "")] <- mag_sum_23
    #mag_df[1,paste("D",j,"_magFC_93", sep = "")] <- mag_sum_93
    #mag_df[1,paste("D",j,"_magFC_92", sep = "")] <- mag_sum_92
    mag_df$numprots <- numprots
  }
  go_term_magFC <- rbind(go_term_magFC,mag_df)
}

nrow(go_term_magFC) #for all silos: 90; this is how many slim terms are associated with 223 proteins
#87

go_term_magFC$term <- sig0.1_GOslims_terms_uniq
#find min and max magnitude fold change
min(unlist(go_term_magFC[,grep("magFC", colnames(go_term_magFC))]))
#[1] 0.004883431
max(unlist(go_term_magFC[,grep("magFC", colnames(go_term_magFC))]))
#[1] 53.24725

#there is quite a difference between max and min magnitude foldchanges

#find the spread of the number of proteins per term
summary(go_term_magFC$numprots)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   2.000   5.000   7.425  10.500  42.000

#find the spread of the magnitude foldchanges
all_mags <- data.frame(unlist(go_term_magFC[,grep("magFC", colnames(go_term_magFC))]))
summary(all_mags[,1])
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.00488  0.60436  1.91399  3.90607  4.95192 53.24725 

#based on the above summary statistics, I will normalize the magnitude foldchanges by the number of proteins in each term

norm2prot <- data.frame()
for(i in 1:nrow(go_term_magFC)){
  row <- go_term_magFC[i,-grep("numprots|term", colnames(go_term_magFC))]/go_term_magFC$numprots[i]
  norm2prot <- rbind(norm2prot,row)
}
norm2prot$term <- go_term_magFC$term

go_term_magFC <- merge(go_term_magFC[,c("term", "numprots")], norm2prot, by = "term")

#write.table(go_term_magFC,"~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1_ASCA/GOnode_attb_semsim0.5_sig0.1ASCA_EValcutoff.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(go_term_magFC,"~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1_ASCA/GOnode_attb_semsim0.5_sig0.1ASCA_SILOS23ONLY_EValcutoff.txt", sep = "\t", row.names = FALSE, quote = FALSE)


