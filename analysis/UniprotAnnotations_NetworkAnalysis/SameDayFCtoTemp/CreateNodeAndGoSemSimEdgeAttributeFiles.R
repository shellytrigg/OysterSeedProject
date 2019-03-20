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
uniprot <- read.csv("/Volumes/web/metacarcinus/Cgigas/all_giga-uniprot-blastP-out.nopipe.annotations.tab", sep ="\t", header = FALSE, stringsAsFactors = FALSE)
#select only some uniprot columns
colnames(uniprot) <- c("protein_ID","Entry", "Entry_name", "perc_ident_match", "align_len", "num_mismatch", "num_gaps","querStart", "querEnd", "subjStart", "subjEnd", "evalue", "bitscore","Entry.1","Entry_name.1", "Protein_names", "Gene_names", "Organism", "Protein_length","Pathway", "GO_bp", "GO","GO_IDs", "Protein_fams")

#read in same day log FC and pval data
sameday_logFC_pval <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/TotNumSpecRatio_FC_Pval/sumNUMSPECSTOT_plus1_ratioFC_logFC_pval_DAYSCOMPARED.csv", stringsAsFactors = FALSE)
colnames(sameday_logFC_pval)[1] <- "protein_ID"
#combine uniprot and foldchange data
sameday_logFC_pval_uniprot<- merge(sameday_logFC_pval, uniprot, by = "protein_ID", all.x = TRUE)
#exclude proteins that didn't map to uniprot DB
sameday_logFC_pval_uniprot_mapped <- sameday_logFC_pval_uniprot[-grep("unmapped", sameday_logFC_pval_uniprot$Entry),]

#####select only proteins with adj Chi sq. pvalue <= 0.1####

#create a list of all column names with adj.Chisq.pval 
adjChiSqpvalColumns <- colnames(sameday_logFC_pval_uniprot_mapped)[grep("adj.ChiSq.pval",colnames(sameday_logFC_pval_uniprot_mapped))]

#build a list of protiens with adj Chi sq. pvalue <= 0.1
#create empty data frame the loop will add too
all_sig_pro <- data.frame()
for (i in 1:length(adjChiSqpvalColumns)){ # for each name in adj.Chisq.pval column name list
  column <- adjChiSqpvalColumns[i] # create a variable for indexed column name
  #make a data frame containing protein IDs for all proteins in indexed column that have adj.Chisq.pval <=0.1
  sig_pro <- data.frame(sameday_logFC_pval_uniprot_mapped[which(sameday_logFC_pval_uniprot_mapped[,column] <= 0.1),1])
  #iteratively add protein lists to initial data frame
  all_sig_pro <- rbind(all_sig_pro, sig_pro)
}

#count how many unique proteins are in the list of proteins with adj Chi sq. pvalue <= 0.1
nrow(unique(all_sig_pro))
#[1] 153

#make a data frame of just unique proteins so we can select these from the foldchange/pvalue data
all_sig0.1_pro <- unique(all_sig_pro)
#select sig proteins @ p.adj 0.1 from logFC list
all_sig0.1_pro_logFC_pval <- sameday_logFC_pval_uniprot_mapped[which(sameday_logFC_pval_uniprot_mapped$protein_ID %in% all_sig0.1_pro[,1]),]

##Export protein list as node attributes table to upload in cytoscape
#remove extra columns (e.g. bit score, map length, GO terms)
all_sig0.1_pro_logFC_pval_abbrv <- all_sig0.1_pro_logFC_pval[,-c(22:29,31:ncol(all_sig0.1_pro_logFC_pval))]
#add column to describe node type for cytoscape
all_sig0.1_pro_logFC_pval_abbrv$type <- "protein"
#save node attribute file
write.csv(all_sig0.1_pro_logFC_pval_abbrv,"~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1proteins_GOsemsim_edges/all_sig0.1_pro_logFC_pval_abbrv_NodeAttb.csv", quote = FALSE, row.names = FALSE)

############################################################

##########################################################
###Getting edge attributes; these are the GO term relationships
###########################################################

#Get GO IDs by subsetting all_sig0.1_pro_logFC_pval for "protein_ID" column and column with list of "GO_IDs" separated by ";"
sig0.1_pro_GO <- all_sig0.1_pro_logFC_pval[,c("protein_ID","GO_IDs")]
#make empty data frame that will get filled in by for loop
sig0.1_pro_GOid_term <- data.frame()
#loop through each line of data frame containing one column with "protein_ID" and one column with list of "GO_IDs" separated by ";"
for (i in 1:nrow(sig0.1_pro_GO)){
  #create a row that lists all GO IDs associated with one protein listed across rows, where each different GO ID gets put in it's own column
  sig0.1_pro_GOid_term_row <- data.frame(t(data.frame(strsplit(as.character(sig0.1_pro_GO$GO_IDs[i]),'; ', fixed = TRUE))), stringsAsFactors = FALSE)
  #add each row created in the line above to the empty data frame
  #this will add NAs in columns where less terms exist; e.g. if a protein only has two terms and another has 10, there will be 8 NAs added to the row for the protein with 2 terms
  sig0.1_pro_GOid_term <- rbind.fill(sig0.1_pro_GOid_term,sig0.1_pro_GOid_term_row)
}

#add protein IDs back to GO IDs
sig0.1_pro_GOid_term <- cbind(all_sig0.1_pro_logFC_pval[,"protein_ID"], sig0.1_pro_GOid_term)
#this results in a table with protein ID listed in one column and each next column contains a GO ID that was listed with the protein in Uniprot DB.
str(sig0.1_pro_GOid_term)

#reshape data so that all GO ID columns are gathered in one column called "GO" 
STACKED_sig0.1_pro_GOid_term <- tidyr::gather(sig0.1_pro_GOid_term,"protein_ID","GO", 2:ncol(sig0.1_pro_GOid_term))
#exlude middle column which just contains the string "protein_ID" in each row
STACKED_sig0.1_pro_GOid_term <- STACKED_sig0.1_pro_GOid_term[,c(1,3)]
#remove duplicate rows
STACKED_sig0.1_pro_GOid_term <- unique(STACKED_sig0.1_pro_GOid_term)
colnames(STACKED_sig0.1_pro_GOid_term)[1] <- "protein_ID"
#remove any rows where GO column has NA value. 
STACKED_sig0.1_pro_GOid_term <- STACKED_sig0.1_pro_GOid_term[which(!is.na(STACKED_sig0.1_pro_GOid_term$GO)),]
#this resulting data frame has two columns "protein_ID" and "GO"

###Next map all GO IDs to GO slim terms
#make list of unique GO terms without a protein ID column
sig0.1_sig_GOids <- unique(STACKED_sig0.1_pro_GOid_term$GO)
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
#this results in a data frame with columns: 
# "GOid" which is the GO slim GO ID
# "Term" which is the GO slim term corresponding to the GO ID
# "Percent" which is the percent of GO IDs in my list that mapped to each GO slim term
# "Count" which is the number of GO IDs in my list that mapped to each GO slim term
# I couldn't figure out how to extract the original GO ID and slim GO ID mapping from this function,
# So I went on to search all ancestor, parent, and child terms for each original GO ID in my list
# for GO slim IDs


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

#add original GO IDs back to ancestor GO IDs
ances_test <- cbind(sig0.1_sig_GOids, ances_test)

###list GO IDs not in 'go' object (don't have ancestors)
ances_test[which(is.na(ances_test$X1)),1]
#[1] GO:0062023 GO:0103025 GO:0102102 GO:0102131 GO:0090736 GO:1905905 GO:0061844 GO:1905907 GO:0106036

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

child_test <- cbind(sig0.1_sig_GOids, child_test)

###count GO IDs not in 'go' object (these GO IDs don't have child terms)
length(child_test[which(is.na(child_test$X1)),1])
#366

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

par_test <- cbind(sig0.1_sig_GOids, par_test)

###count GO IDs not in 'go' object (these GO IDs don't have parent terms)
length(par_test[which(is.na(par_test$X1)),1])
#405

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

nrow(par_ch_anc)
#[1] 12429


#pull out all GO IDs from 'par_ch_anc' list that are GO slim IDs
par_ch_anc_slimBP <- par_ch_anc[which(par_ch_anc$term %in% slims$GOid),]

#how many proteins have "biological process" term (which is a really vague, non-descriptive term)
nrow(par_ch_anc_slimBP[grep("GO:0008150", par_ch_anc_slimBP$term),])
#[1] 379

#checking how GSEA and ontologyX GOid frequncies compare 
par_ch_anc_slim_freq <- data.frame(table(par_ch_anc_slimBP$term)) #create a frequency table of GO slim IDs in par_ch_anc_slimBP data frame
colnames(par_ch_anc_slim_freq)[1] <- "GOid" #rename the first column
par_ch_anc_slim_freq$GOid <- as.character(par_ch_anc_slim_freq$GOid) #change the class of the first column from factor to character

#view the par_ch_anc_slim_freq table with the slims table to compare how many original GO ids map to GO Slim IDs
View(merge(slims[,c("GOid", "Count")], par_ch_anc_slim_freq, by = "GOid", all = TRUE))
###There is a difference with how GSEAbase calculates counts and the GOid frequencies from ontologyX, but the frequencies are relatively similar
#this provides a little explanation as to how GSEA is calculating frquencies: https://support.bioconductor.org/p/100403/
#for now, I'm going to go with ontologyX counts

#make a list of protein IDs and their slim terms
colnames(par_ch_anc_slimBP) <- c("GO", "GOslim") #first rename the columns in par_ch_anc_slimBP so it can be merged
#merge par_ch_anc_slimBP and STACKED_sig0.1_pro_GOid_term and remove proteins with 'Biological Process' GO slim term since these are uninformative
goslim_protein <- merge(STACKED_sig0.1_pro_GOid_term,par_ch_anc_slimBP[-grep("GO:0008150",par_ch_anc_slimBP$GOslim),], by = "GO")
#count the number of unique proteins that have GO Slim terms
length(unique(goslim_protein$protein_ID))
#[1] 102
###There are 51 proteins that are excluded. They seem to be excluded because they don't have GO BP terms, only cell component or molec. function terms

#add slim terms to table in addition to their IDs
colnames(goslim_protein)[3] <- "GOid"
goslim_protein <- merge(goslim_protein,slims[,3:4], by = "GOid")
colnames(goslim_protein)[1] <- "GOslim"
colnames(goslim_protein)[5] <- "GOslimTerm"
#this results in a table with 5 columns:
# GOslim ID
# Original GO ID 
# protein_ID
# term type (orig, ancestor, child, parent)
# GO slim term

########
##Using GO semantic similarity to define relationships among GO slim terms. This is similar methodology to REVIGO
#####

#create a list of just unique GO Slim IDs excluding the biological process ID GO:0008150
sig0.1_GOslims <- unique(par_ch_anc_slim_freq[-grep("GO:0008150",par_ch_anc_slim_freq$GOid),"GOid"])
#output this list so that we can upload it to REVIGO to see how these terms can be further slimmed/categorized
write.csv(sig0.1_GOslims, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1proteins_GOsemsim_edges/sig0.1_GOslimIDs.csv", quote = FALSE, row.names = FALSE)

#make a list of GO slim IDs in a format can be used in the OntologySimilarity function
beach <- list()
for(i in 1:length(sig0.1_GOslims)){
  temp_beach <- try(go$id[[sig0.1_GOslims[i]]], TRUE)
  if(isTRUE(class(temp_beach)=="try-error")) {next} else {beach[[i]] = temp_beach}
}

#creating a GO_IC formatted file for the OntologySimilarity function
test <- get_term_info_content(go, term_sets = beach)

# create a GO semantic similarity matrix from our GO slim IDs using get_sim_grid function in the OntologySimilarity package 
sim_matrix <- get_sim_grid(
  ontology=go, 
  information_content=test,
  term_sets=beach)
#add column and row names to the semantic similarity matrix
rownames(sim_matrix) <- sig0.1_GOslims
colnames(sim_matrix) <- sig0.1_GOslims
#convert lower triangle of matrix to NA val including the diagonal
sim_matrix[lower.tri(sim_matrix, diag = TRUE)] <- NA
#reshape the data so that each GO ID combination and semantic similarity value is listed on a different row 
term_term <- melt(sim_matrix)
#make a new data frame with GO ID combinations with semantic similarity values greater than 0.7. This also excludes NAs.
term_term_0.7 <- term_term[which(term_term$value > 0.7),]

#####convert GO IDs to terms####

#first create an empty data frame the length of the term_term_0.7 data frame and with 3 columns
longterm_term <- data.frame(matrix(0,nrow(term_term_0.7),3))
#name the columns
colnames(longterm_term) <- c("Var1","Var2", "value")

#Loop through each row of the term_term_0.7 data frame
for(i in 1:nrow(term_term_0.7)){
  longterm_term$Var1[i] <- slims[which(slims$GOid == term_term_0.7$Var1[i]),"Term"] #fill in new data frame with column 1 GO ID's GO term
  longterm_term$Var2[i] <- slims[which(slims$GOid == term_term_0.7$Var2[i]),"Term"]#fill in new data frame with column 2 GO ID's GO term
  longterm_term$value[i] <- term_term_0.7$value[i]#fill in new data frame with the semantic similarity value for the GO combination
}
#create column with interaction type information for edge attribute file
longterm_term$type <- "term-term"
#rename column
colnames(longterm_term)[3] <- "semsimValue"

#create data frame containing protein-term data to merge with term-term data for edge attribute file
prot_term <- unique(goslim_protein[,c("protein_ID","GOslim","GOslimTerm")])
prot_term$type <- "protein-term"
prot_term <- prot_term[,-2]
colnames(prot_term) <- c("Var1", "Var2", "type")
#merge term-term and protein-term data frames
edge_attb <- rbind.fill(longterm_term,prot_term)
#write out edge attribute file
write.csv(edge_attb, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1proteins_GOsemsim_edges/edge_attb_semsim0.7_sig0.1.csv", row.names = FALSE, quote = FALSE)
