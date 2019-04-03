library(tidyr)
library(GSEABase)
library(plyr)

#read in data
NSAF <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/silo3and9_nozerovals_AVGs.csv", stringsAsFactors = FALSE)

#I had already changed the zeros to 0.1, so I will sub that back with zero for this analysis
NSAF[NSAF == 0.1] <- 0

# #order NSAF data by temp and remove day 0
# NSAF <- NSAF[which(NSAF$day != 0),]
# NSAF <- NSAF[order(NSAF$temp),]
# 
# NSAF_t <- data.frame(t(NSAF[,-c(1,2)]))
# colnames(NSAF_t) <- paste(NSAF$day, NSAF$temp, sep = "_")
# 

#separate 3 and 9

silo3 <- NSAF[grep("23", NSAF$temp),]
silo9 <- NSAF[grep("29", NSAF$temp),]

#reformat to long list
silo3_t <- data.frame(t(silo3[,-c(1,2)]))
colnames(silo3_t) <- silo3$day

#remove proteins that were only in silo9
View(silo3_t[which(apply(silo3_t, 1, var) == 0),])
no_val_proteins <- rownames(silo3_t[which(apply(silo3_t, 1, var) == 0),])

View(silo3_t[-c(which(rownames(silo3_t) %in% no_val_proteins)),])
silo3_t_nozero <- silo3_t[-c(which(rownames(silo3_t) %in% no_val_proteins)),]
#confirm it worked
nrow(silo3_t)-nrow(silo3_t_nozero)
length(no_val_proteins)

#add protein column
silo3_t_nozero$protein_ID <- rownames(silo3_t_nozero)
silo3_t_nozero <- silo3_t_nozero[grep("CHOYP", silo3_t_nozero$protein_ID),]

#merge with uniprot
#uniprot <- read.csv("/Volumes/web/metacarcinus/Cgigas/all_giga-uniprot-blastP-out.nopipe.annotations.tab", sep ="\t", header = FALSE, stringsAsFactors = FALSE)
#combine uniprot and NSAF data
silo3_t_nozero_unip<- merge(silo3_t_nozero, uniprot, by = "protein_ID", all.x = TRUE)
#exclude proteins that didn't map to uniprot DB
silo3_t_nozero_unip_mapped <- silo3_t_nozero_unip[-grep("unmapped", silo3_t_nozero_unip$Entry),]
nrow(silo3_t_nozero_unip_mapped)
silo3_t_nozero_unip_mapped <- silo3_t_nozero_unip_mapped[which(silo3_t_nozero_unip_mapped$evalue <= 10^-10),]
nrow(silo3_t_nozero_unip_mapped)

#go terms
#Subset "protein_ID" column and "GO IDs" column containing list of GO_IDs separated by ";" from all_sig0.1_pro_logFC_pval 
pro_GO <- silo3_t_nozero_unip_mapped[,c("protein_ID","GO_IDs")]

pro_GOid_term <- data.frame()
#loop through each line of data frame containing one column with "protein_ID" and one column with list of "GO_IDs" separated by ";"
for (i in 1:nrow(pro_GO)){
  #create a row that lists all GO IDs associated with one protein listed across rows, where each different GO ID gets put in it's own column
  pro_GOid_term_row <- data.frame(t(data.frame(strsplit(as.character(pro_GO$GO_IDs[i]),'; ', fixed = TRUE))), stringsAsFactors = FALSE)
  #add each row created in the line above to the empty data frame
  #this will add NAs in columns where less terms exist; e.g. if a protein only has two terms and another has 10, there will be 8 NAs added to the row for the protein with 2 terms
  pro_GOid_term <- rbind.fill(pro_GOid_term,pro_GOid_term_row)
}

#add protein IDs back to GO IDs
pro_GOid_term <- cbind(silo3_t_nozero_unip_mapped[,"protein_ID"], pro_GOid_term)
#this results in a table with protein ID listed in one column and each next column contains a GO ID that was listed with the protein in Uniprot DB.
str(pro_GOid_term)

#reshape data so that all GO ID columns are gathered in one column called "GO" 
STACKED_pro_GOid_term <- tidyr::gather(pro_GOid_term,"protein_ID","GO", 2:ncol(pro_GOid_term))
#exlude middle column which just contains the string "protein_ID" in each row
STACKED_pro_GOid_term <- STACKED_pro_GOid_term[,c(1,3)]
#remove duplicate rows
STACKED_pro_GOid_term <- unique(STACKED_pro_GOid_term)
colnames(STACKED_pro_GOid_term)[1] <- "protein_ID"
#remove any rows where GO column has NA value. 
STACKED_pro_GOid_term <- STACKED_pro_GOid_term[which(!is.na(STACKED_pro_GOid_term$GO)),]
#this resulting data frame has two columns "protein_ID" and "GO"

###Next map all GO IDs to GO slim terms
#make list of unique GO terms without a protein ID column
GOids <- unique(STACKED_pro_GOid_term$GO)
#goslims with GSEA
myCollection <- GOCollection(GOids)
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

write.csv(slims, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/silo3_goslims.csv", row.names = FALSE, quote = FALSE)
