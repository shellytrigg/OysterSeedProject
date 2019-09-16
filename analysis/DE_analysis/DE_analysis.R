library(heatmap3)
library(ggplot2)
library(colorRamps)
library(plotrix)
library(randomcoloR)
library(plyr)

#read in avg NSAF same day logFC data
avgNSAF_logFC <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/TotNumSpecRatio_FC_Pval/avgADJNSAF_logFC_DAYSCOMPARED.csv", stringsAsFactors = FALSE)
avgNSAF_logFC <- avgNSAF_logFC[grep("CHOYP",avgNSAF_logFC$protein_ID),]


#find proteins with 0 FC:
zeroFC <- avgNSAF_logFC[which(apply(avgNSAF_logFC[,-7],1,sd) == 0),]

sd_df <- data.frame(apply(avgNSAF_logFC[,-7],1,sd))
colnames(sd_df)[1] <- "sd"
sd_df$protein <- avgNSAF_logFC$protein_ID
ggplot(sd_df) + geom_density(aes(sd))
ggplot(sd_df) + geom_histogram(aes(sd), binwidth = 0.1)

nearZeroFC <- avgNSAF_logFC[which(apply(avgNSAF_logFC[,-7],1,sd) <= 0.2),]
nearZeroFC_STACKED <- tidyr::gather(nearZeroFC, "day", "logFC_NSAF", 1:6)
nearZeroFC_STACKED$day <- gsub("D_","", nearZeroFC_STACKED$day)
nearZeroFC_STACKED$day <- gsub("_logFC_NSAF","", nearZeroFC_STACKED$day)
nearZeroFC_STACKED$day <- as.numeric(as.character(nearZeroFC_STACKED$day))


ggplot(data.frame(nearZeroFC_STACKED),aes(day,logFC_NSAF)) + geom_line(aes(group = protein_ID)) + scale_x_continuous(breaks = c(3,5,7,9,11,13), labels = c(3,5,7,9,11,13)) + theme_bw()


l#remove protein ID and save matrix
NSAF_matrix <- as.matrix(avgNSAF_logFC[which(!(avgNSAF_logFC$protein_ID %in% zeroFC$protein_ID)),-7])
rownames(NSAF_matrix) <- avgNSAF_logFC[which(!(avgNSAF_logFC$protein_ID %in% zeroFC$protein_ID)),7]


heatmap3(NSAF_matrix,cexRow = 0.1, cexCol = 0.5, Colv = NA, method = "average")
hm <- as.dist(1-cor(t(NSAF_matrix), use="pa"))
hm2 <- hclust(hm, method = 'average')
plot(hm2)
mycl <- cutree(hm2, h=0.4)
clusterCols <- distinctColorPalette(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
heatmap3(NSAF_matrix,cexRow = 0.5, cexCol = 0.5, RowSideColors=myClusterSideBar,RowAxisColors=1,hclustfun=hclust, distfun = function(x) as.dist(1 - cor(t(x), use = "pa")), method = "average", scale = "row", Colv = NA)


#get color name
clusterColor <- lapply(myClusterSideBar, color.id)
clusterColor <- lapply(clusterColor, `[[`,1)
clusterColor <- data.frame(unlist(clusterColor), stringsAsFactors = FALSE)

foo <- cbind(data.frame(mycl), clusterColor)
#add protein column to data frame for merging
foo$protein_ID <- rownames(foo)


clade.prot <- merge(data.frame(NSAF_matrix), foo, by = "row.names")
clade.prot <- clade.prot[,-1]


colnames(clade.prot)[8] <- "ClusterColor"
colnames(clade.prot)[7] <- "ClusterID"

#prepare table for ggplotting

clade.prot_STACKED <- tidyr::gather(clade.prot, "day", "logFC_NSAF", 1:6)
clade.prot_STACKED$day <- gsub("D_","", clade.prot_STACKED$day)
clade.prot_STACKED$day <- gsub("_logFC_NSAF","", clade.prot_STACKED$day)
clade.prot_STACKED$day <- as.numeric(as.character(clade.prot_STACKED$day))

#plot facetted trend line plots
ggplot(clade.prot_STACKED,aes(day,logFC_NSAF)) + geom_line(aes(group = protein_ID, color  = "red"))  + geom_boxplot(aes(group = day, alpha = 0.01)) + scale_x_continuous(breaks = c(3,5,7,9,11,13), labels = c(3,5,7,9,11,13)) + facet_wrap(~ClusterColor) + theme_bw()



###explore genes within clades

#indianred3
IR3_prots <- data.frame(clade.prot[which(clade.prot$ClusterColor == "indianred3"),"protein_ID"])
colnames(IR3_prots) <- "protein_ID"

#read in uniprot mappings
uniprot <- read.csv("/Volumes/web/metacarcinus/Cgigas/all_giga-uniprot-blastP-out.nopipe.annotations.tab", sep ="\t", header = FALSE, stringsAsFactors = FALSE)
colnames(uniprot) <- c("protein_ID","Entry", "Entry_name", "perc_ident_match", "align_len", "num_mismatch", "num_gaps","querStart", "querEnd", "subjStart", "subjEnd", "evalue", "bitscore","Entry.1","Entry_name.1", "Protein_names", "Gene_names", "Organism", "Protein_length","Pathway", "GO_bp", "GO","GO_IDs", "Protein_fams")
uniprot_mapped <- data.frame(uniprot[which(uniprot$Entry != "unmapped"),], stringsAsFactors = FALSE)
class(uniprot_mapped$evalue) <- "numeric"
all_giga_prots_mapped <- uniprot_mapped[which(uniprot_mapped$evalue <= 10^-10),]


IR3_prot_GO <- merge(IR3_prots, all_giga_prots_mapped[,c(1,23)], by = "protein_ID",all.x = TRUE)

IR3_GO_IDs <- data.frame()

#loop through each line of data frame containing one column with "protein_ID" and one column with list of "GO_IDs" separated by ";"
for (i in 1:nrow(IR3_prot_GO)){
  #create a row that lists all GO IDs associated with one protein listed across rows, where each different GO ID gets put in it's own column
  IR3_GO_ID_row <- data.frame(t(data.frame(strsplit(as.character(IR3_prot_GO$GO_IDs[i]),'; ', fixed = TRUE))), stringsAsFactors = FALSE)
  #add each row created in the line above to the empty data frame
  #this will add NAs in columns where less terms exist; e.g. if a protein only has two terms and another has 10, there will be 8 NAs added to the row for the protein with 2 terms
  IR3_GO_IDs <- rbind.fill(IR3_GO_IDs,IR3_GO_ID_row)
}

#add protein IDs back to GO IDs
IR3_GO_IDs <- cbind(IR3_prot_GO[,"protein_ID"], IR3_GO_IDs)
#this results in a table with protein ID listed in one column and each next column contains a GO ID that was listed with the protein in Uniprot DB.

#reshape data so that all GO ID columns are gathered in one column called "GO" 
STACKED_IR3_GO_IDs <- tidyr::gather(IR3_GO_IDs,"protein_ID","GO", 2:ncol(IR3_GO_IDs))
#exlude middle column which just contains the string "protein_ID" in each row
STACKED_IR3_GO_IDs <- STACKED_IR3_GO_IDs[,c(1,3)]
#remove duplicate rows
STACKED_IR3_GO_IDs <- unique(STACKED_IR3_GO_IDs)
colnames(STACKED_IR3_GO_IDs)[1] <- "protein_ID"
#remove any rows where GO column has NA value. 
STACKED_IR3_GO_IDs <- STACKED_IR3_GO_IDs[which(!is.na(STACKED_IR3_GO_IDs$GO)),]
#this resulting data frame has two columns "protein_ID" and "GO"

