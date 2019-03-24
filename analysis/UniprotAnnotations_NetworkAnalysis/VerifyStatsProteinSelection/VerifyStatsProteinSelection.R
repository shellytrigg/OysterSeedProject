#Determining which stats methods lead to selection of proteins that are changing based on temperature

#load libraries
library(plyr)
library(tidyr)
library(reshape2)
library(ggplot2)


#read in uniprot mapping
#uniprot <- read.csv("/Volumes/web/metacarcinus/Cgigas/all_giga-uniprot-blastP-out.nopipe.annotations.tab", sep ="\t", header = FALSE, stringsAsFactors = FALSE)
#select only some uniprot columns
#colnames(uniprot) <- c("protein_ID","Entry", "Entry_name", "perc_ident_match", "align_len", "num_mismatch", "num_gaps","querStart", "querEnd", "subjStart", "subjEnd", "evalue", "bitscore","Entry.1","Entry_name.1", "Protein_names", "Gene_names", "Organism", "Protein_length","Pathway", "GO_bp", "GO","GO_IDs", "Protein_fams")

#read in same day log FC and pval data
sameday_logFC_pval <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/TotNumSpecRatio_FC_Pval/sumNUMSPECSTOT_plus1_ratioFC_logFC_pval_DAYSCOMPARED.csv", stringsAsFactors = FALSE)
colnames(sameday_logFC_pval)[1] <- "protein_ID"

#####select only proteins with adj Chi sq. pvalue <= 0.1####

#create a list of all column names with adj.Chisq.pval 
adjChiSqpvalColumns <- colnames(sameday_logFC_pval[grep("adj.ChiSq.pval",colnames(sameday_logFC_pval))])

#build a list of protiens with adj Chi sq. pvalue <= 0.1
#create empty data frame the loop will add too
all_sig_pro <- data.frame()
for (i in 1:length(adjChiSqpvalColumns)){ # for each name in adj.Chisq.pval column name list
  column <- adjChiSqpvalColumns[i] # create a variable for indexed column name
  #make a data frame containing protein IDs for all proteins in indexed column that have adj.Chisq.pval <=0.1
  sig_pro <- data.frame(sameday_logFC_pval[which(sameday_logFC_pval[,column] <= 0.1),1],stringsAsFactors = FALSE)
  #iteratively add protein lists to initial data frame
  all_sig_pro <- rbind(all_sig_pro, sig_pro)
}

#count how many unique proteins are in the list of proteins with adj Chi sq. pvalue <= 0.1
nrow(unique(all_sig_pro))
#[1] 153

#make a data frame of just unique proteins so we can select these from the foldchange/pvalue data
all_sig0.1_pro <- unique(all_sig_pro)
colnames(all_sig0.1_pro)[1]<- "protein_ID"
all_sig0.1_pro$Chi <- "PropTest"

#read in ASCA data
ASCA_tempdata <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_all_proteins_avgADJNSAF/ASCA_avgNSAFvals_allProteins_noDay0/ASCA_TempAffectedProteins_loadings.csv", stringsAsFactors = FALSE)
colnames(ASCA_tempdata)[1]<- "protein_ID"
ASCA_tempdata$ASCA <- "ASCA"
#rank ASCA proteins by their magnitude loadings value and make a column that contains this info
ASCA_tempdata$rank<- NA
ASCA_tempdata <- ASCA_tempdata[order(ASCA_tempdata$PC1loadings),]
ASCA_tempdata$magloading <- abs(ASCA_tempdata$PC1loadings)
ASCA_tempdata <- ASCA_tempdata[order(ASCA_tempdata$magloading, decreasing = TRUE),]
ASCA_tempdata$rank <- 1:nrow(ASCA_tempdata)


#read in clustering data
clust_data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/250unique-prot.csv ", stringsAsFactors = FALSE)
clust_data <- data.frame(unique(clust_data[,"ID"]),stringsAsFactors = FALSE)
colnames(clust_data)[1] <- "protein_ID"
clust_data$clust <- "clustering"

#make list of ChiSq, ASCA, and Clustering proteins to pull out
all_sig0.1_ASCA_clust_pro <- merge(all_sig0.1_pro, ASCA_tempdata, by = "protein_ID", all = TRUE)
all_sig0.1_ASCA_clust_pro <- merge(all_sig0.1_ASCA_clust_pro, clust_data, by = "protein_ID", all = TRUE)
all_sig0.1_ASCA_clust_pro$method <- paste(all_sig0.1_ASCA_clust_pro$Chi, all_sig0.1_ASCA_clust_pro$ASCA,all_sig0.1_ASCA_clust_pro$clust, sep = "")
all_sig0.1_ASCA_clust_pro$method <- gsub("NA","",all_sig0.1_ASCA_clust_pro$method)



#read in list of proteins that mapped to uniprot db
proteins_evalpass <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/proteins_evalpass.csv", stringsAsFactors = FALSE)
#exclude proteins that don't pass eval cut-off
all_sig0.1_ASCA_clust_pro <- all_sig0.1_ASCA_clust_pro[which(all_sig0.1_ASCA_clust_pro$protein_ID %in% proteins_evalpass[,1]),]
nrow(all_sig0.1_ASCA_clust_pro)
#132
write.csv(all_sig0.1_ASCA_clust_pro, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/all_sig0.1_ASCA_clust_pro.csv", quote = FALSE, row.names = FALSE)

###look at total spectral data

TOTSpecNum_data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/TotNumSpecData.csv", stringsAsFactors = FALSE)
#exclude day 0
TOTSpecNum_data <- TOTSpecNum_data[-1,]
STACKED_TOTSpecNum_data <- gather(TOTSpecNum_data, "protein_ID", "TotNumSpec", 5:ncol(TOTSpecNum_data))
#merge with all_sig0.1_ASCA_clust_pro
all_sig0.1_ASCA_clust_pro_totnumspec <- merge(all_sig0.1_ASCA_clust_pro, STACKED_TOTSpecNum_data, by = "protein_ID")
all_sig0.1_ASCA_clust_pro_totnumspec$temp <- as.character(all_sig0.1_ASCA_clust_pro_totnumspec$temp)


#plot abundances for ASCA and order facets by rank
d <- ggplot(all_sig0.1_ASCA_clust_pro_totnumspec[grep("ASCA", all_sig0.1_ASCA_clust_pro_totnumspec$method),], aes(x = day, y = TotNumSpec, color = temp)) + geom_line() + facet_wrap(~rank, scale = "free") + ggtitle("Total Spectral Abundance for proteins selected by ASCA (ordered by rank)")
ggsave("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection//plots/orderedASCA_selects_TotNumSpec.pdf", d, width = 19, height = 13 )
d <- ggplot(all_sig0.1_ASCA_clust_pro_totnumspec[grep("ASCA", all_sig0.1_ASCA_clust_pro_totnumspec$method),], aes(x = day, y = TotNumSpec, color = temp)) + geom_line() + facet_wrap(~protein_ID, scale = "free") + ggtitle("Total Spectral Abundance for proteins selected by ASCA")
ggsave("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/ASCA_selects_TotNumSpec.pdf", d, width = 19, height = 13 )
d <- ggplot(all_sig0.1_ASCA_clust_pro_totnumspec[grep("PropTest", all_sig0.1_ASCA_clust_pro_totnumspec$method),], aes(x = day, y = TotNumSpec, color = temp)) + geom_line() + facet_wrap(~protein_ID, scale = "free") + ggtitle("Total Spectral Abundance for proteins selected by ChiSq prop. test pval 0.1") 
ggsave("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/chisqpval0.1_selects_TotNumSpec.pdf", d, width = 19, height = 13 )
d <- ggplot(all_sig0.1_ASCA_clust_pro_totnumspec[grep("clustering", all_sig0.1_ASCA_clust_pro_totnumspec$method),], aes(x = day, y = TotNumSpec, color = temp)) + geom_line() + facet_wrap(~protein_ID, scale = "free") + ggtitle("Total Spectral Abundance for proteins selected by clustering")
ggsave("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/clustering_selects_TotNumSpec.pdf", d, width = 19, height = 13 )

#plotting abundance with avg NSAF values
avgNSAF_data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/silo3and9_nozerovals_AVGs.csv", stringsAsFactors = FALSE)
#exclude day 0
avgNSAF_data <- avgNSAF_data[-1,]

STACKED_avgNSAF_data <- gather(avgNSAF_data, "protein_ID", "avgNSAF", 3:ncol(avgNSAF_data))
all_sig0.1_ASCA_clust_pro_avgNSAF <- merge(all_sig0.1_ASCA_clust_pro, STACKED_avgNSAF_data, by = "protein_ID")
#convert temp from integer to character so plot will have two distinctly colored lines
all_sig0.1_ASCA_clust_pro_avgNSAF$temp <- as.character(all_sig0.1_ASCA_clust_pro_avgNSAF$temp)
#plot avg NSAF abundances for ASCA and order facets by rank
d <- ggplot(all_sig0.1_ASCA_clust_pro_avgNSAF[grep("ASCA", all_sig0.1_ASCA_clust_pro_avgNSAF$method),], aes(x = day, y = avgNSAF, color = temp)) + geom_line() + facet_wrap(~rank, scale = "free") + ggtitle("Avg NSAF Abundance for proteins selected by ASCA (facets ordered by rank)")
ggsave("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/orderedASCA_selects_avgNSAF.pdf", d, width = 19, height = 13 )
d <- ggplot(all_sig0.1_ASCA_clust_pro_avgNSAF[grep("ASCA", all_sig0.1_ASCA_clust_pro_avgNSAF$method),], aes(x = day, y = avgNSAF, color = temp)) + geom_line() + facet_wrap(~protein_ID, scale = "free") + ggtitle("Avg NSAF Abundance for proteins selected by ASCA")
ggsave("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/ASCA_selects_avgNSAF.pdf", d, width = 19, height = 13 )
d <- ggplot(all_sig0.1_ASCA_clust_pro_avgNSAF[grep("PropTest", all_sig0.1_ASCA_clust_pro_avgNSAF$method),], aes(x = day, y = avgNSAF, color = temp)) + geom_line() + facet_wrap(~protein_ID, scale = "free") + ggtitle("Avg NSAF Abundance for proteins selected by ChiSq. prop. test pval 0.1")
ggsave("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/chisqpval0.1_selects_avgNSAF.pdf", d, width = 19, height = 13 )
d <- ggplot(all_sig0.1_ASCA_clust_pro_avgNSAF[grep("clustering", all_sig0.1_ASCA_clust_pro_avgNSAF$method),], aes(x = day, y = avgNSAF, color = temp)) + geom_line() + facet_wrap(~protein_ID, scale = "free") + ggtitle("Avg NSAF Abundance for proteins selected by clustering")
ggsave("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/clustering_selects_avgNSAF.pdf", d, width = 19, height = 13 )

######NSAF values seem to show a different pattern than TotNumSpec values in about half the proteins, so it's possible the clustering and ASCA are asking different questions than the prop test.
####Proteins identified by clustering have a 93% overlap with proteins identified by other methods, only one protein was uniquely identified.
###
#can order proteins by method and turn cluster off; or make 3 heatmaps


###plot heatmap of selected proteins (including unmapped)
all_sig0.1_ASCA_clust_pro_avgNSAF$daytemp <- paste("Day",all_sig0.1_ASCA_clust_pro_avgNSAF$day, all_sig0.1_ASCA_clust_pro_avgNSAF$temp, "C", sep = "_")

data4heatmap <- tidyr::spread(all_sig0.1_ASCA_clust_pro_avgNSAF[,c("method","protein_ID","avgNSAF", "daytemp")], daytemp, avgNSAF)
data4heatmap <- data4heatmap[order(data4heatmap$method),]
data4heatmap3 <- data4heatmap[,mixedsort(colnames(data4heatmap[,-c(1:2)]))]
rownames(data4heatmap3) <- data4heatmap[,2]

#make data frames to store heatmap color values
#columns are days/temp; color 23C blue, 29C red, each day a different shade of gray
ColSideColors<-cbind(Temp=c(rep(c("steelblue2","red"),6)), Day=c(rep("#D9D9D9",2),rep("#BDBDBD",2),rep("#969696",2),rep("#737373",2),rep("#525252",2),rep("#252525",2)))

#rows are proteins ordered by method; color by method
#first figure out the frequency of each method
table(data4heatmap$method)

RowSideColors <- c(rep("lightslateblue",101),rep("lightseagreen",8),rep("gold",8),rep("palegreen",129), rep("maroon3",17),rep("chocolate4",17))

heatmap3(data4heatmap3,Colv = NA, cexRow = 0.1, cexCol = 0.5, RowSideColors=RowSideColors,ColSideColors = ColSideColors,RowAxisColors=1,ColAxisColors=1)
#code below replaces original legend with methods legend, would be nice to figure out how to keep both

png("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/heatmap_NSAF_ColorByMethod.png",width = 1000, height = 1000,bg="transparent")
heatmap3(data4heatmap3,Colv = NA, cexRow = 0.5, cexCol = 1, RowSideColors=RowSideColors,ColSideColors = ColSideColors,RowAxisColors=1,ColAxisColors=1,legendfun=function() showLegend(legend=c("ASCA","ASCA_Clustering", "Clustering", "Proportions_Test", "Proportions_Test_ASCA","All methods"),col=c("lightslateblue","lightseagreen","gold","palegreen","maroon3","chocolate4"), title = "Method", cex = 1.5))
dev.off()

#can't figure out how to fix the resolution of the axis labels



