#Determining which stats methods lead to selection of proteins that are changing based on silo

#load libraries
library(plyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(gtools)

#read in uniprot mapping
#uniprot <- read.csv("/Volumes/web/metacarcinus/Cgigas/all_giga-uniprot-blastP-out.nopipe.annotations.tab", sep ="\t", header = FALSE, stringsAsFactors = FALSE)
#select only some uniprot columns
#colnames(uniprot) <- c("protein_ID","Entry", "Entry_name", "perc_ident_match", "align_len", "num_mismatch", "num_gaps","querStart", "querEnd", "subjStart", "subjEnd", "evalue", "bitscore","Entry.1","Entry_name.1", "Protein_names", "Gene_names", "Organism", "Protein_length","Pathway", "GO_bp", "GO","GO_IDs", "Protein_fams")

#read in same day log FC and pval data
sameday_logFC_pval <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/TotNumSpecRatio_FC_Pval/sumNUMSPECSTOT_plus1_AllSilos_ratioFC_logFC_pval_DAYSCOMPARED.csv", stringsAsFactors = FALSE)
colnames(sameday_logFC_pval)[1] <- "protein_ID"
sameday_logFC_pval$protein_ID <- gsub("\\|","\\.",sameday_logFC_pval$protein_ID)

#####select only proteins with adj Chi sq. pvalue <= 0.1####

#create a list of all column names with adj.Chisq.pval 
adjChiSqpvalColumns <- colnames(sameday_logFC_pval[grep("adj.ChiSq.pval_23",colnames(sameday_logFC_pval))])

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
#[1] 227

#make a data frame of just unique proteins so we can select these from the foldchange/pvalue data
all_sig0.1_pro <- unique(all_sig_pro)
colnames(all_sig0.1_pro)[1]<- "protein_ID"
all_sig0.1_pro$Chi <- "PropTest"

#read in ASCA data
ASCA_tempdata <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_all_proteins_avgADJNSAF/Silo2vs3/ASCA_SiloAffectedProteins_noDay0or15.csv", stringsAsFactors = FALSE)
colnames(ASCA_tempdata)[1]<- "protein_ID"
ASCA_tempdata$ASCA <- "ASCA"


#make list of ChiSq, ASCA, and Clustering proteins to pull out
all_sig0.1_ASCA_pro <- merge(all_sig0.1_pro, ASCA_tempdata, by = "protein_ID", all = TRUE)
all_sig0.1_ASCA_pro$method <- paste(all_sig0.1_ASCA_pro$Chi, all_sig0.1_ASCA_pro$ASCA, sep = "")
all_sig0.1_ASCA_pro$method <- gsub("NA","",all_sig0.1_ASCA_pro$method)
data.frame(table(all_sig0.1_ASCA_pro$method))


#read in list of proteins that mapped to uniprot db
proteins_evalpass <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/proteins_evalpass.csv", stringsAsFactors = FALSE)
#exclude proteins that don't pass eval cut-off
all_sig0.1_ASCA_pro <- all_sig0.1_ASCA_pro[which(all_sig0.1_ASCA_pro$protein_ID %in% proteins_evalpass[,1]),]
nrow(all_sig0.1_ASCA_pro)
#172
write.csv(all_sig0.1_ASCA_pro, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/allS2S3_sig0.1_ASCA_pro.csv", quote = FALSE, row.names = FALSE)


#plotting abundance with avg NSAF values
avgNSAF_data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/silo3and2_nozerovals_AVGs.csv", stringsAsFactors = FALSE)
#exclude day 0 and 15
avgNSAF_data <- avgNSAF_data[which(avgNSAF_data$day !=0 & avgNSAF_data$day !=15 ),]

STACKED_avgNSAF_data <- gather(avgNSAF_data, "protein_ID", "avgNSAF", 4:ncol(avgNSAF_data))
all_sig0.1_ASCA_pro_avgNSAF <- merge(all_sig0.1_ASCA_pro, STACKED_avgNSAF_data, by = "protein_ID")

#plot avg NSAF abundances for ASCA and order facets by rank
d <- ggplot(all_sig0.1_ASCA_pro_avgNSAF[grep("ASCA", all_sig0.1_ASCA_pro_avgNSAF$method),], aes(x = day, y = avgNSAF, color = silo)) + geom_line() + facet_wrap(~protein_ID, scale = "free") + ggtitle("Avg NSAF Abundance for proteins selected by ASCA")
ggsave("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/silo2v3/plots/ASCA_selects_avgNSAF.pdf", d, width = 19, height = 13 )
d <- ggplot(all_sig0.1_ASCA_pro_avgNSAF[grep("PropTest", all_sig0.1_ASCA_pro_avgNSAF$method),], aes(x = day, y = avgNSAF, color = silo)) + geom_line() + facet_wrap(~protein_ID, scale = "free") + ggtitle("Avg NSAF Abundance for proteins selected by ChiSq. prop. test pval 0.1")
ggsave("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/silo2v3/plots/chisqpval0.1_selects_avgNSAF.pdf", d, width = 19, height = 13 )

#order proteins by method and turn cluster off; or make 3 heatmaps
###plot heatmap of selected proteins (including unmapped)
all_sig0.1_ASCA_pro_avgNSAF$daysilo <- paste("Day",all_sig0.1_ASCA_pro_avgNSAF$day, "silo",substr(all_sig0.1_ASCA_pro_avgNSAF$silo,2,2), sep = "_")

data4heatmap <- tidyr::spread(all_sig0.1_ASCA_pro_avgNSAF[,c("method","protein_ID","avgNSAF", "daysilo")], daysilo, avgNSAF)
data4heatmap <- data4heatmap[order(data4heatmap$method),]
data4heatmap3 <- data4heatmap[,mixedsort(colnames(data4heatmap[,-c(1:2)]))]
rownames(data4heatmap3) <- data4heatmap[,2]

#make data frames to store heatmap color values
#columns are days/temp; color 23C blue, 29C red, each day a different shade of gray
ColSideColors<-cbind(Temp=c(rep(c("darkorange","blueviolet"),6)), Day=c(rep("#D9D9D9",2),rep("#BDBDBD",2),rep("#969696",2),rep("#737373",2),rep("#525252",2),rep("#252525",2)))

#rows are proteins ordered by method; color by method
#first figure out the frequency of each method
table(data4heatmap$method)

RowSideColors <- c(rep("lightslateblue",46),rep("palegreen",105),rep("chocolate4",20))

heatmap3(data4heatmap3,Colv = NA, cexRow = 0.1, cexCol = 0.5, RowSideColors=RowSideColors,ColSideColors = ColSideColors,RowAxisColors=1,ColAxisColors=1)
#code below replaces original legend with methods legend, would be nice to figure out how to keep both

png("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/silo2v3/plots/heatmap_NSAF_ColorByMethod.png",width = 1000, height = 1000,bg="transparent")
heatmap3(as.matrix(data4heatmap3),Colv = NA, cexRow = 0.5, cexCol = 1, RowSideColors=RowSideColors,ColSideColors = ColSideColors,RowAxisColors=1,ColAxisColors=1,legendfun=function() showLegend(legend=c("ASCA", "Clustering", "Proportions_Test", "Proportions_Test_ASCA","All methods"),col=c("lightslateblue","lightseagreen","gold","palegreen","maroon3","chocolate4"), title = "Method", cex = 1.5))
dev.off()

#can't figure out how to fix the resolution of the axis labels



