#Determining which stats methods lead to selection of proteins that are changing based on temperature

#load libraries
library(plyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(gtools)
library(heatmap3)
library(devtools)
#install_github("js229/Vennerable")
library(Vennerable)
library(colorRamps)
library(plotrix)


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
#[1] 163

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
#clust_data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/250unique-prot.csv", stringsAsFactors = FALSE)
clust_data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/kmeans/Silo3_and_9/silo3_9-NSAF_euclidean/150unique-prot.csv", stringsAsFactors = FALSE)
clust_data <- data.frame(unique(clust_data[,"ID"]),stringsAsFactors = FALSE)
colnames(clust_data)[1] <- "protein_ID"
clust_data$clust <- "clustering"

#make list of ChiSq, ASCA, and Clustering proteins to pull out
all_sig0.1_ASCA_clust_pro <- merge(all_sig0.1_pro, ASCA_tempdata, by = "protein_ID", all = TRUE)
all_sig0.1_ASCA_clust_pro <- merge(all_sig0.1_ASCA_clust_pro, clust_data, by = "protein_ID", all = TRUE)
all_sig0.1_ASCA_clust_pro$method <- paste(all_sig0.1_ASCA_clust_pro$Chi, all_sig0.1_ASCA_clust_pro$ASCA,all_sig0.1_ASCA_clust_pro$clust, sep = "")
all_sig0.1_ASCA_clust_pro$method <- gsub("NA","",all_sig0.1_ASCA_clust_pro$method)
nrow(all_sig0.1_ASCA_clust_pro)
#341
write.csv(all_sig0.1_ASCA_clust_pro, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/all_sig0.1_ASCA_clust_pro.csv", quote = FALSE, row.names = FALSE)


#read in list of proteins that mapped to uniprot db
proteins_evalpass <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/proteins_evalpass.csv", stringsAsFactors = FALSE)
#exclude proteins that don't pass eval cut-off
all_sig0.1_ASCA_clust_pro_evalpass <- all_sig0.1_ASCA_clust_pro[which(all_sig0.1_ASCA_clust_pro$protein_ID %in% proteins_evalpass[,1]),]
nrow(all_sig0.1_ASCA_clust_pro_evalpass)
#284
write.csv(all_sig0.1_ASCA_clust_pro_evalpass, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/all_sig0.1_ASCA_clust_pro_evalpass.csv", quote = FALSE, row.names = FALSE)


###make venn
selects <- list(ASCA = all_sig0.1_ASCA_clust_pro[grep("ASCA",all_sig0.1_ASCA_clust_pro$method),"protein_ID"], Clustering = all_sig0.1_ASCA_clust_pro[grep("clustering",all_sig0.1_ASCA_clust_pro$method),"protein_ID"], Prop_test = all_sig0.1_ASCA_clust_pro[grep("Prop",all_sig0.1_ASCA_clust_pro$method),"protein_ID"])
V <- Venn(selects)
V
pdf("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/Venn_NSAF_Allprots_ColorByMethod.pdf",width = 10, height = 10)
plot(V,doWeights = TRUE)
dev.off()


selects_evalpass <- list(ASCA = all_sig0.1_ASCA_clust_pro_evalpass[grep("ASCA",all_sig0.1_ASCA_clust_pro_evalpass$method),"protein_ID"], Clustering = all_sig0.1_ASCA_clust_pro_evalpass[grep("clustering",all_sig0.1_ASCA_clust_pro_evalpass$method),"protein_ID"], Prop_test = all_sig0.1_ASCA_clust_pro_evalpass[grep("Prop",all_sig0.1_ASCA_clust_pro_evalpass$method),"protein_ID"])
Veval <- Venn(selects_evalpass)
Veval
pdf("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/Venn_NSAF_evalpassProts_ColorByMethod.pdf",width = 10, height = 10)
plot(Veval,doWeights = TRUE)
dev.off()


###########
#PLotting NSAF values of selected proteins for each method; this is to get an idea of how well our methods are selecting proteins with abundances changes related to temperature
##########

#read in avg NSAF data
avgNSAF_data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/silo3and9_nozerovals_AVGs.csv", stringsAsFactors = FALSE)
#exclude day 0
avgNSAF_data <- avgNSAF_data[-1,]
#format data into long format
STACKED_avgNSAF_data <- gather(avgNSAF_data, "protein_ID", "avgNSAF", 3:ncol(avgNSAF_data))
#subset NSAF data by merging it with selected protein data
all_sig0.1_ASCA_clust_pro_avgNSAF <- merge(all_sig0.1_ASCA_clust_pro, STACKED_avgNSAF_data, by = "protein_ID")
#convert temp from integer to character so plot will have two distinctly colored lines
all_sig0.1_ASCA_clust_pro_avgNSAF$temp <- as.character(all_sig0.1_ASCA_clust_pro_avgNSAF$temp)

###CREATE HEATMAP OF ALL SELECTED PROTEINS
###ordered by day
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
RowSideColors <- c(rep("lightgoldenrod1",70),rep("steelblue1",39),rep("mediumpurple1",69),rep("mediumaquamarine",127), rep("coral2",9),rep("olivedrab1",25), rep("sandybrown",2))
heatmap3(data4heatmap3,Colv = NA, cexRow = 0.1, cexCol = 0.5, RowSideColors=RowSideColors,ColSideColors = ColSideColors,RowAxisColors=1,ColAxisColors=1)
#code below replaces original legend with methods legend, would be nice to figure out how to keep both
pdf("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/heatmap_NSAF_noD0_ColorByMethod.pdf",width = 10, height = 10)
heatmap3(data4heatmap3,Colv = NA, cexRow = 0.1, cexCol = 0.5, RowSideColors=RowSideColors,ColSideColors = ColSideColors,RowAxisColors=1,ColAxisColors=1,legendfun=function() showLegend(legend=c("ASCA","ASCA_Clustering", "Clustering", "Proportions_Test", "Proportions_Test_ASCA","All methods"),col=c("lightgoldenrod1","steelblue1","mediumpurple1","mediumaquamarine","coral2","olivedrab1","sandybrown"), title = "Method", cex = 0.8))
dev.off()

###ordered by temperature then day
data4heatmap3_orderT <- data4heatmap[,c(7,9,11,13,3,5,8,10,12,14,4,6)]
rownames(data4heatmap3_orderT) <- rownames(data4heatmap3)
ColSideColors_orderT <- ColSideColors[order(ColSideColors[,1], decreasing = TRUE),]
#pdf("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/heatmap_NSAF_noD0_ordT_ColorByMethod.pdf",width = 10, height = 10)
#heatmap3(data4heatmap3_orderT,Colv = NA, cexRow = 0.1, cexCol = 0.5, RowSideColors=RowSideColors,ColSideColors = ColSideColors_orderT,RowAxisColors=1,ColAxisColors=1,legendfun=function() showLegend(legend=c("ASCA","ASCA_Clustering", "Clustering", "Proportions_Test", "Proportions_Test_ASCA","All methods"),col=c("lightgoldenrod1","steelblue1","mediumpurple1","mediumaquamarine","coral2","olivedrab1","sandybrown"), title = "Method", cex = 0.8))
pdf("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/heatmap_NSAF_noD0_ordT_ColorByMethod_clustByClade.pdf",width = 10, height = 10)
heatmap3(data4heatmap3_orderT,cexRow = 0.1, cexCol = 0.5, RowSideColors=RowSideColors,ColSideColors = ColSideColors_orderT,RowAxisColors=1,ColAxisColors=1,legendfun=function() showLegend(legend=c("ASCA","ASCA_Clustering", "Clustering","Clustering_Proportions_Test", "Proportions_Test", "Proportions_Test_ASCA","All methods"),col=c("lightgoldenrod1","steelblue1","mediumpurple1","sandybrown","mediumaquamarine","coral2","olivedrab1"), title = "Method", cex = 0.8),hclustfun=hclust, distfun = function(x) as.dist(1 - cor(t(x), use = "pa")), method = "average", scale = "row", Colv=NA)
dev.off()

#extracting clades
heatmap3(data4heatmap3_orderT,cexRow = 0.1, cexCol = 0.5, RowSideColors=RowSideColors,ColSideColors = ColSideColors_orderT,RowAxisColors=1,ColAxisColors=1,legendfun=function() showLegend(legend=c("ASCA","ASCA_Clustering", "Clustering", "Proportions_Test", "Proportions_Test_ASCA","All methods"),col=c("lightgoldenrod1","steelblue1","mediumpurple1","mediumaquamarine","coral2","olivedrab1","sandybrown"), title = "Method", cex = 0.8),hclustfun=hclust, distfun = function(x) as.dist(1 - cor(t(x), use = "pa")), method = "average", scale = "row", Colv=NA)
hm <- as.dist(1-cor(t(data4heatmap3_orderT), use="pa"))
hm2 <- hclust(hm, method = 'average')
plot(hm2)

#Identify clades in heatmap
# define some clusters
mycl <- cutree(hm2, h=0.7)
clusterCols <- colorRamps::primary.colors(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
pdf("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/heatmap_NSAF_noD0_ordT_ColorByMethodAndClade.pdf",width = 10, height = 10)
heatmap3(data4heatmap3_orderT,cexRow = 0.5, cexCol = 0.5, RowSideColors=myClusterSideBar,ColSideColors = ColSideColors_orderT,RowAxisColors=1,ColAxisColors=1,legendfun=function() showLegend(legend=c("ASCA","ASCA_Clustering", "Clustering", "Proportions_Test", "Proportions_Test_ASCA","All methods"),col=c("lightgoldenrod1","steelblue1","mediumpurple1","mediumaquamarine","coral2","olivedrab1","sandybrown"), title = "Method", cex = 0.8),hclustfun=hclust, distfun = function(x) as.dist(1 - cor(t(x), use = "pa")), method = "average", scale = "row", Colv=NA)
dev.off()

#get color name
clusterColor <- lapply(myClusterSideBar, color.id)
clusterColor <- lapply(clusterColor, `[[`,1)
clusterColor <- data.frame(unlist(clusterColor), stringsAsFactors = FALSE)



foo <- cbind(data4heatmap[,2],data4heatmap3_orderT, clusterID=mycl, clusterColor)
#foo <- cbind(data4heatmap3_orderT, clusterID=mycl)

clade.prot <- foo[hm2$order,]
colnames(clade.prot)[1] <- "protein_ID"
colnames(clade.prot)[15] <- "ClusterColor"

write.csv(clade.prot, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/ASCA_Clust_ChiSq_data4heatmap3_orderT.csv", row.names= FALSE, quote = FALSE)


###CREATE HEATMAP OF ASCA SELECTED PROTEINS
ASCAdata4heatmap <- data4heatmap[grep("ASCA", data4heatmap$method),]

#ordered by day
ASCAdata4heatmap3 <- ASCAdata4heatmap[,mixedsort(colnames(ASCAdata4heatmap[,-c(1:2)]))]
rownames(ASCAdata4heatmap3) <- ASCAdata4heatmap[,2]
#make data frames to store heatmap color values
#columns are days/temp; color 23C blue, 29C red, each day a different shade of gray
ColSideColors<-cbind(Temp=c(rep(c("steelblue2","red"),6)), Day=c(rep("#D9D9D9",2),rep("#BDBDBD",2),rep("#969696",2),rep("#737373",2),rep("#525252",2),rep("#252525",2)))
pdf("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/heatmap_ASCA_NSAF_noD0_ColorByMethod.pdf",width = 10, height = 10)
heatmap3(ASCAdata4heatmap3,Colv = NA, cexRow = 0.1, cexCol = 0.5,ColSideColors = ColSideColors,ColAxisColors=1)
dev.off()
#ordered by temp then day
ASCAdata4heatmap3_orderT <- ASCAdata4heatmap[,c(7,9,11,13,3,5,8,10,12,14,4,6)]
ColSideColors_orderT <- ColSideColors[order(ColSideColors[,1], decreasing = TRUE),]
pdf("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/heatmap_ASCA_NSAF_noD0_ordT_ColorByMethod.pdf",width = 10, height = 10)
heatmap3(ASCAdata4heatmap3_orderT,Colv = NA, cexRow = 0.1, cexCol = 0.5,ColSideColors = ColSideColors_orderT,ColAxisColors=1)
dev.off()

#Abundance plots for ASCA selected proteins
d <- ggplot(all_sig0.1_ASCA_clust_pro_avgNSAF[grep("ASCA", all_sig0.1_ASCA_clust_pro_avgNSAF$method),], aes(x = day, y = avgNSAF, color = temp)) + geom_line() + facet_wrap(~protein_ID, scale = "free") + scale_x_continuous(breaks = c(3,5,7,9,11,13), labels = c(3,5,7,9,11,13)) +   theme(strip.text.x = element_text(size = 3)) + ggtitle("Avg NSAF Abundance for proteins selected by ASCA")
ggsave("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/ASCA_selects_avgNSAF.pdf", d, width = 19, height = 13 )

###CREATE HEATMAP OF CLUSTERING SELECTED PROTEINS
clusterdata4heatmap <- data4heatmap[grep("cluster", data4heatmap$method),]

#ordered by day
clusterdata4heatmap3 <- clusterdata4heatmap[,mixedsort(colnames(clusterdata4heatmap[,-c(1:2)]))]
rownames(clusterdata4heatmap3) <- clusterdata4heatmap[,2]
#make data frames to store heatmap color values
#columns are days/temp; color 23C blue, 29C red, each day a different shade of gray
ColSideColors<-cbind(Temp=c(rep(c("steelblue2","red"),6)), Day=c(rep("#D9D9D9",2),rep("#BDBDBD",2),rep("#969696",2),rep("#737373",2),rep("#525252",2),rep("#252525",2)))
pdf("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/heatmap_clustering_NSAF_noD0_ColorByMethod.pdf",width = 10, height = 10)
heatmap3(clusterdata4heatmap3,Colv = NA, cexRow = 0.1, cexCol = 0.5,ColSideColors = ColSideColors,ColAxisColors=1)
dev.off()
#ordered by temp then day
clusterdata4heatmap3_orderT <- clusterdata4heatmap[,c(7,9,11,13,3,5,8,10,12,14,4,6)]
ColSideColors_orderT <- ColSideColors[order(ColSideColors[,1], decreasing = TRUE),]
pdf("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/heatmap_clusterign_NSAF_noD0_ordT_ColorByMethod.pdf",width = 10, height = 10)
heatmap3(clusterdata4heatmap3_orderT,Colv = NA, cexRow = 0.1, cexCol = 0.5,ColSideColors = ColSideColors_orderT,ColAxisColors=1)
dev.off()

#Abundance plots for cluster selected proteins
d <- ggplot(all_sig0.1_ASCA_clust_pro_avgNSAF[grep("cluster", all_sig0.1_ASCA_clust_pro_avgNSAF$method),], aes(x = day, y = avgNSAF, color = temp)) + geom_line() + facet_wrap(~protein_ID, scale = "free") + scale_x_continuous(breaks = c(3,5,7,9,11,13), labels = c(3,5,7,9,11,13)) +   theme(strip.text.x = element_text(size = 3)) + ggtitle("Avg NSAF Abundance for proteins selected by clustering")
ggsave("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/clustering_selects_avgNSAF.pdf", d, width = 19, height = 13 )


###CREATE HEATMAP OF PROP TEST SELECTED PROTEINS
### *** Since prop test was on TOTNUMSPEC, I will plot those instead of NSAFs

TOTSpecNum_data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/TotNumSpecData.csv", stringsAsFactors = FALSE)
#exclude day 0
TOTSpecNum_data <- TOTSpecNum_data[-1,]
STACKED_TOTSpecNum_data <- gather(TOTSpecNum_data, "protein_ID", "TotNumSpec", 5:ncol(TOTSpecNum_data))
#merge with all_sig0.1_ASCA_clust_pro
all_sig0.1_ASCA_clust_pro_totnumspec <- merge(all_sig0.1_ASCA_clust_pro, STACKED_TOTSpecNum_data, by = "protein_ID")
all_sig0.1_ASCA_clust_pro_totnumspec$temp <- as.character(all_sig0.1_ASCA_clust_pro_totnumspec$temp)


all_sig0.1_ASCA_clust_pro_totnumspec$daytemp <- paste("Day",all_sig0.1_ASCA_clust_pro_totnumspec$day, all_sig0.1_ASCA_clust_pro_totnumspec$temp, "C", sep = "_")
data4heatmap <- tidyr::spread(all_sig0.1_ASCA_clust_pro_totnumspec[,c("method","protein_ID","TotNumSpec", "daytemp")], daytemp, TotNumSpec)
data4heatmap <- data4heatmap[order(data4heatmap$method),]
data4heatmap3 <- data4heatmap[,mixedsort(colnames(data4heatmap[,-c(1:2)]))]
rownames(data4heatmap3) <- data4heatmap[,2]

proptestdata4heatmap <- data4heatmap[grep("Prop", data4heatmap$method),]

#ordered by day
proptestdata4heatmap3 <- proptestdata4heatmap[,mixedsort(colnames(proptestdata4heatmap[,-c(1:2)]))]
rownames(proptestdata4heatmap3) <- proptestdata4heatmap[,2]
#make data frames to store heatmap color values
#columns are days/temp; color 23C blue, 29C red, each day a different shade of gray
ColSideColors<-cbind(Temp=c(rep(c("steelblue2","red"),6)), Day=c(rep("#D9D9D9",2),rep("#BDBDBD",2),rep("#969696",2),rep("#737373",2),rep("#525252",2),rep("#252525",2)))
pdf("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/heatmap_proptest_NSAF_noD0_ColorByMethod.pdf",width = 10, height = 10)
heatmap3(proptestdata4heatmap3,Colv = NA, cexRow = 0.1, cexCol = 0.5,ColSideColors = ColSideColors,ColAxisColors=1)
dev.off()
#ordered by temp then day
proptestdata4heatmap3_orderT <- proptestdata4heatmap[,c(7,9,11,13,3,5,8,10,12,14,4,6)]
ColSideColors_orderT <- ColSideColors[order(ColSideColors[,1], decreasing = TRUE),]
pdf("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/heatmap_proptest_NSAF_noD0_ordT_ColorByMethod.pdf",width = 10, height = 10)
heatmap3(proptestdata4heatmap3_orderT,Colv = NA, cexRow = 0.1, cexCol = 0.5,ColSideColors = ColSideColors_orderT,ColAxisColors=1)
dev.off()

d <- ggplot(all_sig0.1_ASCA_clust_pro_totnumspec[grep("PropTest", all_sig0.1_ASCA_clust_pro_totnumspec$method),], aes(x = day, y = TotNumSpec, color = temp)) + geom_line() + facet_wrap(~protein_ID, scale = "free") + scale_x_continuous(breaks = c(3,5,7,9,11,13), labels = c(3,5,7,9,11,13)) +   theme(strip.text.x = element_text(size = 3)) + ggtitle("Total Spectral Abundance for proteins selected by ChiSq prop. test pval 0.1")
ggsave("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/plots/chisqpval0.1_selects_TotNumSpec.pdf", d, width = 19, height = 13 )

