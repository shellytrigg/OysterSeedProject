library(heatmap3)
library(reshape2)
library(tidyr)
library(plyr)
library(ontologyIndex)
library(ontologySimilarity)
library(colorRamps)
library(plotrix)


NSAF_for_big_heatmap <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/silo3and9_nozerovals_AVGs.csv", stringsAsFactors = FALSE)
NSAF_for_big_heatmap <- NSAF_for_big_heatmap[,c(1,2,grep("CHOYP", colnames(NSAF_for_big_heatmap)))]
NSAF_for_big_heatmap_t <- t(NSAF_for_big_heatmap[,-c(1:2)])
colnames(NSAF_for_big_heatmap_t) <- paste0("day",NSAF_for_big_heatmap$day,"_", NSAF_for_big_heatmap$temp,"C")


ColSideColors<-cbind(Temp=c("purple",rep(c("steelblue2","red"),6)), Day=c("white",rep("#D9D9D9",2),rep("#BDBDBD",2),rep("#969696",2),rep("#737373",2),rep("#525252",2),rep("#252525",2)))


pdf("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/General_proteome_characterization/NSAF_all_CHOYP_proteins.pdf",width = 10, height = 10)
heatmap3(NSAF_for_big_heatmap_t, cexRow = 0.1, cexCol = 0.8, Colv = NA, method = "average",ColSideColors = ColSideColors,ColAxisColors=1)
dev.off()

jpeg("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/General_proteome_characterization/NSAF_all_CHOYP_proteins.jpg", width = 2000, height = 2000)
heatmap3(NSAF_for_big_heatmap_t, cexRow = 0.1, cexCol = 0.8, Colv = NA, method = "average", ColSideColors = ColSideColors,ColAxisColors=1)
dev.off()




###extracting proteins that show abundance differences in heatmap
NSAF_for_big_heatmap_t_scaled <- data.frame(t(apply(data.frame(NSAF_for_big_heatmap_t),1,scale)))

colnames(NSAF_for_big_heatmap_t_scaled) <- colnames(NSAF_for_big_heatmap_t)



up <- data.frame()
down <- data.frame()

for(i in (1:ncol(NSAF_for_big_heatmap_t_scaled))){
  up_prots <- data.frame(rownames(NSAF_for_big_heatmap_t_scaled[which(NSAF_for_big_heatmap_t_scaled[,i] > 3),]))
  down_prots <- data.frame(rownames(NSAF_for_big_heatmap_t_scaled[which(NSAF_for_big_heatmap_t_scaled[,i] < -2),]))
  
  colnames(up_prots) <- "protein_ID"
  colnames(down_prots) <- "protein_ID"
  
  up_prots$daytemp <- colnames(NSAF_for_big_heatmap_t_scaled)[i]
  down_prots$daytemp <- colnames(NSAF_for_big_heatmap_t_scaled)[i]
  
  up <- rbind(up, up_prots)
  down <- rbind(down, down_prots)
}


table(up$daytemp)
table(down$daytemp)



GO_up_BPs <- all_giga_prots_mapped[which(all_giga_prots_mapped$protein_ID %in% up$protein_ID),]
GO_BPs <- data.frame()

#get GO BP terms
for (i in 1:nrow(GO_up_BPs)){
  #create a row that lists all GO IDs associated with one protein listed across rows, where each different GO ID gets put in it's own column
  all_giga_prots_mapped_GObp_row <- data.frame(t(data.frame(strsplit(as.character(GO_up_BPs$GO_bp[i]),'; ', fixed = TRUE))), stringsAsFactors = FALSE)
  #add each row created in the line above to the empty data frame
  #this will add NAs in columns where less terms exist; e.g. if a protein only has two terms and another has 10, there will be 8 NAs added to the row for the protein with 2 terms
  GO_BPs <- rbind.fill(GO_BPs,all_giga_prots_mapped_GObp_row)
}

#add protein IDs back to GO IDs
GO_BPs <- cbind(GO_up_BPs[,"protein_ID"], GO_BPs)
colnames(GO_BPs)[1] <- "protein_ID"

up_bps <- merge(up,GO_BPs, by = "protein_ID")

up_bps_STACKED <- tidyr::gather(up_bps,"term","GO_term", 3:ncol(up_bps))
up_bps_STACKED$term <- NULL
up_bps_STACKED$GO_ID <- gsub(".* \\[","",up_bps_STACKED$GO_term)
up_bps_STACKED$GO_ID <- gsub("\\]","",up_bps_STACKED$GO_ID)
length(up_bps_STACKED$GO_ID)
#[1] 70720
#[1] 68530???

###remove proteins that don't have GO BP terms
up_bps_STACKED <- up_bps_STACKED[which(!(is.na(up_bps_STACKED$GO_ID))),]
nrow(up_bps_STACKED)
#[1] 2456
##[1] 3611???


## make down BP dataframes

down_bps <- merge(down,GO_BPs, by = "protein_ID")

down_bps_STACKED <- tidyr::gather(down_bps,"term","GO_term", 3:ncol(down_bps))
down_bps_STACKED$term <- NULL
down_bps_STACKED$GO_ID <- gsub(".* \\[","",down_bps_STACKED$GO_term)
down_bps_STACKED$GO_ID <- gsub("\\]","",down_bps_STACKED$GO_ID)
length(down_bps_STACKED$GO_ID)
#[1] 70720
#[1] 68530???

###remove proteins that don't have GO BP terms
down_bps_STACKED <- down_bps_STACKED[which(!(is.na(down_bps_STACKED$GO_ID))),]
nrow(down_bps_STACKED)
#[1] 2456
##[1] 3611???

###remove proteins/term pairs that are not in the OntologyIndex data file
data(go)
data(GO_IC)

up_bps_STACKED <- up_bps_STACKED[which(up_bps_STACKED$GO_ID %in% go$id),]

up_bps_STACKED <- up_bps_STACKED[which(up_bps_STACKED$GO_ID %in% names(GO_IC)),]
nrow(up_bps_STACKED)
#[1] 3415

###make list of proteins and their go IDs like beach example

up_bps_spread_day0 <- tidyr::spread(up_bps_STACKED[grep("day0", up_bps_STACKED$daytemp),], protein_ID, GO_ID)
up_bps_spread_day0 <- up_bps_spread_day0[,-c(1:2)]

up_bps_day0_list <- lapply(as.list(up_bps_spread_day0),function(x) x[!is.na(x)])


#run get sim matrix
sim_matrix <- get_sim_grid(ontology = go, information_content = GO_IC, term_sets = up_bps_day0_list)

#need to create a dissimilarity matrix so need to transform similarity matrix
dist_mat <- max(sim_matrix) - sim_matrix

#plot heatmap
heatmap3(dist_mat, scale = "none")
## most of the data is not similar to each other (mostly blue)
## try a cutoff and re-plot


heatmap3(sim_matrix, scale = "none")


hm <- as.dist(1-cor(t(sim_matrix), use="pa"))

hm2 <- hclust(hm, method = 'complete')
plot(hm2)

#Identify clades in heatmap
# define some clusters

#make a list of proteins with cluster IDs based on where the tree is cut

#for complete clustering
mycl <- cutree(hm2,h=1)

#assign colors to cluster IDs
clusterCols <- colorRamps::primary.colors(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
heatmap3(sim_matrix,cexRow = 0.5, cexCol = 0.5, RowSideColors=myClusterSideBar,RowAxisColors=1,hclustfun=hclust, distfun = function(x) as.dist(1 - cor(t(x), use = "pa")), method = "complete", scale = "none")

#get color name
clusterColor <- lapply(myClusterSideBar, color.id)
clusterColor <- lapply(clusterColor, `[[`,1)
clusterColor <- data.frame(unlist(clusterColor), stringsAsFactors = FALSE)

foo <- cbind(data.frame(mycl), clusterColor)
#add protein column to data frame for merging
foo$protein_ID <- rownames(foo)


View(up_bps_STACKED[which(up_bps_STACKED$protein_ID %in% foo[grep("darkslategray1",foo$unlist.clusterColor.),"protein_ID"]),])




#### try term reduction with GO slims using GSEA #####

##find BP slim terms for light blue clade
#first get short list of GO terms
light_blue_GOIDs <- unique(up_bps_STACKED[which(up_bps_STACKED$protein_ID %in% foo[grep("darkslategray1",foo$unlist.clusterColor.),"protein_ID"]),"GO_ID"])
myCollection <- GOCollection(light_blue_GOIDs)
f1 <- system.file("extdata", "goslim_generic.obo", package="GSEABase")
slim <- getOBOCollection(f1)
View(goSlim(myCollection, slim, "BP"))


dark_mag_GOIDs <- unique(up_bps_STACKED[which(up_bps_STACKED$protein_ID %in% foo[grep("darkmagenta",foo$unlist.clusterColor.),"protein_ID"]),"GO_ID"])
myCollection_dm <- GOCollection(dark_mag_GOIDs)
View(goSlim(myCollection_dm, slim, "BP"))

#### THIS ISN"T INFORMATIVE, GO SLIM TERMS ARE TOO GENERAL ... ####

#p-value by "matrix" method
group <- 1:10
get_sim_p(sim_matrix,group=group)





####################################
####DAY 3 23 C
####################################
up_bps_spread_day3_23C <- tidyr::spread(up_bps_STACKED[grep("day3_23C", up_bps_STACKED$daytemp),], protein_ID, GO_ID)
up_bps_spread_day3_23C <- up_bps_spread_day3_23C[,-c(1:2)]

up_bps_day3_23_list <- lapply(as.list(up_bps_spread_day3_23C),function(x) x[!is.na(x)])


#run get sim matrix
sim_matrix_day3_23 <- get_sim_grid(ontology = go, information_content = GO_IC, term_sets = up_bps_day3_23_list)

#plot heatmap
heatmap3(sim_matrix_day3_23, scale = "none")
## most of the data is not similar to each other (mostly blue)
## try a cutoff and re-plot


hm <- as.dist(1-cor(t(sim_matrix_day3_23), use="pa"))
hm2 <- hclust(hm, method = 'complete')
plot(hm2)

#Identify clades in heatmap
# define some clusters

#make a list of proteins with cluster IDs based on where the tree is cut

#for complete clustering
mycl <- cutree(hm2,h=1)

#assign colors to cluster IDs
clusterCols <- colorRamps::primary.colors(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
heatmap3(sim_matrix_day3_23,cexRow = 0.5, cexCol = 0.5, RowSideColors=myClusterSideBar,RowAxisColors=1,hclustfun=hclust, distfun = function(x) as.dist(1 - cor(t(x), use = "pa")), method = "complete", scale = "none")

#get color name
clusterColor <- lapply(myClusterSideBar, color.id)
clusterColor <- lapply(clusterColor, `[[`,1)
clusterColor <- data.frame(unlist(clusterColor), stringsAsFactors = FALSE)

foo_day3_23 <- cbind(data.frame(mycl), clusterColor)
#add protein column to data frame for merging
foo_day3_23$protein_ID <- rownames(foo_day3_23)

############################################
#### Day 3 29 C
###########################################
####################################
up_bps_spread_day3_29C <- tidyr::spread(up_bps_STACKED[grep("day3_29C", up_bps_STACKED$daytemp),], protein_ID, GO_ID)
up_bps_spread_day3_29C <- up_bps_spread_day3_29C[,-c(1:2)]

up_bps_day3_29_list <- lapply(as.list(up_bps_spread_day3_29C),function(x) x[!is.na(x)])


#run get sim matrix
sim_matrix_day3_29 <- get_sim_grid(ontology = go, information_content = GO_IC, term_sets = up_bps_day3_29_list)

#plot heatmap
heatmap3(sim_matrix_day3_29, scale = "none")
## most of the data is not similar to each other (mostly blue)
## try a cutoff and re-plot


hm <- as.dist(1-cor(t(sim_matrix_day3_29), use="pa"))
hm2 <- hclust(hm, method = 'complete')
plot(hm2)

#Identify clades in heatmap
# define some clusters

#make a list of proteins with cluster IDs based on where the tree is cut

#for complete clustering
mycl <- cutree(hm2,h=1)

#assign colors to cluster IDs
clusterCols <- colorRamps::primary.colors(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
heatmap3(sim_matrix_day3_29,cexRow = 0.5, cexCol = 0.5, RowSideColors=myClusterSideBar,RowAxisColors=1,hclustfun=hclust, distfun = function(x) as.dist(1 - cor(t(x), use = "pa")), method = "complete", scale = "none")

#get color name
clusterColor <- lapply(myClusterSideBar, color.id)
clusterColor <- lapply(clusterColor, `[[`,1)
clusterColor <- data.frame(unlist(clusterColor), stringsAsFactors = FALSE)

foo_day3_29 <- cbind(data.frame(mycl), clusterColor)
#add protein column to data frame for merging
foo_day3_29$protein_ID <- rownames(foo_day3_29)



get_sim_p_from_ontology(
  ontology=go,
  information_content=GO_IC,
  term_sets=gene_GO_terms,
  group=names(beach)
)
up_bps_STACKED[which(up_bps_STACKED$protein_ID %in% foo_day3_23[which(foo_day3_23$unlist.clusterColor. == "green"),"protein_ID"]),]


library(reshape2)
melt_day0 <- melt(sim_matrix)
melt_day0 <- melt_day0[which(melt_day0$value >= 0.5),]

day0_0.5_matrix <- tidyr::spread(melt_day0, Var2, value)

rownames(day0_0.5_matrix) <- day0_0.5_matrix$Var1
day0_0.5_matrix$Var1 <- NULL




View(up_bps_STACKED[which(up_bps_STACKED$GO_ID %in% go$id),])


day0_up_bps <- unlist(up_bps[which(up_bps$daytemp=="day0_16C"),3:ncol(up_bps)])

up_uniprotIDs <- merge(up, all_giga_prots_mapped[,c(1:3,17)])
up_uniprotIDs <- up_uniprotIDs[order(up_uniprotIDs$daytemp),]
up_uniprotIDs$Entry_name <- gsub("_.*","", up_uniprotIDs$Entry_name)

write.csv(up_uniprotIDs,"~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/General_proteome_characterization/up_uniprotIDs.csv", quote = FALSE, row.names = FALSE)

hm <- as.dist(1-cor(t(NSAF_for_big_heatmap_t), use="pa"))
hm2 <- hclust(hm, method = 'complete')
plot(hm2)

#Identify clades in heatmap
# define some clusters

#make a list of proteins with cluster IDs based on where the tree is cut


#for complete clustering
mycl <- cutree(hm2,h=1.83)

#assign colors to cluster IDs
clusterCols <- colorRamps::primary.colors(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
heatmap3(NSAF_for_big_heatmap_t,cexRow = 0.5, cexCol = 0.5, RowSideColors=myClusterSideBar,ColSideColors = ColSideColors,RowAxisColors=1,ColAxisColors=1,hclustfun=hclust, distfun = function(x) as.dist(1 - cor(t(x), use = "pa")), method = "complete", scale = "row", Colv = NA)



#ASCA proteins 



heatmap3(cbind(NSAF_for_big_heatmap_t[,grep("16C|23C",colnames(NSAF_for_big_heatmap_t))],NSAF_for_big_heatmap_t[,grep("29C",colnames(NSAF_for_big_heatmap_t))]), cexRow = 0.1, cexCol = 0.8, Colv = NA)



# PCA

pca <- prcomp(log(NSAF_for_big_heatmap[,-c(1:2)],2), center = T, scale = T)
pca_meta <- cbind(NSAF_for_big_heatmap[,c(1:2)],pca$x)
colnames(pca_meta)[1:2] <- c("day","temp")

g <- ggplot(pca_meta, aes(PC1, PC2)) + geom_point(aes(col = factor(pca_meta$day), shape = factor(pca_meta$temp), size = 0.5)) + theme_bw() + ggtitle("PCA of ADJNSAF values where zeros were replaced with 0.1")
ggsave("~/Documents/GitHub/OysterSeedProject/Manuscript/MainText_DraftFigures/Fig2a.jpg",g)


# tSNE

### Run tSNE on uncorrected NSAF data
tsne_out <- Rtsne(NSAF_for_big_heatmap[,-c(1,2)],pca=FALSE,perplexity=3,theta=0.0)
ggplot(data.frame(tsne_out$Y), aes(data.frame(tsne_out$Y)[,1], data.frame(tsne_out$Y)[,2])) + geom_point(aes(col = as.factor(pca_meta$day), shape = as.factor(pca_meta$temp), size = 0.5)) + theme_bw() + ggtitle("tSNE of ADJNSAF values")

### Run tSNE on log2 transformed NSAF data
tsne_out <- Rtsne(log(NSAF_for_big_heatmap[,-c(1,2)],2),pca=FALSE,perplexity=3,theta=0.0)
ggplot(data.frame(tsne_out$Y), aes(data.frame(tsne_out$Y)[,1], data.frame(tsne_out$Y)[,2])) + geom_point(aes(col = as.factor(pca_meta$day), shape = as.factor(pca_meta$temp), size = 0.5)) + theme_bw() + ggtitle("tSNE of log2 ADJNSAF values")

### Run tSNE on scaled NSAF data
set.seed(42)
tsne_out <- Rtsne(scale(NSAF_for_big_heatmap[,-c(1,2)]),pca=FALSE,perplexity=3,theta=0.0)
ggplot(data.frame(tsne_out$Y), aes(data.frame(tsne_out$Y)[,1], data.frame(tsne_out$Y)[,2])) + geom_point(aes(col = as.factor(pca_meta$day), shape = as.factor(pca_meta$temp), size = 0.5)) + theme_bw() + ggtitle("tSNE of scaled ADJNSAF values")

tsne_out <- Rtsne(scale(log(NSAF_for_big_heatmap[,-c(1,2)],2)),pca=FALSE,perplexity=3,theta=0.0)
ggplot(data.frame(tsne_out$Y), aes(data.frame(tsne_out$Y)[,1], data.frame(tsne_out$Y)[,2])) + geom_point(aes(col = as.factor(pca_meta$day), shape = as.factor(pca_meta$temp), size = 0.5)) + theme_bw() + ggtitle("tSNE of scaled log2 ADJNSAF values")


###PLOTTING PROTEIN ABUNDANCE DISTRIBUTIONS 
NSAF_STACKED <-tidyr::gather(NSAF_for_big_heatmap[,-c(1,2)], "protein", "NSAF")
NSAF_STACKED$logNSAF <- log(NSAF_STACKED$NSAF,2)
NSAF_STACKED$scallogNSAF <- tidyr::gather(data.frame(scale(log(NSAF_for_big_heatmap[,-c(1,2)],2))))[2]
NSAF_STACKED$scallogNSAF <- NSAF_STACKED$scallogNSAF$value
NSAF_STACKED$scallogtNSAF <- tidyr::gather(data.frame(scale(log(t(NSAF_for_big_heatmap[,-c(1,2)]),2))))[2]
NSAF_STACKED$scallogtNSAF <- NSAF_STACKED$scallogtNSAF$value

#plot
ggplot(NSAF_STACKED) + geom_density(aes(NSAF , group = protein))+ ggtitle("dist. of NSAF values per protein")
ggplot(NSAF_STACKED) + geom_density(aes(logNSAF , group = protein))+ ggtitle("dist. of logNSAF values per protein")
#column mean = mean of each protein abundance from all samples
ggplot(NSAF_STACKED) + geom_density(aes(scallogNSAF, group = protein))+ ggtitle("dist. of scaled logNSAF values per protein")
#column mean = mean of all protein abundances per sample
ggplot(NSAF_STACKED) + geom_density(aes(scallogtNSAF, group = protein))+ ggtitle("dist. of scaled logtNSAF values per protein")



ggplot(tidyr::gather(log(NSAF_for_big_heatmap[,-c(1,2)],2), "protein", "NSAF")) + geom_density(aes(NSAF , group = protein)) 
ggplot(tidyr::gather(scale(log(NSAF_for_big_heatmap[,-c(1,2)],2)), "protein", "NSAF")) + geom_density(aes(NSAF , group = protein)) + ggtitle("dist. of log2 NSAF values per protein")




Create tempXdays dataframe and merge it with NSAF data
```{r}
tempXdays <- data.frame((NSAF_for_big_heatmap$day + 1) * NSAF_for_big_heatmap$temp)

NSAF <- cbind(tempXdays, NSAF_for_big_heatmap)
colnames(NSAF)[1] <- "tempXdays"
```

normalize all NSAFs by tempXdays
```{r}
NSAF_normTxD <- data.frame()

for(i in 1:nrow(NSAF)){
  row <- NSAF[i,4:ncol(NSAF)]/(NSAF$tempXdays[i])
  NSAF_normTxD <- rbind(NSAF_normTxD, row)
}
```

add temp and day info back
```{r}
NSAF_normTxD <- cbind(NSAF[,2:3], NSAF_normTxD)
```

pca <- prcomp(log(NSAF_normTxD[,-c(1:2)],2), center = T, scale = T)
pca_meta <- cbind(NSAF_normTxD[,c(1:2)],pca$x)
colnames(pca_meta)[1:2] <- c("day","temp")

ggplot(pca_meta, aes(PC1, PC2)) + geom_point(aes(col = factor(pca_meta$day), shape = factor(pca_meta$temp), size = 0.5)) + theme_bw() + ggtitle("PCA of ADJNSAF values where zeros were replaced with 0.1")










#clustering

#use bray-curtis dissimilarity for clustering

silo3_NSAF_log <- silo3_NSAF

silo3_NSAF_log[,1:6] <- log(silo3_NSAF_log[,1:6],2)

library(vegan)
nsaf.bray<-vegdist(NSAF_for_big_heatmap_t[,grep("23C", colnames(NSAF_for_big_heatmap_t))], method='bray')

nsaf.euc<-vegdist(silo3_NSAF_log[,-7], method='euclidean')

#average clustering method to cluster the data
library(cluster)
clust.avg<-hclust(nsaf.bray, method='average')
plot(clust.avg)


clust.avg<-hclust(nsaf.euc, method='average')
plot(clust.avg)


coef.hclust(clust.avg)
#coeff of ~1 means clusters are distinct and dissimilar from each other (silo3_9 = "Error in coef.hclust(clust.avg) : !is.unsorted(ht) is not TRUE")
#[1] 0.9397352
#euc:
#[1] 0.9757751

#cophenetic correlation
#how well cluster hierarchy represents original object-by-object dissimilarity space
cor(nsaf.bray, cophenetic(clust.avg))
#I think you want this to be close-ish to 1 (silo3_9 = 0.7411092)
#[1] 0.7868425
cor(nsaf.euc, cophenetic(clust.avg))
#[1] 0.7624899

#Scree plot
hclus.scree(clust.avg)

jpeg(filename = "s3_9_scree.jpeg", width = 1000, height = 1000)
hclus.scree(clust.avg)
dev.off()

#Look for the elbow/inflection point on the scree plot and you can estimate number of clusters. But  it seems that this information cannot be pulled from the scree plot. (less than 500, maybe around 300?)

#cut dendrogram at selected height (example is given for 0.5) based on what looks reasonable because SCIENCE
plot(clust.avg)
rect.hclust(clust.avg, h=0.6)

jpeg(filename = "s3_9_dendrogram.jpeg", width = 1000, height = 1000)
plot(clust.avg)
rect.hclust(clust.avg, h=0.8)
dev.off()

#this looks reasonable
clust.class<-cutree(clust.avg, h=0.55)
max(clust.class)

clust.class<-cutree(clust.avg, h=7)
max(clust.class)


#Cluster Freq table
silo3.freq <- data.frame(table(clust.class))

#Make df
silo3.clus <- data.frame(clust.class)
silo3.clus$protein_ID <- rownames(silo3.clus)
silo3_NSAF <- data.frame(NSAF_for_big_heatmap_t[,grep("23C", colnames(NSAF_for_big_heatmap_t))])
silo3_NSAF$protein_ID <- rownames(silo3_NSAF)
silo3.clus <- merge(silo3.clus, silo3_NSAF, by = "protein_ID")

silo3.clus_STACKED <- tidyr::gather(silo3.clus, "day", "NSAF", 3:8)
silo3.clus_STACKED$day <- gsub("day","",silo3.clus_STACKED$day)
silo3.clus_STACKED$day <- gsub("_23C","",silo3.clus_STACKED$day)

str(silo3.clus_STACKED)
silo3.clus_STACKED$day <- as.numeric(silo3.clus_STACKED$day)

library(ggplot2)
ggplot(silo3.clus_STACKED, aes(factor(day),NSAF)) + geom_line(aes(group = protein_ID)) + facet_wrap(~clust.class, scale = "free")


colnames(silo3_9.clus)[1] <- "S3_9.Protein"
colnames(silo3_9.clus)[2] <- "Cluster"
silo3_9.all <- merge(silo3_9.clus, silo3_9.detected, by.x = "Protein", by.y = "X")


#this gives matrix of 2 columns, first with proteins second with cluster assignment
#Line plots for each cluster
library(ggthemes)
library(reshape)
library(ggplot2)

melted_all_s3_9<-melt(silo3_9.all, id.vars=c('Protein', 'Cluster'))

ggplot(melted_all_s3_9, aes(x=variable, y=value, group=Protein)) +geom_line(alpha=0.1) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')

jpeg(filename = "silo3_9clus_lineplots.jpeg", width = 1000, height = 1000)
ggplot(melted_all_s3_9, aes(x=variable, y=value, group=Protein)) +geom_line(alpha=0.1) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')
dev.off()

#Merge Silo clusters with Silo annotated and tagged datasheet
silo3_9.annotated <- read.csv("silo3_9_annotated.csv")
silo3_9.final <- merge(silo3_9.clus, silo3_9.annotated, by.x = "Protein", by.y = "Protein")

write.csv(silo3_9.final, file = "silo3_9-anno_clus")
write.csv(silo3_9.freq, file = "silo3_9-clus_freq")
