###You'll notice that it is difficult to discern how many clusters you should choose based on the scree plot, which is to be expected from a dataset with so many variables. 
##What I did was to make a dendrogram and choose a height on the y-axis of the dendrogram that looked reasonable for separating my variables (proteins), 
##and selecting that height yielded 23 clusters.

#Load in NSAF data
silo3 <- read.csv("~/Documents/robertslab/labnotebook/raw_data/silo3.csv", row.names=1)
silo3.detected <- read.csv("~/Documents/robertslab/labnotebook/raw_data/silo3.csv")
colnames(silo3.detected)[1] <- "X"

#use bray-curtis dissimilarity for clustering
library(vegan)
nsaf.bray<-vegdist(silo3, method='bray')

#average clustering method to cluster the data
library(cluster)
clust.avg<-hclust(nsaf.bray, method='average')
plot(clust.avg)

coef.hclust(clust.avg)
#coeff of ~1 means clusters are distinct and dissimilar from each other (silo2 = 0.9471688)(silo3 = 0.9403264) (silo9 = 0.940945)

#cophenetic correlation
#how well cluster hierarchy represents original object-by-object dissimilarity space
cor(nsaf.bray, cophenetic(clust.avg))
#I think you want this to be close-ish to 1 (silo2 = 0.7518792) (silo 3 = 0.7555737) (silo9 = 0.7613414)

#Scree plot
source("biostats.R")
hclus.scree(clust.avg)

jpeg(filename = "s3_scree.jpeg", width = 1000, height = 1000)
hclus.scree(clust.avg)
dev.off()

#Look for the elbow/inflection point on the scree plot and you can estimate number of clusters. But  it seems that this information cannot be pulled from the scree plot. (less than 500, maybe around 300?)

#cut dendrogram at selected height (example is given for 0.5) based on what looks reasonable because SCIENCE
plot(clust.avg)
rect.hclust(clust.avg, h=0.6)

jpeg(filename = "s3_dendrogram.jpeg", width = 1000, height = 1000)
plot(clust.avg)
rect.hclust(clust.avg, h=0.6)
dev.off()

#this looks reasonable, (silo2 = 24 clusters) (silo3 = 23 clus) (silo9 = 16 clus)
clust.class<-cutree(clust.avg, h=0.6)
max(clust.class)

#Cluster Freq table
silo3.freq <- data.frame(table(clust.class))

#Make df
silo3.clus <- data.frame(clust.class)
names <- rownames(silo3.clus)
silo3.clus <- cbind(names, silo3.clus)
rownames(silo3.clus) <- NULL
colnames(silo3.clus)[1] <- "Protein"
colnames(silo3.clus)[2] <- "Cluster"
silo3.clus <- merge(silo3.clus, silo3.detected, by.x = "Protein", by.y = "X")


#this gives matrix of 2 columns, first with proteins second with cluster assignment
#Line plots for each cluster
library(ggthemes)
library(reshape)
library(ggplot2)

melted_all_s3<-melt(silo3.clus, id.vars=c('Protein', 'Cluster'))

ggplot(melted_all_s3, aes(x=variable, y=value, group=Protein)) +geom_line(alpha=0.1) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')

jpeg(filename = "silo3clus_lineplots.jpeg", width = 1000, height = 1000)
ggplot(melted_all, aes(x=variable, y=value, group=Protein)) +geom_line(alpha=0.1) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')
dev.off()


