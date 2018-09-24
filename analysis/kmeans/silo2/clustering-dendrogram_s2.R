###You'll notice that it is difficult to discern how many clusters you should choose based on the scree plot, which is to be expected from a dataset with so many variables. 
##What I did was to make a dendrogram and choose a height on the y-axis of the dendrogram that looked reasonable for separating my variables (proteins), 
##and selecting that height yielded 23 clusters.


setwd('Github/OysterSeedProject/analysis/kmeans/')
source("biostats.R")

setwd('silo2')

#Load in NSAF data
silo2 <- read.csv("silo2.csv", row.names=1)
silo2.detected <- read.csv("silo2.csv")
colnames(silo2.detected)[1] <- "X"

#use bray-curtis dissimilarity for clustering
library(vegan)
nsaf.bray<-vegdist(silo2, method='bray')

#average clustering method to cluster the data
library(cluster)
clust.avg<-hclust(nsaf.bray, method='average')
plot(clust.avg)

coef.hclust(clust.avg)
#coeff of ~1 means clusters are distinct and dissimilar from each other (silo2 = 0.9471688/0.9286963)(silo3 = 0.9403264/0.9284654) (silo9 = 0.940945/0.9250218)

#cophenetic correlation
#how well cluster hierarchy represents original object-by-object dissimilarity space
cor(nsaf.bray, cophenetic(clust.avg))
#I think you want this to be close-ish to 1 (silo2 = 0.7518792/0.742041) (silo 3 = 0.7555737/0.744824) (silo9 = 0.7613414/0.7516957)

#Scree plot
hclus.scree(clust.avg)

jpeg(filename = "s2_scree.jpeg", width = 1000, height = 1000)
hclus.scree(clust.avg)
dev.off()

#Look for the elbow/inflection point on the scree plot and you can estimate number of clusters. But  it seems that this information cannot be pulled from the scree plot. (less than 500, maybe around 300?)

#cut dendrogram at selected height (example is given for 0.5) based on what looks reasonable because SCIENCE
plot(clust.avg)
rect.hclust(clust.avg, h=0.7)

jpeg(filename = "s2_dendrogram.jpeg", width = 1000, height = 1000)
plot(clust.avg)
rect.hclust(clust.avg, h=0.7)
dev.off()

#this looks reasonable, (silo2 = 24 clusters) (silo3 = 23 clus; noexp = 25 clus) (silo9 = 16 clus)
clust.class<-cutree(clust.avg, h=0.7)
max(clust.class)

#Cluster Freq table
silo2.freq <- data.frame(table(clust.class))

#Make df
silo2.clus <- data.frame(clust.class)
names <- rownames(silo2.clus)
silo2.clus <- cbind(names, silo2.clus)
rownames(silo2.clus) <- NULL
colnames(silo2.clus)[1] <- "S2.Protein"
colnames(silo2.clus)[2] <- "Cluster"
silo2.all <- merge(silo2.clus, silo2.detected, by.x = "S2.Protein", by.y = "X")


#this gives matrix of 2 columns, first with proteins second with cluster assignment
#Line plots for each cluster
library(ggthemes)
library(reshape)
library(ggplot2)

melted_all_s2<-melt(silo2.all, id.vars=c('S2.Protein', 'Cluster'))

ggplot(melted_all_s2, aes(x=variable, y=value, group=S2.Protein)) +geom_line(alpha=0.1) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')

jpeg(filename = "silo2clus_lineplots.jpeg", width = 1000, height = 1000)
ggplot(melted_all_s2, aes(x=variable, y=value, group=S2.Protein)) +geom_line(alpha=0.1) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')
dev.off()

#Merge Silo 2 clusters with Silo 2 annotated and tagged datasheet
silo2.annotated <- read.csv("silo2_annotated.csv")
silo2.final <- merge(silo2.clus, silo2.annotated, by.x = "S2.Protein", by.y = "S2.Protein")

write.csv(silo2.final, file = "silo2-anno_clus")
write.csv(silo2.freq, file = "silo2-clus_freq")
