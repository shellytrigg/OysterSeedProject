###You'll notice that it is difficult to discern how many clusters you should choose based on the scree plot, which is to be expected from a dataset with so many variables. 
##What I did was to make a dendrogram and choose a height on the y-axis of the dendrogram that looked reasonable for separating my variables (proteins), 
##and selecting that height yielded 23 clusters.

setwd('Github/OysterSeedProject/analysis/kmeans/')
source("biostats.R")

setwd('silo9')


#Load in NSAF data
silo9 <- read.csv("silo9.csv", row.names=1)
silo9.detected <- read.csv("silo9.csv")
colnames(silo9.detected)[1] <- "X"

#use bray-curtis dissimilarity for clustering
library(vegan)
nsaf.bray<-vegdist(silo9, method='bray')

#average clustering method to cluster the data
library(cluster)
clust.avg<-hclust(nsaf.bray, method='average')
plot(clust.avg)

coef.hclust(clust.avg)
#coeff of ~1 means clusters are distinct and dissimilar from each other (silo2 = 0.9471688)(silo3 = 0.9403264/0.9284654) (silo9 = 0.940945/0.9250218)

#cophenetic correlation
#how well cluster hierarchy represents original object-by-object dissimilarity space
cor(nsaf.bray, cophenetic(clust.avg))
#I think you want this to be close-ish to 1 (silo2 = 0.7518792) (silo 3 = 0.7555737/0.744824) (silo9 = 0.7613414/0.7516957)

#Scree plot
hclus.scree(clust.avg)

jpeg(filename = "s9_scree.jpeg", width = 1000, height = 1000)
hclus.scree(clust.avg)
dev.off()

#Look for the elbow/inflection point on the scree plot and you can estimate number of clusters. But  it seems that this information cannot be pulled from the scree plot. (less than 500, maybe around 300?)

#cut dendrogram at selected height (example is given for 0.5) based on what looks reasonable because SCIENCE
plot(clust.avg)
rect.hclust(clust.avg, h=0.7)

jpeg(filename = "s9_dendrogram.jpeg", width = 1000, height = 1000)
plot(clust.avg)
rect.hclust(clust.avg, h=0.7)
dev.off()

#this looks reasonable
clust.class<-cutree(clust.avg, h=0.7)
max(clust.class)

#Cluster Freq table
silo9.freq <- data.frame(table(clust.class))

#Make df
silo9.clus <- data.frame(clust.class)
names <- rownames(silo9.clus)
silo9.clus <- cbind(names, silo9.clus)
rownames(silo9.clus) <- NULL
colnames(silo9.clus)[1] <- "S9.Protein"
colnames(silo9.clus)[2] <- "Cluster"
silo9.all <- merge(silo9.clus, silo9.detected, by.x = "S9.Protein", by.y = "X")


#this gives matrix of 2 columns, first with proteins second with cluster assignment
#Line plots for each cluster
library(ggthemes)
library(reshape)
library(ggplot2)

melted_all_s9<-melt(silo9.all, id.vars=c('S9.Protein', 'Cluster'))

ggplot(melted_all_s9, aes(x=variable, y=value, group=S9.Protein)) +geom_line(alpha=0.1) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')

jpeg(filename = "silo9clus_lineplots.jpeg", width = 1000, height = 1000)
ggplot(melted_all_s9, aes(x=variable, y=value, group=S9.Protein)) +geom_line(alpha=0.1) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')
dev.off()

#Merge Silo 9 clusters with Silo 9 annotated and tagged datasheet
silo9.annotated <- read.csv("silo9_annotated.csv")
silo9.final <- merge(silo9.clus, silo9.annotated, by.x = "S9.Protein", by.y = "S9.Protein")

write.csv(silo9.final, file = "silo9-anno_clus")
write.csv(silo9.freq, file = "silo9-clus_freq")
