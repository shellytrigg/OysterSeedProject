###You'll notice that it is difficult to discern how many clusters you should choose based on the scree plot, which is to be expected from a dataset with so many variables. 
##What I did was to make a dendrogram and choose a height on the y-axis of the dendrogram that looked reasonable for separating my variables (proteins), 
##and selecting that height yielded 23 clusters.


source("biostats.R")

#Combine silo 3 and silo 9 data
setwd("Documents/Kaitlyn/Github/OysterSeedProject/analysis/kmeans/silo9/")
system("cat ../silo3/silo3.csv >> ../Silo3_and_9/silo3_9.csv")
system("cat silo9.csv >> ../Silo3_and_9/silo3_9.csv")
#restart R

#Load in NSAF data
setwd("Documents/Kaitlyn/Github/OysterSeedProject/analysis/kmeans/Silo3_and_9/")
silo3_9 <- read.csv("silo3and9.csv", row.names=1)
silo3_9.detected <- read.csv("silo3and9.csv")
colnames(silo3_9.detected)[1] <- "X"

#use bray-curtis dissimilarity for clustering
library(vegan)
nsaf.bray<-vegdist(silo3_9, method='bray')

#average clustering method to cluster the data
library(cluster)
clust.avg<-hclust(nsaf.bray, method='average')
plot(clust.avg)


coef.hclust(clust.avg)
#coeff of ~1 means clusters are distinct and dissimilar from each other (silo3_9 = "Error in coef.hclust(clust.avg) : !is.unsorted(ht) is not TRUE")

#cophenetic correlation
#how well cluster hierarchy represents original object-by-object dissimilarity space
cor(nsaf.bray, cophenetic(clust.avg))
#I think you want this to be close-ish to 1 (silo3_9 = 0.7411092)

#Scree plot
hclus.scree(clust.avg)

jpeg(filename = "s3_9_scree.jpeg", width = 1000, height = 1000)
hclus.scree(clust.avg)
dev.off()

#Look for the elbow/inflection point on the scree plot and you can estimate number of clusters. But  it seems that this information cannot be pulled from the scree plot. (less than 500, maybe around 300?)

#cut dendrogram at selected height (example is given for 0.5) based on what looks reasonable because SCIENCE
plot(clust.avg)
rect.hclust(clust.avg, h=0.8)

jpeg(filename = "s3_9_dendrogram.jpeg", width = 1000, height = 1000)
plot(clust.avg)
rect.hclust(clust.avg, h=0.8)
dev.off()

#this looks reasonable
clust.class<-cutree(clust.avg, h=0.8)
max(clust.class)

#Cluster Freq table
silo3_9.freq <- data.frame(table(clust.class))

#Make df
silo3_9.clus <- data.frame(clust.class)
names <- rownames(silo3_9.clus)
silo3_9.clus <- cbind(names, silo3_9.clus)
rownames(silo3_9.clus) <- NULL
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
