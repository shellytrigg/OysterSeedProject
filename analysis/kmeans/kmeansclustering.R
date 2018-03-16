#K MEANS CLUSTERING

#install packages for vegan and cluster
install.packages("vegan")
library(permute, lattice)
library(vegan)
install.packages("cluster")
library(cluster)

#read in file with things you want to cluster (i.e., protein names) as rows and variables you want to cluster across (i.e., NSAF for each day & temperature) as columns
allsilos <- read.csv('Documents/robertslab/labnotebook/raw_data/ABACUSdata.csv', header=T, row.names=1 )

#if all data are on the same scle, you do not need to standardize (https://stackoverflow.com/questions/14438942/standardization-of-data-in-r)

#Make individual silos
silo2 <- allsilos[, (seq(from = 2, to = 22, by = 3))]
silo3 <- allsilos[, (seq(from = 3, to = 22, by = 3))]
silo9 <- allsilos[, (seq(from = 4, to = 22, by = 3))]

#make a euclidean distance matrix
silo9.eucd <- vegdist(silo9, method='euclidean')

#look at the scree plot to determine the number of clusters you want to choose. Generally, choose the highest scree value >2

#PCA and Scree plot
pca_silo9 <- princomp(silo9)
summary(pca_silo9)
screeplot(pca_silo9)

#set wd
setwd('Documents/robertslab/labnotebook/analysis/kmeans/')
getwd()

#Save scree plot
jpeg(filename = "scree_s9.jpeg", width = 1000, height = 1000)
screeplot(pca_silo9)
dev.off()

#let's say you choose 4 clusters based on the scree plot

silo9.4clus <- kmeans(silo9.eucd, centers=4, iter.max=10000, nstart=25)
silo9.4clus


#this will give you an output with the cluter assignments for each protein

#plot things
jpeg(filename = "kmeans_s2.jpeg", width = 1000, height = 1000)
plot(silo2[c("X2_15", "X2_3")], col = silo2.4clus$cluster)
dev.off()
