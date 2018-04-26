#K MEANS CLUSTERING

#install packages
# install.packages("cluster", "pvclust", "fpc", "BioStatR", "vegan")
library(vegan, cluster, pvclust, fpc, Biostats)
library(vegan)
library(cluster)
library(fpc)
library(BioStatR)
library(pvclust)
setwd("Documents/kmeans/")
source("biostats.R")

#load silo in
silo2 <- read.csv("silo2.csv", row.names = 1)

#make a euclidean distance matrix
silo2.eucd <- vegdist(silo2, method='euclidean')

#look at the scree plot to determine the number of clusters you want to choose. Generally, choose the highest scree value >2

#PCA and Scree plot
pca_silo9 <- princomp(silo9)
summary(pca_silo9)
screeplot(pca_silo9)

jpeg(filename = "scree_s9.jpeg", width = 1000, height = 1000)
screeplot(pca_silo9)
dev.off()

#scree code
nhclus.scree(silo2.eucd, max.k=7047)

#remove undetected proteins, look at elbow, run kmeans
nsafcl.s2kmeans <- kmeans(nsaf, tra, centers=3, iter.max=1000, nstart=25)

setwd('Documents/robertslab/labnotebook/analysis/kmeans/')
getwd()
jpeg(filename = "nhclus-scree_s2.jpeg", width = 1000, height = 1000)
nsafcl.kmeans <- kmeans(nsaf, tra, centers=3, iter.max=1000, nstart=25)
dev.off()

#let's say you choose 4 clusters based on the scree plot
silo9.4clus <- kmeans(silo9.eucd, centers=4, iter.max=10000, nstart=25)
silo9.4clus


#this will give you an output with the cluter assignments for each protein

#plot things
jpeg(filename = "kmeans_s2.jpeg", width = 1000, height = 1000)
plot(silo2[c("X2_15", "X2_3")], col = silo2.4clus$cluster)
dev.off()


df$newcolumnname
