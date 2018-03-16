#K MEANS CLUSTERING

#install packages for vegan and cluster
install.packages("vegan")
library(permute)
library(lattice)
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
silo2.eucd <- vegdist(silo2, method='euclidean')

      #look at the scree plot to determine the number of clusters you want to choose. Generally, choose the highest scree value >2

#PCA and Scree plot
pca_silo2 <- princomp(silo2)
summary(pca_silo2)
screeplot(pca_silo2)

#set wd
setwd('Documents/robertslab/labnotebook/analysis/kmeans/')
getwd()

#Save scree plot
jpeg(filename = "scree_s2.jpeg", width = 1000, height = 1000)
screeplot(pca_silo2)
dev.off()

      #let's say you choose 4 clusters based on the scree plot

silo2.3clus <- kmeans(silo2.eucd, centers=3, iter.max=10000, nstart=25)
silo2.2clus


#this will give you an output with the cluter assignments for each protein

#-------------------------

#EIGENVECTORS
#this assumes an NMDS object (called "nmds") created in R with the vegan package and a data file, mine are usually log(x+1) transformed (called "dat.trans")

dat.eigen<-envfit(nmds$points, dat.trans, perm=1000)
dat.eigen
#this will spit out a few columns of numbers in your R console. You can copy and paste to excel and then play with them there. I haven't found a good way of exporting or manipulating in R. One column will have protein names, another with have eigenvector weights on NMDS1 (value of 0-1), weights on NMDS2, and a p-value.