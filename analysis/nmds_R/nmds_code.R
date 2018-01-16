#Update R (do on Rgui)
install.packages("installr"); 
library(installr)
updateR()

#Set WD
getwd()
setwd("C:/Users/Katty/Documents/robertslab/oysterseedproject/R_work/")

#Install and load packages
install.packages("vegan")
library(vegan)
install.packages("raster")
library(raster)
install.packages("BioStatR")
library(BioStatR)
source("biostats.R")

########Create NMDS plot with silo 3 and 9################

#upload file
ABACUSdata <- read.csv("~/Documents/robertslab/oysterseedproject/raw_os_data/ABACUSdata.csv", header=TRUE)
View(ABACUSdata)

#Remove NA
ABACUSdata[is.na(ABACUSdata)] <- 0

#Transpose- switch rows and columns
tABACUSdata <- t(ABACUSdata)
View(tABACUSdata)

#Rename Columns and remove row
colnames(tABACUSdata) <- tABACUSdata[1,]
tABACUSdata = tABACUSdata[-1,]
Head(tABACUSdata)

#Convert to numeric
is.numeric(tABACUSdata)

#Remove Silo 2
silo3and9 <- tABACUSdata[-(seq(from = 2, to = 22, by = 3)), ]

#Make MDS dissimilarity matrix
nmds.silo3and9 <- metaMDS(silo3and9, distance = 'euclidean', k = 2, trymax = 3000, autotransform = FALSE)

#Save as jpeg
jpeg(filename = "nmdss3ands9.jpeg", width = 1000, height = 1000)
ordiplot(nmds.silo3and9, choices=c(1,2), type="text", display="sites")
dev.off()



