#load packages
library(vegan)
library(raster)
library(BioStatR)
source("biostats.R")

########Create NMDS plot with silo 3 and 9################

#upload data file
ABACUSdata <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/nb2017_sr320_ABACUS_output021417NSAF.csv", header=TRUE, stringsAsFactors = FALSE)
View(ABACUSdata)
## change column names in ABACUSdata to just sampleID
colnames(ABACUSdata) <- gsub(pattern = "X20161205_SAMPLE_", "", colnames(ABACUSdata))
colnames(ABACUSdata) <- gsub(pattern = "_ADJNSAF", "", colnames(ABACUSdata))

#upload meta data
meta_data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/Rhonda_new_sample_names.csv", header = TRUE, stringsAsFactors = FALSE)
meta_data$silo <- substr(meta_data$Contents,5,5)
meta_data$day <- substr(meta_data$SampleName,5,6)
meta_data$SampleName <- gsub(pattern = "H","",meta_data$SampleName)
meta_data$SampleName <- gsub(pattern = "C","",meta_data$SampleName)
#create a temperature column
meta_data$temp <- "temp"
for(i in 1:nrow(meta_data)){
  if(meta_data$silo[i] == "2"){
    meta_data$temp[i] <- "23"
  }
  if(meta_data$silo[i] == "3"){
    meta_data$temp[i] <- "23"
  }
  if(meta_data$silo[i] == "9"){
    meta_data$temp[i] <- "29"
  }
  if(meta_data$silo[i] == "e"){
    meta_data$temp[i] <- "16"
  }
}

#Remove NA
ABACUSdata[is.na(ABACUSdata)] <- 0

#Transpose- switch rows and columns
tABACUSdata <- t.data.frame(ABACUSdata[,-1])
colnames(tABACUSdata) <- ABACUSdata[,1]
tABACUSdata <- cbind(data.frame(rownames(tABACUSdata)),tABACUSdata)
colnames(tABACUSdata)[1] <- "SampleID"
View(tABACUSdata)

tABACUSdata <- merge(meta_data[,c(1,2,7,8)],tABACUSdata, by = "SampleID")

#Remove Silo 2 and day 15
silo3and9 <- tABACUSdata[which(substr(tABACUSdata$SampleName,1,2) != "S2" & tABACUSdata$day != "15"),]
rownames(silo3and9) <- silo3and9$SampleID
silo3and9 <- silo3and9[order(as.numeric(silo3and9$day),silo3and9$temp),]
#Make MDS dissimilarity matrix
nmds.silo3and9 <- metaMDS(silo3and9[,-c(1:4)], distance = 'euclidean', k = 2, trymax = 3000, autotransform = FALSE)

#Save as jpeg
jpeg(filename = "~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/nmdss3ands9-color.jpeg", width = 670, height = 650)
fig <- ordiplot(nmds.silo3and9, choices=c(1,2),type="none", display="sites", xlab='Axis 1', ylab='Axis 2', cex=0.5, main = "NMDS of tech reps using ABACUS_output021417NSAF.tsv from sr320/nb-2017/C_gigas/data
")
points(nmds.silo3and9, "sites", col=c(rep('black',2), rep('red',4), rep('orange',4), rep('palegreen',4), rep('springgreen4',4),
                                      rep('lightskyblue1',4), rep('blue',4)),
       pch=c(rep(19,2), rep(17,2), rep(15,2), rep(17,2), rep(15,2), rep(17,2), rep(15,2), 
             rep(17,2), rep(15,2), rep(17,2), rep(15,2), rep(17,2), rep(15,2)))
legend("topright", legend=c("pool", "23C-Silo3", "29C-Silo9"), pch=c(19,17,15))
legend("topleft", legend=c("Day 0", "Day 3", "Day 5", "Day 7", "Day 9", "Day 11", "Day 13"), 
       col=c('black', 'red', 'orange', 'palegreen','springgreen4','lightskyblue1','blue'), pch=19)
dev.off()


###########################
#compare to Sean's data Feb 14
########################

#upload data file
ABACUSdata2 <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/scaphapoda_Sean_Rhonda-2016-Oyster-Intermediates_ABACUS_output.tsv", sep = "\t",header=TRUE, stringsAsFactors = FALSE)
ABACUSdata2 <- ABACUSdata2[,c(1,grep("ADJNSAF", colnames(ABACUSdata2)))]
View(ABACUSdata2)

## change column names in ABACUSdata2 to just sampleID
colnames(ABACUSdata2) <- gsub(pattern = "X20161205_SAMPLE_", "", colnames(ABACUSdata2))
colnames(ABACUSdata2) <- gsub(pattern = "_ADJNSAF", "", colnames(ABACUSdata2))

#Remove NA
ABACUSdata2[is.na(ABACUSdata2)] <- 0

#Transpose- switch rows and columns
tABACUSdata2 <- t.data.frame(ABACUSdata2[,-1])
colnames(tABACUSdata2) <- ABACUSdata2[,1]
tABACUSdata2 <- cbind(data.frame(rownames(tABACUSdata2)),tABACUSdata2)
colnames(tABACUSdata2)[1] <- "SampleID"
View(tABACUSdata2)

tABACUSdata2 <- merge(meta_data[,c(1,2,7,8)],tABACUSdata2, by = "SampleID")

#Remove Silo 2 and day 15
silo3and9v2 <- tABACUSdata2[which(substr(tABACUSdata2$SampleName,1,2) != "S2" & tABACUSdata2$day != "15"),]
rownames(silo3and9v2) <- silo3and9v2$SampleID
silo3and9v2 <- silo3and9v2[order(as.numeric(silo3and9v2$day),silo3and9v2$temp),]

#Make MDS dissimilarity matrix
nmds.silo3and9v2 <- metaMDS(silo3and9v2[,-c(1:4)], distance = 'euclidean', k = 2, trymax = 3000, autotransform = FALSE)

#Save as jpeg
jpeg(filename = "~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/nmdss3ands9_SEAN-color.jpeg", width = 670, height = 650)
fig <- ordiplot(nmds.silo3and9v2, choices=c(1,2),type="none", display="sites", xlab='Axis 1', ylab='Axis 2', cex=0.5, main = "NMDS of tech reps ADJNSAF vals from \nhttp://owl.fish.washington.edu/\nscaphapoda/Sean/Rhonda-2016-Oyster-Intermediates/ABACUS_output.tsv
")
points(nmds.silo3and9v2, "sites", col=c(rep('black',2), rep('red',4), rep('orange',4), rep('palegreen',4), rep('springgreen4',4),
                                      rep('lightskyblue1',4), rep('blue',4)),
       pch=c(rep(19,2), rep(17,2), rep(15,2), rep(17,2), rep(15,2), rep(17,2), rep(15,2), 
             rep(17,2), rep(15,2), rep(17,2), rep(15,2), rep(17,2), rep(15,2)))
legend("topright", legend=c("pool", "23C-Silo3", "29C-Silo9"), pch=c(19,17,15))
legend("topleft", legend=c("Day 0", "Day 3", "Day 5", "Day 7", "Day 9", "Day 11", "Day 13"), 
       col=c('black', 'red', 'orange', 'palegreen','springgreen4','lightskyblue1','blue'), pch=19)
dev.off()



###################################
#Sean's data from Mar 1; this is the same as sean's other file
###################################

#upload data file
ABACUSdata3 <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output.tsv", sep = "\t",header=TRUE, stringsAsFactors = FALSE)
ABACUSdata3 <- ABACUSdata3[,c(1,grep("ADJNSAF", colnames(ABACUSdata3)))]
View(ABACUSdata3)

## change column names in ABACUSdata3 to just sampleID
colnames(ABACUSdata3) <- gsub(pattern = "X20161205_SAMPLE_", "", colnames(ABACUSdata3))
colnames(ABACUSdata3) <- gsub(pattern = "_ADJNSAF", "", colnames(ABACUSdata3))

#Remove NA
ABACUSdata3[is.na(ABACUSdata3)] <- 0

#Transpose- switch rows and columns
tABACUSdata3 <- t.data.frame(ABACUSdata3[,-1])
colnames(tABACUSdata3) <- ABACUSdata3[,1]
tABACUSdata3 <- cbind(data.frame(rownames(tABACUSdata3)),tABACUSdata3)
colnames(tABACUSdata3)[1] <- "SampleID"
View(tABACUSdata3)

tABACUSdata3 <- merge(meta_data[,c(1,2,7,8)],tABACUSdata3, by = "SampleID")

#Remove Silo 2 and day 15
silo3and9v3 <- tABACUSdata3[which(substr(tABACUSdata3$SampleName,1,2) != "S2" & tABACUSdata3$day != "15"),]
rownames(silo3and9v3) <- silo3and9v3$SampleID
silo3and9v3 <- silo3and9v3[order(as.numeric(silo3and9v3$day),silo3and9v3$temp),]

#Make MDS dissimilarity matrix
nmds.silo3and9v3 <- metaMDS(silo3and9v3[,-c(1:4)], distance = 'euclidean', k = 2, trymax = 3000, autotransform = FALSE)

#Save as jpeg
jpeg(filename = "~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/nmdss3ands9_seanMAR1-color.jpeg", width = 670, height = 650)
fig <- ordiplot(nmds.silo3and9v3, choices=c(1,2),type="none", display="sites", xlab='Axis 1', ylab='Axis 2', cex=0.5, main = "NMDS of tech reps ADJNSAF vals from \nhttp://owl.fish.washington.edu/\nscaphapoda/Sean/Rhonda-2016-Oyster-Intermediates/ABACUS_output.tsv
                ")
points(nmds.silo3and9v3, "sites", col=c(rep('black',2), rep('red',4), rep('orange',4), rep('palegreen',4), rep('springgreen4',4),
                                        rep('lightskyblue1',4), rep('blue',4)),
       pch=c(rep(19,2), rep(17,2), rep(15,2), rep(17,2), rep(15,2), rep(17,2), rep(15,2), 
             rep(17,2), rep(15,2), rep(17,2), rep(15,2), rep(17,2), rep(15,2)))
legend("topright", legend=c("pool", "23C-Silo3", "29C-Silo9"), pch=c(19,17,15))
legend("topleft", legend=c("Day 0", "Day 3", "Day 5", "Day 7", "Day 9", "Day 11", "Day 13"), 
       col=c('black', 'red', 'orange', 'palegreen','springgreen4','lightskyblue1','blue'), pch=19)
dev.off()


###Stevens data from 02132017
#upload data file
ABACUSdata4 <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output_sr320nb02132017.tsv", sep = "\t" , header=TRUE, stringsAsFactors = FALSE)
ABACUSdata4 <- ABACUSdata4[,c(1,grep("ADJNSAF", colnames(ABACUSdata4)))]

## change column names in ABACUSdata4 to just sampleID
colnames(ABACUSdata4) <- gsub(pattern = "X20161205_SAMPLE_", "", colnames(ABACUSdata4))
colnames(ABACUSdata4) <- gsub(pattern = "_ADJNSAF", "", colnames(ABACUSdata4))
View(ABACUSdata4)

###Stevens data from 03042017

#upload data file
ABACUSdata5 <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/headABACUS_output_sr320nb03042017.tsv", sep = "\t" , header=TRUE, stringsAsFactors = FALSE)
ABACUSdata5 <- ABACUSdata5[,c(1,grep("ADJNSAF", colnames(ABACUSdata5)))]

## change column names in ABACUSdata4 to just sampleID
colnames(ABACUSdata4) <- gsub(pattern = "X20161205_SAMPLE_", "", colnames(ABACUSdata4))
colnames(ABACUSdata4) <- gsub(pattern = "_ADJNSAF", "", colnames(ABACUSdata4))
View(ABACUSdata4)


