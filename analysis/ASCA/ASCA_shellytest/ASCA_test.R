library(dplyr)
library(tidyr)
library(MetStaT)
library(MetaboAnalystR)
library(ggplot2)
library(heatmap3)

#Messing around with ASCA example:
## use the data matrix, 'ASCAX', and an experimental design matrix, 'ASCAF'
data(ASCAdata)
ASCA <- ASCA.Calculate(ASCAX, ASCAF, equation.elements = "1,2,12", scaling = FALSE)

#Check format of example data file
View(ASCAX)
#There are 60 rows and 4 columns

#Check format of example levels file
View(ASCAF)
#there are 60 rows and two factors; (i.e. first column could be temperature
#and second column could be timepoint)

## plotting the results
## plot the loadings of the first two principal components of the first factor (i.e. temp.)
ASCA.PlotLoadings(ASCA, ee = "1", pcs="1,2")
#plot single score plot for the second factor on the first two PCs
ASCA.PlotScores(ASCA, ee = "2", PCs = "1,2")
## plot the scores for the first two principal components and the projections of 
## the data for the second factor (i.e. time)
ASCA.PlotScoresPerLevel(ASCA, ee = "2", pcs = "1,2")
## the data for the first factor (i.e. temp.)
ASCA.PlotScoresPerLevel(ASCA, ee = "1", pcs = "1,2")
## the data for the interaction of factors (i.e. temp. and time interactive effect)
ASCA.PlotScoresPerLevel(ASCA, ee = "12", pcs = "1,2")

#print ASCA summary showing variance explained per component
ASCA.GetSummary(ASCA)

## Do a permutation test to evaluate the significance to the two factors and the interaction.
ASCA.DoPermutationTest(ASCA, perm=1000)
#output
#1     2    12
#0.001 0.003 0.191
#means that factor 1 and factor 2 both have significant effects on molecules, but the interaction of the factors does not have a significant effect



#read in Oyster Temp. protein data frame
data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/kmeans/Silo3_and_9/silo3and9.csv")
#create silo column
data$silo <- substr(data$Protein,1,1)
#change class of protein column from factor to character
data$Protein <- as.character(unlist(data$Protein))
#remove the silo identifier from the protein names
data$Protein <- substr(data$Protein, 3, nchar(data$Protein))
#Change table orientation so that all abundances are listed in one column 
#and time points are listed down a column
#In this command, we only want to apply it to the timepoint columns that 
#have abundance values and keep the protein and silo columns unchanged
#so we use c(2:9) to select for timepoint columns with abundance values only
data <- tidyr::gather(data, "Time", "abundance", c(2:9))
#Change the table orientation again so that proteins are listed across the top
#This command will condense the proteins so that only unique proteins are listed
#across the top
data <- tidyr::spread(data, "Protein", "abundance")
#change silo column class to numeric
data$silo <- as.numeric(data$silo)
#change time column class to numeric; exclude the X with the substr command
data$Time <- as.numeric(substr(data$Time,2,nchar(data$Time)))

#format and export table for metaboanalyst
data$sample <- paste("S",data$silo,"T",data$Time, sep = "")
#move sample columns to first column
data <- data[,c(ncol(data),1:ncol(data)-1)]
#replace NAs with non-zero value so we can run ASCA
data[is.na(data)] <- 0.1

write.csv(data, "~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_shellytest/silo3_9_reformat4MetabA.csv", row.names = FALSE, quote = FALSE)

#find max abundance value
max(data[,4:ncol(data)], na.rm = TRUE)
# 379

#find min abundance value
min(data[data > 0], na.rm = TRUE)
#0.5


#create matrix to pass to ASCA command, excluding the silo and time info
ASCA_X <- as.matrix(data[,-c(1:3)])
#create matrix to pass to ASCA command with only the silo and time info
ASCA_F <- as.matrix(data[,2:3])


#Run ASCA command
ASCA <- ASCA.Calculate(ASCA_X, ASCA_F, equation.elements = "1,2,12", scaling = FALSE)

#make data frames of all arrays in ASCA containing protein names.  This way we can compare them and see if they are all in the same order.
#If they are, we can assume the data in svd$v matches the order of the column names
d1 <- data.frame(colnames(ASCA$data))
d2 <- data.frame(colnames(ASCA$remainder))
d3 <- data.frame(colnames(ASCA$`1`$reduced.matrix))
d4 <- data.frame(colnames(ASCA$`2`$reduced.matrix))
d5 <- data.frame(colnames(ASCA$`12`$reduced.matrix))
#binding all data frames of protein names together
test <- cbind(d1,d2,d3,d4,d5)
colnames(test) <- c("d1","d2", "d3","d4","d5")
View(test)
#these all appear the same, no order is changed 
#test if columns are identical
identical(test$d1,test$d2)
identical(test$d1,test$d3)
identical(test$d1,test$d4)
identical(test$d1,test$d5)
identical(test$d2,test$d3)
identical(test$d2,test$d4)
identical(test$d2,test$d5)
identical(test$d3,test$d4)
identical(test$d3,test$d5)
identical(test$d4,test$d5)
#all are true; protein orders are identical

#create a dataframe with protein names and loadings
Temp_loadings <- cbind(d1,ASCA$`1`$svd$v[,1])
#rename columns
colnames(Temp_loadings) <- c("protein", "PC1_loadings")
#output loadings table
write.csv(Temp_loadings, "~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_shellytest/ACSAr_temp_loadings.csv", row.names = FALSE, quote = FALSE)

#make a short list of proteins with cut-offs based on plotting loadings curve
Temp_loadings_cut <- data.frame(Temp_loadings[which(Temp_loadings$PC1_loadings >= 0.05 | Temp_loadings$PC1_loadings <= -0.05),])
#count the number of proteins remaining after cut-off
nrow(Temp_loadings_cut)
#48

#select proteins that made the loadings cutoff out of the data
data_PC1_0.05_selects <- data[,c(1:3,which(colnames(data) %in% Temp_loadings_cut$protein))]
#count the number of proteins to verify selection worked; first 3 columns are labels
ncol(data_PC1_0.05_selects[,4:ncol(data_PC1_0.05_selects)])
#48

#output table for uploading into metaboanalyst and seeing what is going on with these proteins
write.csv(data_PC1_0.05_selects, "~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_shellytest/data_PC1_0.05_selects.csv", row.names = FALSE, quote = FALSE)


###########EDIT DATA FOR METABOANALYST ANALYSIS##############
#read in Oyster Temp. protein data frame
data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/kmeans/Silo3_and_9/silo3and9.csv")
#create silo column
data$silo <- substr(data$Protein,1,1)
#change class of protein column from factor to character
data$Protein <- as.character(unlist(data$Protein))
#remove the silo identifier from the protein names
data$Protein <- substr(data$Protein, 3, nchar(data$Protein))
#Change table orientation so that all abundances are listed in one column 
#and time points are listed down a column
#In this command, we only want to apply it to the timepoint columns that 
#have abundance values and keep the protein and silo columns unchanged
#so we use c(2:9) to select for timepoint columns with abundance values only
data <- tidyr::gather(data, "Time", "abundance", c(2:9))
#Change the table orientation again so that proteins are listed across the top
#This command will condense the proteins so that only unique proteins are listed
#across the top
data <- tidyr::spread(data, "Protein", "abundance")
#change silo column class to numeric
data$silo <- as.numeric(data$silo)
#change time column class to numeric; exclude the X with the substr command
data$Time <- as.numeric(substr(data$Time,2,nchar(data$Time)))

#format and export table for metaboanalyst
data$sample <- paste("S",data$silo,"T",data$Time, sep = "")
#move sample columns to first column
data <- data[,c(ncol(data),1:ncol(data)-1)]
data_choyp <- data[,c(1:3,grep("CHOYP", colnames(data)))]
data_std <- data[,-grep("CHOYP", colnames(data))]
write.csv(data_choyp, "~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_shellytest/silo3_9_CHOYP_reformat4MetabA.csv", row.names = FALSE, quote = FALSE)
write.csv(data_std, "~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_shellytest/silo3_9_STD_reformat4MetabA.csv", row.names = FALSE, quote = FALSE)

#running metaboanalyst in R
#read in timeseries data
mSet<-InitDataObjects("conc", "ts", FALSE)
mSet<-SetDesignType(mSet, "time")
mSet<-Read.TextData(mSet, "~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_shellytest/silo3_9_CHOYP_reformat4MetabA.csv", "rowts", "disc")
mSet<-SanityCheckData(mSet)
#exclude any peptides with missing or zero values
mSet<-RemoveMissingPercent(mSet, percent=0.5)
mSet<-ImputeVar(mSet, method="exclude")
#filter by interquartile range to exclude proteins with near constant values throughout the experiment: 
#mSet<-FilterVariable(mSet, "iqr", "F", 25)

#do mean-centering
mSet<-Normalization(mSet, "NULL", "NULL", "MeanCenter", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_CHOYP_nozeros_", "png", 72, width=NA)
#make a data frame out of mean-centered peptides with zero-value peptides excluded
normData <- data.frame(mSet$dataSet$norm)

#create matrix to pass to ASCA command, excluding the silo and time info
ASCA_normData_X <- as.matrix(normData)
#create matrix to pass to ASCA command with only the silo and time info
ASCA_normData_F <- cbind(as.matrix(as.numeric(substr(rownames(normData),2,2))), as.matrix(as.numeric(substr(rownames(normData),4,length(rownames(normData))))))

#Run ASCA command
ASCA <- ASCA.Calculate(ASCA_normData_X, ASCA_normData_F, equation.elements = "1,2,12", scaling = FALSE)
#plot loadings for PC1 and PC2 for temperature, time, and their interaction
ASCA.PlotLoadings(ASCA, ee = "1", pcs="1,2")
ASCA.PlotLoadings(ASCA, ee = "2", pcs="1,2")
ASCA.PlotLoadings(ASCA, ee = "12", pcs="1,2")
ASCA.PlotScores(ASCA, ee = "2", PCs = "1,2")
## plot the scores for the first two principal components and the projections of 
ASCA.PlotScoresPerLevel(ASCA, ee = "1", pcs = "1,2")
ASCA.PlotScoresPerLevel(ASCA, ee = "2", pcs = "1,2")
ASCA.PlotScoresPerLevel(ASCA, ee = "12", pcs = "1,2")

ASCA.GetSummary(ASCA)
ASCA.DoPermutationTest(ASCA, perm=1000)
# 1     2    12
# 0.139 0.009 0.995
#means that time has a significant effect on molecules, but temperature or time x temp interaction of the factors does not have a significant effect

plot(ASCA$`1`$svd$v[,1],ASCA$`1`$svd$v[,2])

protnames <- data.frame(colnames(ASCA$data))
d <- cbind(protnames, ASCA$`1`$svd$v[,1])
colnames(d) <- c("protein", "PC1loadings")
d_great <- d[which(d$PC1loadings > 0),]
d_less <- d[which(d$PC1loadings < 0),]

#These next plots are to observe a point with diminishing returns
#print plot PC1 loadings values for proteins with loadings values > 0
plot(d_great[order(d_great$PC1loadings, decreasing = TRUE),2])
#print plot PC1 loadings values for proteins with loadings values < 0
plot(d_less[order(d_less$PC1loadings, decreasing = TRUE),2])

#make list of proteins with temperature PC1 loadings values >= 0.04
cutd <- d[which(abs(d$PC1loadings) >= 0.04),]
View(normData)
cutnormD <-normData[,which(colnames(normData) %in% cutd$protein)]
#first remove pipe
cutDnonNorm <- data_choyp
colnames(cutDnonNorm) <- gsub("\\|","",colnames(cutDnonNorm))
cutDnonNorm <- cutDnonNorm[,c(2:3,which(colnames(cutDnonNorm) %in% cutd$protein))]


cutnormD$silo <- as.numeric(substr(rownames(cutnormD),2,2))
cutnormD$time <- as.matrix(as.numeric(substr(rownames(cutnormD),4,length(rownames(cutnormD)))))
#reorder columns
cutnormD <- cutnormD[,c(60,61,1:59)]
#order table by silo and time
ord_cutnormD <- cutnormD[order(cutnormD$silo,cutnormD$time),]
#make a heatmap to see if proteins cluster based on abundance over time
heatmap3(as.matrix(t(ord_cutnormD[-c(1,9),3:61])), Colv = NA)

ord_cutnormD_long <- tidyr::gather(ord_cutnormD, "protein","abundance", 3:ncol(ord_cutnormD))
class(ord_cutnormD_long$silo) <- "character"
ord_cutnormD_long <- merge(ord_cutnormD_long, d)
#make line plots of abundance over time for individual proteins
ggplot(ord_cutnormD_long, aes(x= time, y = abundance)) + geom_line(aes(colour = silo)) + facet_wrap(~protein, ncol = 10)
ggplot(ord_cutnormD_long[which(ord_cutnormD_long$PC1loadings > 0.1),], aes(x= time, y = abundance)) + geom_line(aes(colour = silo)) + facet_wrap(~protein)
ggplot(ord_cutnormD_long[which(ord_cutnormD_long$PC1loadings < -0.05),], aes(x= time, y = abundance)) + geom_line(aes(colour = silo)) + facet_wrap(~protein)

#make abundance plots of non-normalized values
cutDnonNorm_long <- tidyr::gather(cutDnonNorm, "protein","abundance", 3:ncol(cutDnonNorm))
cutDnonNorm_long <- merge(cutDnonNorm_long, d)
ggplot(cutDnonNorm_long, aes(x= Time, y = abundance)) + 
  geom_line(aes(colour = as.character(silo))) + facet_wrap(~protein, ncol = 10, scales = "free") +
  theme(strip.text.x = element_text(size = 5)) + ggtitle("all proteins with temperature PC1 loadings absolute values >= 0.04")

ggplot(cutDnonNorm_long[which(cutDnonNorm_long$PC1loadings > 0),], aes(x= Time, y = abundance)) + 
  geom_line(aes(colour = as.character(silo))) + facet_wrap(~protein, ncol = 10, scales = "free") +
  theme(strip.text.x = element_text(size = 5)) + ggtitle("all proteins with temperature PC1 loadings values >= 0.04")

ggplot(cutDnonNorm_long[which(cutDnonNorm_long$PC1loadings < 0),], aes(x= Time, y = abundance)) + 
  geom_line(aes(colour = as.character(silo))) + facet_wrap(~protein, ncol = 10, scales = "free") +
  theme(strip.text.x = element_text(size = 5)) + ggtitle("all proteins with temperature PC1 loadings values =< -0.04")

#uniprot accessions
uniprot <- read.table("~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_shellytest/background-proteome-accession-no-pipes_UniProtIDs.txt", header = FALSE, sep = "\t")
#check for no duplicates
#> nrow(uniprot)
#[1] 26535
#> length(unique(uniprot$V1))
#[1] 26535

colnames(uniprot) <- c("protein", "uniprot")

uniprot_PC1_0.04 <- merge(cutd, uniprot, all.x = TRUE)
write.csv(uniprot_PC1_0.04,"~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_shellytest/uniprot_PC1_0.04.csv", 
          row.names = FALSE, quote = FALSE)

#make a list of peptides with the same uniprot id
dups <- data.frame(uniprot_PC1_0.04[duplicated(uniprot_PC1_0.04$uniprot),])
dups <- dups[!is.na(dups$uniprot),]
dups <- uniprot_PC1_0.04[which(uniprot_PC1_0.04$uniprot %in% dups$uniprot),]
dups$rename <- paste(dups$uniprot,dups$protein, sep = "_")

dups_cutDnonNormLong <- merge(cutDnonNorm_long, dups, by = "protein")

ggplot(dups_cutDnonNormLong, aes(x= Time, y = abundance)) + 
  geom_line(aes(colour = as.character(silo))) + facet_wrap(~rename, ncol = 4, scales = "free") +
  theme(strip.text.x = element_text(size = 5)) + ggtitle("all proteins with temperature PC1 loadings values =< -0.04")
