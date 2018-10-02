library(dplyr)
library(tidyr)
library(MetStaT)


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



