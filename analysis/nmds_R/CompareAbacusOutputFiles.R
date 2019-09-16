install.packages("arsenal")
library(arsenal)
library(gtools)


#Compare 02/14/2017 data with Sean's march 1 data
data_SR <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output021417.tsv", sep = "\t" , header=TRUE, stringsAsFactors = FALSE)
data_SB <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_outputMar1.tsv", sep = "\t" , header=TRUE, stringsAsFactors = FALSE)
compare(data_SR,data_SB)
###shows no difference in files
  #Compare Object
  #Function Call: 
  #  compare.data.frame(x = data_SR, y = data_SB)
  
  #Shared: 457 variables and 8443 observations.
  #Not shared: 0 variables and 0 observations.
  
  #Differences found in 0/456 variables compared.
  #0 variables compared have non-identical attributes.

#confirmed by command line diff
  #D-10-18-212-233:Desktop Shelly$ diff ~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_outputMar1.tsv ~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output021417.tsv 
  #D-10-18-212-233:Desktop Shelly$ 
    

# determine what said 'NSAF' values are
data_SR_NSAF <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output021417NSAF.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
data_SB_NUMSPECADJ <- data_SB[,c(1,grep("NUMSPECSADJ", colnames(data_SB)))]
colnames(data_SB_NUMSPECADJ) <- gsub("NUMSPECSADJ","ADJNSAF", colnames(data_SB_NUMSPECADJ))
compare(data_SR_NSAF,data_SB_NUMSPECADJ)
#Compare Object

#Function Call: 
  #compare.data.frame(x = data_SR_NSAF, y = data_SB_NUMSPECADJ)

#Shared: 46 variables and 8443 observations.
#Not shared: 0 variables and 0 observations.

#Differences found in 0/45 variables compared.
#0 variables compared have non-identical attributes.


#What is Kaitlyn's data?
data_KM <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUSdata_only.csv", header = TRUE, stringsAsFactors = FALSE)

#convert steven's column names
data_SR_NSAF_avg <- data_SR_NSAF[,c(1,grep("NSAF", colnames(data_SR_NSAF)))]
colnames(data_SR_NSAF_avg) <- gsub(pattern = "X20161205_SAMPLE_", "", colnames(data_SR_NSAF_avg))
colnames(data_SR_NSAF_avg) <- gsub(pattern = "_ADJNSAF", "", colnames(data_SR_NSAF_avg))

data_SR_NSAF_avg <- data_SR_NSAF_avg[,c("PROTID",mixedsort(colnames(data_SR_NSAF_avg[,-1])))]
#find avg of each technical rep (https://stackoverflow.com/questions/13739243/average-pairs-of-columns-in-r)
data_SR_NSAFavg <- data.frame(sapply(seq(2,ncol(data_SR_NSAF_avg),2), function(i) {
    rowMeans(data_SR_NSAF_avg[,c(i, i+1)], na.rm=T)
  }))

#create column names similar to Kaitlyn's file 
colnames(data_SR_NSAFavg) <- mixedsort(colnames(data_SR_NSAF_avg[,c(-1,-grep("A",colnames(data_SR_NSAF_avg)))]))
meta_data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/Rhonda_new_sample_names.csv", header = TRUE, stringsAsFactors = FALSE)
meta_data$silo <- substr(meta_data$Contents,5,5)
meta_data$day <- substr(meta_data$SampleName,5,6)

new_colnames <- list()
for (i in 1:length(colnames(data_SR_NSAFavg))){
  sample <- colnames(data_SR_NSAFavg)[i]
  new_colnames[[i]] <- paste(meta_data[meta_data$SampleID == sample,"silo"],meta_data[meta_data$SampleID == sample, "day"], sep = "_")
  }

colnames(data_SR_NSAFavg) <- new_colnames
  
data_SR_NSAF_avg <- cbind(data.frame(data_SR_NSAF_avg[,"PROTID"],stringsAsFactors = FALSE), data_SR_NSAFavg)

colnames(data_KM)
#[1] "Protein.ID"        "CompetentLarvae_1" "X2_3"              "X3_3"              "X9_3"             
#[6] "X2_5"              "X3_5"              "X9_5"              "X2_7"              "X3_7"             
#[11] "X9_7"              "X2_9"              "X3_9"              "X9_9"              "X2_11"            
#[16] "X3_11"             "X9_11"             "X2_13"             "X3_13"             "X9_13"            
#[21] "X2_15"             "X3_15"             "X9_15"            
colnames(data_SR_NSAF_avg)
#[1] "data_SR_NSAF_avg$PROTID" "e_0"                     "2_3"                     "3_3"                    
#[5] "9_3"                     "2_5"                     "3_5"                     "9_5"                    
#[9] "2_7"                     "3_7"                     "9_7"                     "2_9"                    
#[13] "3_9"                     "9_9"                     "2_11"                    "3_11"                   
#[17] "9_11"                    "2_13"                    "3_13"                    "9_13"                   
#[21] "2_15"                    "3_15"                    "9_15"                   

colnames(data_SR_NSAF_avg) <- colnames(data_KM)

str(data_KM)
str(data_SR_NSAF_avg)

#sort files so they are in the same order
data_KM <- data_KM[order(data_KM$Protein.ID),]
data_SR_NSAF_avg <- data_SR_NSAF_avg[order(data_SR_NSAF_avg$Protein.ID),]


compare(data_KM,data_SR_NSAF_avg)
