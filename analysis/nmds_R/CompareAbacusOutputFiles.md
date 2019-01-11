CompareAbacusOutputFiles
================
Shelly Trigg
1/11/2019

``` r
install.packages("arsenal")
```

``` r
library(arsenal)
library(gtools)
```

    ## Warning: package 'gtools' was built under R version 3.4.4

Compare 02/14/2017 data with Sean's march 1 data

``` r
#load data
data_SR <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output021417.tsv", sep = "\t" , header=TRUE, stringsAsFactors = FALSE)
data_SB <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output.tsv", sep = "\t" , header=TRUE, stringsAsFactors = FALSE)
```

``` r
compare(data_SR,data_SB)
```

    ## Compare Object
    ## 
    ## Function Call: 
    ## compare.data.frame(x = data_SR, y = data_SB)
    ## 
    ## Shared: 457 variables and 8443 observations.
    ## Not shared: 0 variables and 0 observations.
    ## 
    ## Differences found in 0/456 variables compared.
    ## 0 variables compared have non-identical attributes.

R compare function shows no difference between files and I confirmed this with the command line function diff, which gave no output

D-10-18-212-233:Desktop Shelly$ diff ~/Documents/GitHub/OysterSeedProject/raw\_data/ABACUS\_outputMar1.tsv ~/Documents/GitHub/OysterSeedProject/raw\_data/ABACUS\_output021417.tsv
D-10-18-212-233:Desktop Shelly$

Determine what ABACUS\_output021417NSAF.tsv 'NSAF' values are

``` r
data_SR_NSAF <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output021417NSAF.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
data_SB_NUMSPECADJ <- data_SB[,c(1,grep("NUMSPECSADJ", colnames(data_SB)))]
colnames(data_SB_NUMSPECADJ) <- gsub("NUMSPECSADJ","ADJNSAF", colnames(data_SB_NUMSPECADJ))
compare(data_SR_NSAF,data_SB_NUMSPECADJ)
```

    ## Compare Object
    ## 
    ## Function Call: 
    ## compare.data.frame(x = data_SR_NSAF, y = data_SB_NUMSPECADJ)
    ## 
    ## Shared: 46 variables and 8443 observations.
    ## Not shared: 0 variables and 0 observations.
    ## 
    ## Differences found in 0/45 variables compared.
    ## 0 variables compared have non-identical attributes.

***NSAF values are in fact NUMSPECSADJ values***

Determine what the values are in Kaitlyn's data

``` r
#load data
data_KM <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUSdata_only.csv", header = TRUE, stringsAsFactors = FALSE)

#convert steven's column names to samlpe number only
data_SR_NSAF_avg <- data_SR_NSAF[,c(1,grep("NSAF", colnames(data_SR_NSAF)))]
colnames(data_SR_NSAF_avg) <- gsub(pattern = "X20161205_SAMPLE_", "", colnames(data_SR_NSAF_avg))
colnames(data_SR_NSAF_avg) <- gsub(pattern = "_ADJNSAF", "", colnames(data_SR_NSAF_avg))

#sort the data frame by order of samples (e.g. 1, 1A, 3, 3A, etc.) but keeping ProteinID as column 1
data_SR_NSAF_avg <- data_SR_NSAF_avg[,c("PROTID",mixedsort(colnames(data_SR_NSAF_avg[,-1])))]
```

Based on a few values in Kaitlyn's data, we think they are averages of the values in ABACUS\_output021417NSAF.tsv, so to confirm this we'll calculate the avg between technical reps in ABACUS\_output021417NSAF.tsv

``` r
#calculate avg between technical reps. This makes a new data frame of just sample averages (https://stackoverflow.com/questions/13739243/average-pairs-of-columns-in-r)
data_SR_NSAFavg <- data.frame(sapply(seq(2,ncol(data_SR_NSAF_avg),2), function(i) {
    rowMeans(data_SR_NSAF_avg[,c(i, i+1)], na.rm=T)
  }))

#create column names similar to Kaitlyn's file so the two data sets can be compared
colnames(data_SR_NSAFavg) <- mixedsort(colnames(data_SR_NSAF_avg[,c(-1,-grep("A",colnames(data_SR_NSAF_avg)))]))

#load meta data to convert sample number to silo_day format
meta_data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/Rhonda_new_sample_names.csv", header = TRUE, stringsAsFactors = FALSE)
#make column for silo
meta_data$silo <- substr(meta_data$Contents,5,5)
#make column for day
meta_data$day <- substr(meta_data$SampleName,5,6)
#create a new list of column names by converting from sample number to silo_day
new_colnames <- list()
#This loops through each sample number, translates it into silo_day, and adds it to the new list of column names
for (i in 1:length(colnames(data_SR_NSAFavg))){
  sample <- colnames(data_SR_NSAFavg)[i]
  new_colnames[[i]] <- paste(meta_data[meta_data$SampleID == sample,"silo"],meta_data[meta_data$SampleID == sample, "day"], sep = "_")
  }
#add new silo_day column names to replicate averages
colnames(data_SR_NSAFavg) <- new_colnames
  
#make a new data frame with ProteinID and replicate averages
data_SR_NSAF_avg <- cbind(data.frame(data_SR_NSAF_avg[,"PROTID"],stringsAsFactors = FALSE), data_SR_NSAFavg)
```

compare Kaitlyn's column names with data\_SR\_NSAF\_avg

``` r
colnames(data_KM)
```

    ##  [1] "Protein.ID"        "CompetentLarvae_1" "X2_3"             
    ##  [4] "X3_3"              "X9_3"              "X2_5"             
    ##  [7] "X3_5"              "X9_5"              "X2_7"             
    ## [10] "X3_7"              "X9_7"              "X2_9"             
    ## [13] "X3_9"              "X9_9"              "X2_11"            
    ## [16] "X3_11"             "X9_11"             "X2_13"            
    ## [19] "X3_13"             "X9_13"             "X2_15"            
    ## [22] "X3_15"             "X9_15"

``` r
colnames(data_SR_NSAF_avg)
```

    ##  [1] "data_SR_NSAF_avg....PROTID.." "e_0"                         
    ##  [3] "2_3"                          "3_3"                         
    ##  [5] "9_3"                          "2_5"                         
    ##  [7] "3_5"                          "9_5"                         
    ##  [9] "2_7"                          "3_7"                         
    ## [11] "9_7"                          "2_9"                         
    ## [13] "3_9"                          "9_9"                         
    ## [15] "2_11"                         "3_11"                        
    ## [17] "9_11"                         "2_13"                        
    ## [19] "3_13"                         "9_13"                        
    ## [21] "2_15"                         "3_15"                        
    ## [23] "9_15"

They are in the same order so I added Kaitlyn's column names to data\_SR\_NSAF\_avg for easier comparison of files

``` r
colnames(data_SR_NSAF_avg) <- colnames(data_KM)
```

Compare column classes between files so files can be easily compared

``` r
str(data_KM)
```

    ## 'data.frame':    8443 obs. of  23 variables:
    ##  $ Protein.ID       : chr  "CHOYP_EF1A.2.4|m.50858" "CHOYP_MYS.4.7|m.18188" "CHOYP_PHUM_PHUM226120.7.7|m.54053" "CHOYP_MATN1.2.5|m.45361" ...
    ##  $ CompetentLarvae_1: num  268 254 179 186 210 ...
    ##  $ X2_3             : num  292 311 222 161 230 ...
    ##  $ X3_3             : num  296 302 250 156 222 ...
    ##  $ X9_3             : num  349 318 258 211 91 ...
    ##  $ X2_5             : num  254 298 222 194 180 ...
    ##  $ X3_5             : num  332 310 251 202 91 ...
    ##  $ X9_5             : num  320 298 244 188 268 ...
    ##  $ X2_7             : num  300 284 237 181 228 ...
    ##  $ X3_7             : num  294 288 240 179 172 ...
    ##  $ X9_7             : num  281.5 259.5 201.5 160.5 88.5 ...
    ##  $ X2_9             : num  262 252 186 177 222 ...
    ##  $ X3_9             : num  298 276 215 203 245 ...
    ##  $ X9_9             : num  356 325 276 236 192 ...
    ##  $ X2_11            : num  144 138 60 120 178 ...
    ##  $ X3_11            : num  358 293.5 245.5 226 80.5 ...
    ##  $ X9_11            : num  325 344 278 206.5 83.5 ...
    ##  $ X2_13            : num  156 160 128 121 185 ...
    ##  $ X3_13            : num  294.5 243.5 42.5 175.5 254.5 ...
    ##  $ X9_13            : num  341 379 307 200 130 ...
    ##  $ X2_15            : num  250 116 188 196 236 ...
    ##  $ X3_15            : num  288 224 130 189 258 ...
    ##  $ X9_15            : num  349 360 300 230 224 ...

``` r
str(data_SR_NSAF_avg)
```

    ## 'data.frame':    8443 obs. of  23 variables:
    ##  $ Protein.ID       : chr  "CHOYP_1433G.2.2|m.63450" "CHOYP_AADAT.1.1|m.12672" "CHOYP_ALF.3.3|m.66837" "CHOYP_FUT8.1.1|m.9231" ...
    ##  $ CompetentLarvae_1: num  27.5 3 0 3.5 11.5 20 1 0 5 4.5 ...
    ##  $ X2_3             : num  31 2.5 0 1.5 8 8 2.5 1 2.5 0.5 ...
    ##  $ X3_3             : num  28 2 0 0.5 7.5 21.5 2.5 0 4.5 1 ...
    ##  $ X9_3             : num  26 2.5 0 1 10 18 1.5 1.5 5.5 1 ...
    ##  $ X2_5             : num  27 3.5 0 0.5 8.5 13.5 1.5 1.5 2.5 1 ...
    ##  $ X3_5             : num  27 2.5 0 0.5 8.5 18 1 2 3 0 ...
    ##  $ X9_5             : num  30 3 0 0 8 16 2 1 3 0 ...
    ##  $ X2_7             : num  29.5 2.5 0 3.5 7 15.5 1.5 0.5 3 2.5 ...
    ##  $ X3_7             : num  27 2.5 0 0.5 9 14 1.5 2 5.5 1 ...
    ##  $ X9_7             : num  28 2 0.5 1.5 7.5 8.5 0.5 0.5 1.5 0 ...
    ##  $ X2_9             : num  26.5 2 0.5 2 7.5 12.5 0.5 0 1 2 ...
    ##  $ X3_9             : num  37.5 3 0.5 0 9 17 1.5 1 4.5 0.5 ...
    ##  $ X9_9             : num  30.5 3.5 1.5 0 9.5 14 1.5 0 4 1 ...
    ##  $ X2_11            : num  12.5 1 0 0 1.5 2.5 0 0 0 0 ...
    ##  $ X3_11            : num  30.5 4.5 2.5 0.5 9 14 2 0 2.5 1 ...
    ##  $ X9_11            : num  31.5 2 1 0 11 16.5 2 0 4.5 1 ...
    ##  $ X2_13            : num  24.5 1.5 0 0 7.5 6 0.5 0 0 0 ...
    ##  $ X3_13            : num  23.5 3 2.5 3 9 16 0 0 0.5 4 ...
    ##  $ X9_13            : num  31.5 2 2.5 1.5 8.5 22.5 2.5 0 2 2.5 ...
    ##  $ X2_15            : num  27 3 3 0.5 8.5 11 0 0 2.5 0.5 ...
    ##  $ X3_15            : num  27.5 2.5 4 1.5 11.5 19.5 1.5 0 0.5 1.5 ...
    ##  $ X9_15            : num  35 2 2 0.5 10 22 2.5 0 4.5 1 ...

sort files so they are in the same order

``` r
data_KM <- data_KM[order(data_KM$Protein.ID),]
data_SR_NSAF_avg <- data_SR_NSAF_avg[order(data_SR_NSAF_avg$Protein.ID),]
```

finally, compare files

``` r
compare(data_KM,data_SR_NSAF_avg)
```

    ## Compare Object
    ## 
    ## Function Call: 
    ## compare.data.frame(x = data_KM, y = data_SR_NSAF_avg)
    ## 
    ## Shared: 24 variables and 8443 observations.
    ## Not shared: 0 variables and 0 observations.
    ## 
    ## Differences found in 0/23 variables compared.
    ## 0 variables compared have non-identical attributes.

***files are the same so Kaitlyn's data is in fact the average of the NUMSPECSADJ***
