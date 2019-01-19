Clustering technical replicates
================
Shelly Trigg
1/11/2019

Load packages

``` r
library(vegan)
```

    ## Warning: package 'vegan' was built under R version 3.4.4

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-3

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 3.4.4

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 3.4.4

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(gtools)
```

    ## Warning: package 'gtools' was built under R version 3.4.4

    ## 
    ## Attaching package: 'gtools'

    ## The following object is masked from 'package:permute':
    ## 
    ##     permute

Load Abacus data, parse out ADJNSAF values, and simplify column names to just sample number

``` r
#upload data file
ABACUSdata <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output021417.tsv", sep = "\t", header=TRUE, stringsAsFactors = FALSE)
#select only columns containing ADJNSAF and Protein ID
ABACUSdata <- ABACUSdata[,c(1,grep("ADJNSAF", colnames(ABACUSdata)))]

## change column names in ABACUSdata to just sampleID
colnames(ABACUSdata) <- gsub(pattern = "X20161205_SAMPLE_", "", colnames(ABACUSdata))
colnames(ABACUSdata) <- gsub(pattern = "_ADJNSAF", "", colnames(ABACUSdata))
```

Load meta data file with temperature and day information

``` r
#upload meta data; this was a csv file I create from Rhonda's notebook entry: https://github.com/Ellior2/Ellior2.github.io/blob/master/_posts/2017-3-11-NMDS-analysis.md
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
```

Reformat Abacus data for NMDS

``` r
#Transpose- switch rows and columns
tABACUSdata <- t.data.frame(ABACUSdata[,-1])
colnames(tABACUSdata) <- ABACUSdata[,1]
tABACUSdata <- cbind(data.frame(rownames(tABACUSdata)),tABACUSdata)
colnames(tABACUSdata)[1] <- "SampleID"

#add meta data to abacus data
tABACUSdata <- merge(meta_data[,c(1,2,7,8)],tABACUSdata, by = "SampleID")

#Remove Silo 9 
silo3and2 <- tABACUSdata[which(substr(tABACUSdata$SampleName,1,2) != "S9"),]
#make rownames from Sample ID column so that the NMDS knows what's what
rownames(silo3and2) <- silo3and2$SampleID
#order the data frame by day and temperature so coloring the points on the plot is easier
silo3and2 <- silo3and2[order(as.numeric(silo3and2$day),silo3and2$temp),]
```

Determine if any proteins have zero ADJNSAF vals for all samples; this would be because they were in Silo 2, but not in Silo 3 or 9

``` r
no_val_proteins <- silo3and2[,which(apply(silo3and2, 2, var) == 0)]
```

    ## Warning in FUN(newX[, i], ...): NAs introduced by coercion

    ## Warning in FUN(newX[, i], ...): NAs introduced by coercion

``` r
ncol(no_val_proteins)
```

    ## [1] 369

Remove proteins if they have a zero value in all samples

``` r
silo3and2_nozerovar <- silo3and2[,-c(1:4,which(colnames(silo3and2) %in% colnames(no_val_proteins)))]
#check to make sure it worked
ncol(silo3and2)-ncol(silo3and2_nozerovar)
```

    ## [1] 373

For proteins with a zero value in any sample, replace with very small value

``` r
silo3and2_nozerovar[silo3and2_nozerovar == 0.0000] <- 0.1000
```

try PCA

``` r
#first convert day and temp data to factors for plotting
silo3and2$day <- factor(silo3and2$day, levels = unique(silo3and2$day))
silo3and2$temp <- factor(silo3and2$temp, levels = unique(silo3and2$temp))
silo3and2$silo <- factor(substr(silo3and2$SampleName,0,2))

#run PCA
pca <- prcomp(silo3and2_nozerovar, center = T, scale = T)
pca_meta <- data.frame(cbind(silo3and2[,c("day","temp","silo")],pca$x))
ggplot(pca_meta, aes(PC1, PC2)) + geom_point(aes(col = day, shape = silo)) + theme_bw() + ggtitle("PCA of ADJNSAF values where zeros were replaced with 0.1")
```

![](ClusteringTechnicalReplicatesSilo2vs3_files/figure-markdown_github/unnamed-chunk-9-1.png)

try PCA on log transformed values

``` r
silo3and2_log <- log(silo3and2_nozerovar,2)
pca_log <- prcomp(silo3and2_log, center = F, scale = F)
pca_log_meta <-  data.frame(cbind(silo3and2[,c("day","temp","silo")],pca$x))
colnames(pca_log_meta)[1:3] <- c("day","temp","silo")
ggplot(pca_log_meta, aes(PC1, PC2)) + geom_point(aes(col = day, shape = silo)) + theme_bw() + ggtitle("PCA of log ADJNSAF values with zeros replaced with 0.1")
```

![](ClusteringTechnicalReplicatesSilo2vs3_files/figure-markdown_github/unnamed-chunk-10-1.png)

boxplots of common contaminant proteins; trypsin shows generally consistent levels as we'd expect

``` r
boxplot(silo3and2_nozerovar[,-grep("CHOYP", colnames(silo3and2_nozerovar))])
```

![](ClusteringTechnicalReplicatesSilo2vs3_files/figure-markdown_github/unnamed-chunk-11-1.png)

Average NSAF values for all proteins

``` r
df_all_avg <- data.frame()
#loop through the data and calculate the mean between replicates for each protein
for (i in seq(1,nrow(silo3and2_nozerovar),2)){
  #this calculates the mean for each odd number row and the row following it
  df_all_avg_row <- apply(silo3and2_nozerovar[c(i,i+1),],2,mean)
  #this sequencially combines rows of data together after the SD is generated
  df_all_avg <- rbind(df_all_avg, df_all_avg_row)
}
#add column names to SD data
colnames(df_all_avg) <- colnames(silo3and2_nozerovar)
rownames(df_all_avg) <- rownames(silo3and2_nozerovar[-grep("A",rownames(silo3and2_nozerovar)),])
```

export data set avg tech. rep. NSAF for all proteins
