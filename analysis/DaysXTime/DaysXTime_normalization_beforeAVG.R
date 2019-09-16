Load Abacus data, parse out ADJNSAF values, and simplify column names to just sample number
```{r}
#upload data file
ABACUSdata <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output021417.tsv", sep = "\t", header=TRUE, stringsAsFactors = FALSE)
#select only columns containing ADJNSAF and Protein ID
ABACUSdata <- ABACUSdata[,c(1,grep("ADJNSAF", colnames(ABACUSdata)))]

## change column names in ABACUSdata to just sampleID
colnames(ABACUSdata) <- gsub(pattern = "X20161205_SAMPLE_", "", colnames(ABACUSdata))
colnames(ABACUSdata) <- gsub(pattern = "_ADJNSAF", "", colnames(ABACUSdata))
```

Load meta data file with temperature and day information
```{r}
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
```{r}
#Transpose- switch rows and columns
tABACUSdata <- t.data.frame(ABACUSdata[,-1])
colnames(tABACUSdata) <- ABACUSdata[,1]
tABACUSdata <- cbind(data.frame(rownames(tABACUSdata)),tABACUSdata)
colnames(tABACUSdata)[1] <- "SampleID"

#add meta data to abacus data
tABACUSdata <- merge(meta_data[,c(1,2,7,8)],tABACUSdata, by = "SampleID")
#add 0.1 to everything
tABACUSdata[,5:ncol(tABACUSdata)] <- tABACUSdata[,5:ncol(tABACUSdata)] + 0.1


tABACUSdata$day <- as.numeric(tABACUSdata$day)
tABACUSdata$temp <- as.numeric(tABACUSdata$temp)

tempXdays <- data.frame(tABACUSdata$day * tABACUSdata$temp)
tempXdays[tempXdays==0]<- 16


tab_normTxD <- data.frame()

for(i in 1:nrow(tABACUSdata)){
  row <- tABACUSdata[i,5:ncol(tABACUSdata)]/(tempXdays[i,1])
  tab_normTxD <- rbind(tab_normTxD, row)
}

tABACUSdata_norm <- cbind(tABACUSdata[,1:4],tab_normTxD)

pca
```{r}
pca <- prcomp(log(tABACUSdata_norm[,-c(1:4)],2))
pca_meta <- cbind(tABACUSdata_norm[,c(3:4)],pca$x)
pca_meta$day <- factor(pca_meta$day)
pca_meta$temp <- factor(pca_meta$temp)

ggplot(pca_meta, aes(PC1, PC2)) + geom_point(aes(col = day, shape = temp)) + theme_bw() + ggtitle("PCA of ADJNSAF values scaled by day*time")

```

pca <- prcomp(log(tABACUSdata_norm[which(substr(tABACUSdata_norm$SampleName,0,2) != "S2"& substr(tABACUSdata_norm$SampleName,0,2) != "po"),-c(1:4)],2))
pca_meta <- cbind(tABACUSdata_norm[which(substr(tABACUSdata_norm$SampleName,0,2) != "S2"& substr(tABACUSdata_norm$SampleName,0,2) != "po"),c(3:4)],pca$x)
pca_meta$day <- factor(pca_meta$day)
pca_meta$temp <- factor(pca_meta$temp)

ggplot(pca_meta, aes(PC1, PC2)) + geom_point(aes(col = day, shape = temp)) + theme_bw() + ggtitle("PCA of ADJNSAF values scaled by day*time")




pca <- prcomp(tABACUSdata_norm[which(substr(tABACUSdata_norm$SampleName,0,2) != "S2"& substr(tABACUSdata_norm$SampleName,0,2) != "po" & tABACUSdata_norm$day !=3),-c(1:4)])
pca_meta <- cbind(tABACUSdata_norm[which(substr(tABACUSdata_norm$SampleName,0,2) != "S2"& substr(tABACUSdata_norm$SampleName,0,2) != "po" & tABACUSdata_norm$day !=3),c(3:4)],pca$x)
pca_meta$day <- factor(pca_meta$day)
pca_meta$temp <- factor(pca_meta$temp)

ggplot(pca_meta, aes(PC1, PC2)) + geom_point(aes(col = day, shape = temp)) + theme_bw() + ggtitle("PCA of ADJNSAF values scaled by day*time")
