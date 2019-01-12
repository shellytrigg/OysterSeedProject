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

#Remove Silo 2 and day 15
silo3and9 <- tABACUSdata[which(substr(tABACUSdata$SampleName,1,2) != "S2" & tABACUSdata$day != "15"),]
#make rownames from Sample ID column so that the NMDS knows what's what
rownames(silo3and9) <- silo3and9$SampleID
#order the data frame by day and temperature so coloring the points on the plot is easier
silo3and9 <- silo3and9[order(as.numeric(silo3and9$day),silo3and9$temp),]
```

Determine if any proteins have zero ADJNSAF vals for all samples; this would be because they were in Silo 2, but not in Silo 3 or 9

``` r
no_val_proteins <- silo3and9[,which(apply(silo3and9, 2, var) == 0)]
```

    ## Warning in FUN(newX[, i], ...): NAs introduced by coercion

    ## Warning in FUN(newX[, i], ...): NAs introduced by coercion

``` r
ncol(no_val_proteins)
```

    ## [1] 451

Remove proteins if they have a zero value in all samples

``` r
silo3and9_nozerovar <- silo3and9[,-c(1:4,which(colnames(silo3and9) %in% colnames(no_val_proteins)))]
#check to make sure it worked
ncol(silo3and9)-ncol(silo3and9_nozerovar)
```

    ## [1] 455

For proteins with a zero value in any sample, replace with very small value

``` r
silo3and9_nozerovar[silo3and9_nozerovar == 0.0000] <- 0.1000
```

try PCA

``` r
pca <- prcomp(silo3and9_nozerovar, center = T, scale = T)
pca_meta <- cbind(silo3and9$day, silo3and9$temp, data.frame(paste(silo3and9$day,silo3and9$temp, sep = "_")),pca$x)
colnames(pca_meta)[1:3] <- c("day","temp","SampleName")
ggplot(pca_meta, aes(PC1, PC2)) + geom_point(aes(col = day, shape = temp)) + theme_bw() + ggtitle("PCA of ADJNSAF values where zeros were replaced with 0.1")
```

![](ClusteringTechnicalReplicates_files/figure-markdown_github/unnamed-chunk-8-1.png)

try PCA on log transformed values

``` r
silo3and9_log <- log(silo3and9_nozerovar,2)
pca_log <- prcomp(silo3and9_log, center = F, scale = F)
pca_log_meta <- cbind(silo3and9$day, silo3and9$temp, data.frame(paste(silo3and9$day,silo3and9$temp, sep = "_")),pca_log$x)
colnames(pca_log_meta)[1:3] <- c("day","temp","SampleName")
ggplot(pca_log_meta, aes(PC1, PC2)) + geom_point(aes(col = day, shape = temp)) + theme_bw() + ggtitle("PCA of log ADJNSAF values with zeros replaced with 0.1")
```

![](ClusteringTechnicalReplicates_files/figure-markdown_github/unnamed-chunk-9-1.png)

Make MDS dissimilarity matrix

``` r
nmds.silo3and9 <- metaMDS(silo3and9_nozerovar, distance = 'euclidean', k = 2, trymax = 3000, autotransform = FALSE)
```

    ## Run 0 stress 0.1649526 
    ## Run 1 stress 0.1832198 
    ## Run 2 stress 0.1725523 
    ## Run 3 stress 0.1631591 
    ## ... New best solution
    ## ... Procrustes: rmse 0.03627924  max resid 0.130967 
    ## Run 4 stress 0.1649514 
    ## Run 5 stress 0.2229974 
    ## Run 6 stress 0.1808447 
    ## Run 7 stress 0.1765492 
    ## Run 8 stress 0.1649514 
    ## Run 9 stress 0.2060802 
    ## Run 10 stress 0.2374735 
    ## Run 11 stress 0.1763841 
    ## Run 12 stress 0.1725524 
    ## Run 13 stress 0.2313873 
    ## Run 14 stress 0.1649514 
    ## Run 15 stress 0.1649515 
    ## Run 16 stress 0.1763242 
    ## Run 17 stress 0.1763798 
    ## Run 18 stress 0.1807039 
    ## Run 19 stress 0.1631588 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001771228  max resid 0.000482403 
    ## ... Similar to previous best
    ## Run 20 stress 0.1631587 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002149395  max resid 0.0007396403 
    ## ... Similar to previous best
    ## *** Solution reached

``` r
#make data frame of NMDS scores
nmds.silo3and9.scores <- cbind(silo3and9$day, silo3and9$temp,data.frame(scores(nmds.silo3and9)))
colnames(nmds.silo3and9.scores)[1:2] <- c("day","temp")
ggplot(nmds.silo3and9.scores, aes(NMDS1, NMDS2)) + geom_point(aes(col = day, shape = temp)) + theme_bw() + ggtitle("NMDS of ADJNSAF values with zeros replaced with 0.1")
```

![](ClusteringTechnicalReplicates_files/figure-markdown_github/unnamed-chunk-10-1.png)

Make MDS dissimilarity matrix with log tranformed ADJNSAF values

``` r
nmds.silo3and9_log <- metaMDS(silo3and9_log, distance = 'euclidean', k = 2, trymax = 3000, autotransform = FALSE)
```

    ## 'comm' has negative data: 'autotransform', 'noshare' and 'wascores' set to FALSE

    ## Run 0 stress 0.1123882 
    ## Run 1 stress 0.1229455 
    ## Run 2 stress 0.1273254 
    ## Run 3 stress 0.1229465 
    ## Run 4 stress 0.1229456 
    ## Run 5 stress 0.1124372 
    ## ... Procrustes: rmse 0.004791233  max resid 0.01675094 
    ## Run 6 stress 0.1301972 
    ## Run 7 stress 0.1123886 
    ## ... Procrustes: rmse 0.0001823376  max resid 0.0005722986 
    ## ... Similar to previous best
    ## Run 8 stress 0.1124373 
    ## ... Procrustes: rmse 0.004845089  max resid 0.01677786 
    ## Run 9 stress 0.1124367 
    ## ... Procrustes: rmse 0.004808159  max resid 0.01674914 
    ## Run 10 stress 0.1124367 
    ## ... Procrustes: rmse 0.004751571  max resid 0.01673026 
    ## Run 11 stress 0.1155856 
    ## Run 12 stress 0.1124365 
    ## ... Procrustes: rmse 0.004782271  max resid 0.01674964 
    ## Run 13 stress 0.125918 
    ## Run 14 stress 0.1124368 
    ## ... Procrustes: rmse 0.004855167  max resid 0.01696226 
    ## Run 15 stress 0.1229455 
    ## Run 16 stress 0.1124367 
    ## ... Procrustes: rmse 0.004813387  max resid 0.01677408 
    ## Run 17 stress 0.1285269 
    ## Run 18 stress 0.1229459 
    ## Run 19 stress 0.112438 
    ## ... Procrustes: rmse 0.004817579  max resid 0.01680181 
    ## Run 20 stress 0.125918 
    ## *** Solution reached

``` r
#make data frame of NMDS scores
nmds.silo3and9_log.scores <- cbind(silo3and9$day, silo3and9$temp,data.frame(scores(nmds.silo3and9_log)))
colnames(nmds.silo3and9_log.scores)[1:2] <- c("day","temp")
ggplot(nmds.silo3and9_log.scores, aes(NMDS1, NMDS2)) + geom_point(aes(col = day, shape = temp)) + theme_bw() + ggtitle("NMDS of log ADJNSAF values with zeros replaced with 0.1")
```

![](ClusteringTechnicalReplicates_files/figure-markdown_github/unnamed-chunk-11-1.png)

Make MDS dissimilarity matrix with log transformed ADJNSAF values and bray curtis distance

``` r
nmds.silo3and9_log_bray <- metaMDS(silo3and9_log, distance = 'bray', k = 2, trymax = 3000, autotransform = FALSE)
```

    ## 'comm' has negative data: 'autotransform', 'noshare' and 'wascores' set to FALSE

    ## Warning in distfun(comm, method = distance, ...): results may be
    ## meaningless because data have negative entries in method "bray"

    ## Run 0 stress 0.2389081 
    ## Run 1 stress 0.2370844 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1201977  max resid 0.4665152 
    ## Run 2 stress 0.2375143 
    ## ... Procrustes: rmse 0.09418388  max resid 0.2744619 
    ## Run 3 stress 0.2418758 
    ## Run 4 stress 0.2362341 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1131596  max resid 0.4707707 
    ## Run 5 stress 0.231722 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1375828  max resid 0.3820759 
    ## Run 6 stress 0.2394878 
    ## Run 7 stress 0.2351632 
    ## Run 8 stress 0.2378052 
    ## Run 9 stress 0.2387708 
    ## Run 10 stress 0.2361597 
    ## Run 11 stress 0.2387661 
    ## Run 12 stress 0.2403012 
    ## Run 13 stress 0.2332155 
    ## Run 14 stress 0.2386243 
    ## Run 15 stress 0.2410579 
    ## Run 16 stress 0.2371749 
    ## Run 17 stress 0.2404129 
    ## Run 18 stress 0.2362204 
    ## Run 19 stress 0.238277 
    ## Run 20 stress 0.2399971 
    ## Run 21 stress 0.2389445 
    ## Run 22 stress 0.2311166 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1324409  max resid 0.3420728 
    ## Run 23 stress 0.2322179 
    ## Run 24 stress 0.2392819 
    ## Run 25 stress 0.2420451 
    ## Run 26 stress 0.2341674 
    ## Run 27 stress 0.2410797 
    ## Run 28 stress 0.2360711 
    ## Run 29 stress 0.2430208 
    ## Run 30 stress 0.2341083 
    ## Run 31 stress 0.2367261 
    ## Run 32 stress 0.2351233 
    ## Run 33 stress 0.2402958 
    ## Run 34 stress 0.2371514 
    ## Run 35 stress 0.2376205 
    ## Run 36 stress 0.2357487 
    ## Run 37 stress 0.2375027 
    ## Run 38 stress 0.2389348 
    ## Run 39 stress 0.238102 
    ## Run 40 stress 0.2401938 
    ## Run 41 stress 0.2406571 
    ## Run 42 stress 0.2381084 
    ## Run 43 stress 0.2360213 
    ## Run 44 stress 0.237131 
    ## Run 45 stress 0.2373157 
    ## Run 46 stress 0.2435319 
    ## Run 47 stress 0.2390894 
    ## Run 48 stress 0.2371829 
    ## Run 49 stress 0.2384527 
    ## Run 50 stress 0.2376988 
    ## Run 51 stress 0.2391981 
    ## Run 52 stress 0.2393474 
    ## Run 53 stress 0.2389696 
    ## Run 54 stress 0.239224 
    ## Run 55 stress 0.2367726 
    ## Run 56 stress 0.2324585 
    ## Run 57 stress 0.2375605 
    ## Run 58 stress 0.2350286 
    ## Run 59 stress 0.2417161 
    ## Run 60 stress 0.2351907 
    ## Run 61 stress 0.2362617 
    ## Run 62 stress 0.2380472 
    ## Run 63 stress 0.2386046 
    ## Run 64 stress 0.2342985 
    ## Run 65 stress 0.23602 
    ## Run 66 stress 0.2367636 
    ## Run 67 stress 0.2386814 
    ## Run 68 stress 0.2381655 
    ## Run 69 stress 0.2361144 
    ## Run 70 stress 0.2376042 
    ## Run 71 stress 0.2374542 
    ## Run 72 stress 0.2392629 
    ## Run 73 stress 0.2369018 
    ## Run 74 stress 0.2329158 
    ## Run 75 stress 0.236133 
    ## Run 76 stress 0.2356022 
    ## Run 77 stress 0.2434633 
    ## Run 78 stress 0.2351909 
    ## Run 79 stress 0.2342961 
    ## Run 80 stress 0.2364713 
    ## Run 81 stress 0.242677 
    ## Run 82 stress 0.2394049 
    ## Run 83 stress 0.2334115 
    ## Run 84 stress 0.2336446 
    ## Run 85 stress 0.2409735 
    ## Run 86 stress 0.2374493 
    ## Run 87 stress 0.236863 
    ## Run 88 stress 0.2636069 
    ## Run 89 stress 0.2373976 
    ## Run 90 stress 0.2336447 
    ## Run 91 stress 0.2337848 
    ## Run 92 stress 0.2358668 
    ## Run 93 stress 0.2394398 
    ## Run 94 stress 0.2371692 
    ## Run 95 stress 0.2390253 
    ## Run 96 stress 0.2415529 
    ## Run 97 stress 0.2370376 
    ## Run 98 stress 0.2365855 
    ## Run 99 stress 0.2389674 
    ## Run 100 stress 0.2356966 
    ## Run 101 stress 0.2387204 
    ## Run 102 stress 0.2351315 
    ## Run 103 stress 0.2383464 
    ## Run 104 stress 0.2338752 
    ## Run 105 stress 0.2359935 
    ## Run 106 stress 0.235421 
    ## Run 107 stress 0.235818 
    ## Run 108 stress 0.2336445 
    ## Run 109 stress 0.2384478 
    ## Run 110 stress 0.2361334 
    ## Run 111 stress 0.2373454 
    ## Run 112 stress 0.2405019 
    ## Run 113 stress 0.2388798 
    ## Run 114 stress 0.237591 
    ## Run 115 stress 0.2355297 
    ## Run 116 stress 0.2333431 
    ## Run 117 stress 0.2367944 
    ## Run 118 stress 0.2353658 
    ## Run 119 stress 0.2338618 
    ## Run 120 stress 0.2413707 
    ## Run 121 stress 0.2364599 
    ## Run 122 stress 0.2387326 
    ## Run 123 stress 0.2397677 
    ## Run 124 stress 0.2382369 
    ## Run 125 stress 0.2353468 
    ## Run 126 stress 0.2425622 
    ## Run 127 stress 0.2365 
    ## Run 128 stress 0.2387528 
    ## Run 129 stress 0.2324366 
    ## Run 130 stress 0.2401792 
    ## Run 131 stress 0.2360555 
    ## Run 132 stress 0.2389988 
    ## Run 133 stress 0.2403115 
    ## Run 134 stress 0.2382962 
    ## Run 135 stress 0.2328418 
    ## Run 136 stress 0.2376678 
    ## Run 137 stress 0.2344606 
    ## Run 138 stress 0.2341201 
    ## Run 139 stress 0.2402366 
    ## Run 140 stress 0.2316003 
    ## ... Procrustes: rmse 0.1193441  max resid 0.3539816 
    ## Run 141 stress 0.233599 
    ## Run 142 stress 0.2402079 
    ## Run 143 stress 0.2389244 
    ## Run 144 stress 0.2392029 
    ## Run 145 stress 0.2392403 
    ## Run 146 stress 0.2365203 
    ## Run 147 stress 0.2350124 
    ## Run 148 stress 0.2370147 
    ## Run 149 stress 0.2380703 
    ## Run 150 stress 0.2379464 
    ## Run 151 stress 0.2374762 
    ## Run 152 stress 0.2343422 
    ## Run 153 stress 0.2369017 
    ## Run 154 stress 0.2368365 
    ## Run 155 stress 0.2425894 
    ## Run 156 stress 0.2391152 
    ## Run 157 stress 0.2336638 
    ## Run 158 stress 0.2362786 
    ## Run 159 stress 0.2370388 
    ## Run 160 stress 0.2361013 
    ## Run 161 stress 0.2347758 
    ## Run 162 stress 0.2378131 
    ## Run 163 stress 0.2371681 
    ## Run 164 stress 0.2375853 
    ## Run 165 stress 0.2375665 
    ## Run 166 stress 0.2347395 
    ## Run 167 stress 0.236478 
    ## Run 168 stress 0.239176 
    ## Run 169 stress 0.2370182 
    ## Run 170 stress 0.2414272 
    ## Run 171 stress 0.2378327 
    ## Run 172 stress 0.2365169 
    ## Run 173 stress 0.2373149 
    ## Run 174 stress 0.2366425 
    ## Run 175 stress 0.2356716 
    ## Run 176 stress 0.2375459 
    ## Run 177 stress 0.3682573 
    ## Run 178 stress 0.2395351 
    ## Run 179 stress 0.2410627 
    ## Run 180 stress 0.2372847 
    ## Run 181 stress 0.2310107 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0605918  max resid 0.150254 
    ## Run 182 stress 0.2324242 
    ## Run 183 stress 0.2323335 
    ## Run 184 stress 0.2403291 
    ## Run 185 stress 0.2335177 
    ## Run 186 stress 0.2410367 
    ## Run 187 stress 0.2362571 
    ## Run 188 stress 0.2342191 
    ## Run 189 stress 0.2382329 
    ## Run 190 stress 0.2354794 
    ## Run 191 stress 0.2385309 
    ## Run 192 stress 0.2384118 
    ## Run 193 stress 0.2350385 
    ## Run 194 stress 0.2333429 
    ## Run 195 stress 0.2383664 
    ## Run 196 stress 0.2364498 
    ## Run 197 stress 0.2371115 
    ## Run 198 stress 0.2350253 
    ## Run 199 stress 0.2337969 
    ## Run 200 stress 0.2344528 
    ## Run 201 stress 0.2377602 
    ## Run 202 stress 0.2310107 
    ## ... New best solution
    ## ... Procrustes: rmse 8.133094e-05  max resid 0.0002660058 
    ## ... Similar to previous best
    ## *** Solution reached

``` r
#make data frame of NMDS scores
nmds.silo3and9_log_bray.scores <- cbind(silo3and9$day, silo3and9$temp,data.frame(scores(nmds.silo3and9_log_bray)))
colnames(nmds.silo3and9_log_bray.scores)[1:2] <- c("day","temp")
ggplot(nmds.silo3and9_log_bray.scores, aes(NMDS1, NMDS2)) + geom_point(aes(col = day, shape = temp)) + theme_bw() + ggtitle("bray curtis NMDS of log ADJNSAF values with zeros replaced with 0.1")
```

![](ClusteringTechnicalReplicates_files/figure-markdown_github/unnamed-chunk-12-1.png)
