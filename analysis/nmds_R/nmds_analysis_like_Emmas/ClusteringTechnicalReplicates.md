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

![](ClusteringTechnicalReplicates_files/figure-markdown_github/unnamed-chunk-8-1.png)

try PCA on log transformed values ![](ClusteringTechnicalReplicates_files/figure-markdown_github/unnamed-chunk-9-1.png)

Make MDS dissimilarity matrix

``` r
nmds.silo3and9 <- metaMDS(silo3and9_nozerovar, distance = 'euclidean', k = 2, trymax = 3000, autotransform = FALSE)
```

    ## Run 0 stress 0.1649526 
    ## Run 1 stress 0.1832198 
    ## Run 2 stress 0.1832206 
    ## Run 3 stress 0.1627165 
    ## ... New best solution
    ## ... Procrustes: rmse 0.03376683  max resid 0.1320568 
    ## Run 4 stress 0.1649515 
    ## Run 5 stress 0.1832208 
    ## Run 6 stress 0.259257 
    ## Run 7 stress 0.1649516 
    ## Run 8 stress 0.1741252 
    ## Run 9 stress 0.1763797 
    ## Run 10 stress 0.1791324 
    ## Run 11 stress 0.1631588 
    ## ... Procrustes: rmse 0.01407279  max resid 0.05865797 
    ## Run 12 stress 0.1741334 
    ## Run 13 stress 0.2030333 
    ## Run 14 stress 0.1865933 
    ## Run 15 stress 0.1627149 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0009731256  max resid 0.003952801 
    ## ... Similar to previous best
    ## Run 16 stress 0.1753259 
    ## Run 17 stress 0.1631588 
    ## ... Procrustes: rmse 0.01483715  max resid 0.06111455 
    ## Run 18 stress 0.2311457 
    ## Run 19 stress 0.1832194 
    ## Run 20 stress 0.1725524 
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
    ## Run 1 stress 0.1229458 
    ## Run 2 stress 0.1263266 
    ## Run 3 stress 0.1229455 
    ## Run 4 stress 0.1260396 
    ## Run 5 stress 0.1124374 
    ## ... Procrustes: rmse 0.004872418  max resid 0.01682082 
    ## Run 6 stress 0.1229456 
    ## Run 7 stress 0.1374773 
    ## Run 8 stress 0.1362539 
    ## Run 9 stress 0.1124366 
    ## ... Procrustes: rmse 0.004796803  max resid 0.01675233 
    ## Run 10 stress 0.1271869 
    ## Run 11 stress 0.1259512 
    ## Run 12 stress 0.1123882 
    ## ... Procrustes: rmse 5.246288e-05  max resid 0.0001863706 
    ## ... Similar to previous best
    ## Run 13 stress 0.1124372 
    ## ... Procrustes: rmse 0.004860352  max resid 0.01681677 
    ## Run 14 stress 0.1123881 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002077532  max resid 0.0005166083 
    ## ... Similar to previous best
    ## Run 15 stress 0.1156691 
    ## Run 16 stress 0.1271873 
    ## Run 17 stress 0.1124365 
    ## ... Procrustes: rmse 0.004689534  max resid 0.01651689 
    ## Run 18 stress 0.1328571 
    ## Run 19 stress 0.1271869 
    ## Run 20 stress 0.1285323 
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
    ## Run 1 stress 0.2313428 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1458068  max resid 0.3633226 
    ## Run 2 stress 0.2364937 
    ## Run 3 stress 0.2402479 
    ## Run 4 stress 0.2363488 
    ## Run 5 stress 0.2370461 
    ## Run 6 stress 0.2340389 
    ## Run 7 stress 0.2367961 
    ## Run 8 stress 0.237666 
    ## Run 9 stress 0.2382266 
    ## Run 10 stress 0.2366783 
    ## Run 11 stress 0.2364045 
    ## Run 12 stress 0.2360764 
    ## Run 13 stress 0.2329337 
    ## Run 14 stress 0.2374006 
    ## Run 15 stress 0.2369469 
    ## Run 16 stress 0.234607 
    ## Run 17 stress 0.3724561 
    ## Run 18 stress 0.2314699 
    ## ... Procrustes: rmse 0.1020664  max resid 0.3460354 
    ## Run 19 stress 0.2370522 
    ## Run 20 stress 0.2394511 
    ## Run 21 stress 0.2390448 
    ## Run 22 stress 0.2465366 
    ## Run 23 stress 0.2351649 
    ## Run 24 stress 0.2390742 
    ## Run 25 stress 0.2368288 
    ## Run 26 stress 0.2411103 
    ## Run 27 stress 0.2360681 
    ## Run 28 stress 0.2377828 
    ## Run 29 stress 0.2346568 
    ## Run 30 stress 0.2397664 
    ## Run 31 stress 0.2359457 
    ## Run 32 stress 0.2359142 
    ## Run 33 stress 0.240143 
    ## Run 34 stress 0.2361214 
    ## Run 35 stress 0.2329272 
    ## Run 36 stress 0.2333391 
    ## Run 37 stress 0.2468811 
    ## Run 38 stress 0.2388654 
    ## Run 39 stress 0.2359586 
    ## Run 40 stress 0.2370228 
    ## Run 41 stress 0.2378328 
    ## Run 42 stress 0.2383541 
    ## Run 43 stress 0.2388285 
    ## Run 44 stress 0.2418964 
    ## Run 45 stress 0.2358165 
    ## Run 46 stress 0.2321129 
    ## Run 47 stress 0.2361022 
    ## Run 48 stress 0.2389842 
    ## Run 49 stress 0.2327497 
    ## Run 50 stress 0.2392564 
    ## Run 51 stress 0.2381602 
    ## Run 52 stress 0.2379955 
    ## Run 53 stress 0.2390927 
    ## Run 54 stress 0.2362435 
    ## Run 55 stress 0.2316022 
    ## ... Procrustes: rmse 0.08588781  max resid 0.3694785 
    ## Run 56 stress 0.2402256 
    ## Run 57 stress 0.2347727 
    ## Run 58 stress 0.2398042 
    ## Run 59 stress 0.234867 
    ## Run 60 stress 0.2425702 
    ## Run 61 stress 0.2365868 
    ## Run 62 stress 0.2377365 
    ## Run 63 stress 0.2377832 
    ## Run 64 stress 0.2343799 
    ## Run 65 stress 0.2333683 
    ## Run 66 stress 0.2348274 
    ## Run 67 stress 0.2360973 
    ## Run 68 stress 0.2363675 
    ## Run 69 stress 0.2434741 
    ## Run 70 stress 0.2385887 
    ## Run 71 stress 0.2372232 
    ## Run 72 stress 0.2342762 
    ## Run 73 stress 0.2315368 
    ## ... Procrustes: rmse 0.07020007  max resid 0.3173675 
    ## Run 74 stress 0.2388284 
    ## Run 75 stress 0.2369661 
    ## Run 76 stress 0.2375142 
    ## Run 77 stress 0.2372546 
    ## Run 78 stress 0.2381769 
    ## Run 79 stress 0.2353885 
    ## Run 80 stress 0.238423 
    ## Run 81 stress 0.2375822 
    ## Run 82 stress 0.2429242 
    ## Run 83 stress 0.2410518 
    ## Run 84 stress 0.236439 
    ## Run 85 stress 0.235861 
    ## Run 86 stress 0.2321361 
    ## Run 87 stress 0.2399469 
    ## Run 88 stress 0.2388185 
    ## Run 89 stress 0.233536 
    ## Run 90 stress 0.2378115 
    ## Run 91 stress 0.234009 
    ## Run 92 stress 0.2320098 
    ## Run 93 stress 0.2365909 
    ## Run 94 stress 0.2347061 
    ## Run 95 stress 0.2387007 
    ## Run 96 stress 0.2398493 
    ## Run 97 stress 0.2347059 
    ## Run 98 stress 0.2333436 
    ## Run 99 stress 0.2370854 
    ## Run 100 stress 0.2382268 
    ## Run 101 stress 0.2352144 
    ## Run 102 stress 0.240805 
    ## Run 103 stress 0.2394035 
    ## Run 104 stress 0.238735 
    ## Run 105 stress 0.2343232 
    ## Run 106 stress 0.2333061 
    ## Run 107 stress 0.2356335 
    ## Run 108 stress 0.2386612 
    ## Run 109 stress 0.2351218 
    ## Run 110 stress 0.2531515 
    ## Run 111 stress 0.2348979 
    ## Run 112 stress 0.2356158 
    ## Run 113 stress 0.2369522 
    ## Run 114 stress 0.2346104 
    ## Run 115 stress 0.233454 
    ## Run 116 stress 0.2370171 
    ## Run 117 stress 0.2380408 
    ## Run 118 stress 0.2355262 
    ## Run 119 stress 0.2372064 
    ## Run 120 stress 0.2399267 
    ## Run 121 stress 0.2337609 
    ## Run 122 stress 0.239003 
    ## Run 123 stress 0.2397839 
    ## Run 124 stress 0.2368674 
    ## Run 125 stress 0.2335189 
    ## Run 126 stress 0.2336484 
    ## Run 127 stress 0.2386741 
    ## Run 128 stress 0.2385608 
    ## Run 129 stress 0.2382124 
    ## Run 130 stress 0.2379988 
    ## Run 131 stress 0.2356661 
    ## Run 132 stress 0.2361931 
    ## Run 133 stress 0.2413376 
    ## Run 134 stress 0.2379521 
    ## Run 135 stress 0.2416726 
    ## Run 136 stress 0.2331773 
    ## Run 137 stress 0.2365435 
    ## Run 138 stress 0.2449909 
    ## Run 139 stress 0.2325837 
    ## Run 140 stress 0.2328663 
    ## Run 141 stress 0.2368792 
    ## Run 142 stress 0.2356577 
    ## Run 143 stress 0.2394755 
    ## Run 144 stress 0.2389805 
    ## Run 145 stress 0.238473 
    ## Run 146 stress 0.2370407 
    ## Run 147 stress 0.2386311 
    ## Run 148 stress 0.2382614 
    ## Run 149 stress 0.236868 
    ## Run 150 stress 0.239128 
    ## Run 151 stress 0.2400806 
    ## Run 152 stress 0.2388831 
    ## Run 153 stress 0.234721 
    ## Run 154 stress 0.2401859 
    ## Run 155 stress 0.2323956 
    ## Run 156 stress 0.2316004 
    ## ... Procrustes: rmse 0.08845173  max resid 0.3586654 
    ## Run 157 stress 0.2416707 
    ## Run 158 stress 0.2385271 
    ## Run 159 stress 0.2336002 
    ## Run 160 stress 0.2376864 
    ## Run 161 stress 0.2333428 
    ## Run 162 stress 0.2379842 
    ## Run 163 stress 0.2365719 
    ## Run 164 stress 0.2373959 
    ## Run 165 stress 0.240606 
    ## Run 166 stress 0.236989 
    ## Run 167 stress 0.237125 
    ## Run 168 stress 0.2367116 
    ## Run 169 stress 0.2366561 
    ## Run 170 stress 0.2334105 
    ## Run 171 stress 0.2365188 
    ## Run 172 stress 0.2376501 
    ## Run 173 stress 0.2406297 
    ## Run 174 stress 0.2391441 
    ## Run 175 stress 0.2378623 
    ## Run 176 stress 0.2432955 
    ## Run 177 stress 0.2351587 
    ## Run 178 stress 0.2395497 
    ## Run 179 stress 0.2396217 
    ## Run 180 stress 0.2362793 
    ## Run 181 stress 0.2375415 
    ## Run 182 stress 0.2398816 
    ## Run 183 stress 0.2356711 
    ## Run 184 stress 0.2377949 
    ## Run 185 stress 0.2360529 
    ## Run 186 stress 0.2404365 
    ## Run 187 stress 0.2373259 
    ## Run 188 stress 0.2363773 
    ## Run 189 stress 0.2374379 
    ## Run 190 stress 0.2358861 
    ## Run 191 stress 0.2395051 
    ## Run 192 stress 0.2382456 
    ## Run 193 stress 0.2368866 
    ## Run 194 stress 0.2374586 
    ## Run 195 stress 0.2360634 
    ## Run 196 stress 0.2386988 
    ## Run 197 stress 0.2341081 
    ## Run 198 stress 0.2359859 
    ## Run 199 stress 0.233491 
    ## Run 200 stress 0.2394199 
    ## Run 201 stress 0.2378654 
    ## Run 202 stress 0.2716536 
    ## Run 203 stress 0.2382932 
    ## Run 204 stress 0.2353914 
    ## Run 205 stress 0.2401152 
    ## Run 206 stress 0.238925 
    ## Run 207 stress 0.2392378 
    ## Run 208 stress 0.2414799 
    ## Run 209 stress 0.234033 
    ## Run 210 stress 0.2361485 
    ## Run 211 stress 0.2420844 
    ## Run 212 stress 0.237379 
    ## Run 213 stress 0.2379621 
    ## Run 214 stress 0.2376154 
    ## Run 215 stress 0.2335991 
    ## Run 216 stress 0.2380525 
    ## Run 217 stress 0.233476 
    ## Run 218 stress 0.2336445 
    ## Run 219 stress 0.2392145 
    ## Run 220 stress 0.2385872 
    ## Run 221 stress 0.2366633 
    ## Run 222 stress 0.2341167 
    ## Run 223 stress 0.2394562 
    ## Run 224 stress 0.2364592 
    ## Run 225 stress 0.2379933 
    ## Run 226 stress 0.2335179 
    ## Run 227 stress 0.2339142 
    ## Run 228 stress 0.2391985 
    ## Run 229 stress 0.2365713 
    ## Run 230 stress 0.2375054 
    ## Run 231 stress 0.2350331 
    ## Run 232 stress 0.2322168 
    ## Run 233 stress 0.234208 
    ## Run 234 stress 0.2406371 
    ## Run 235 stress 0.2367387 
    ## Run 236 stress 0.2390895 
    ## Run 237 stress 0.239314 
    ## Run 238 stress 0.239328 
    ## Run 239 stress 0.2375479 
    ## Run 240 stress 0.2382244 
    ## Run 241 stress 0.2351707 
    ## Run 242 stress 0.2356261 
    ## Run 243 stress 0.2381296 
    ## Run 244 stress 0.2376141 
    ## Run 245 stress 0.2360564 
    ## Run 246 stress 0.238931 
    ## Run 247 stress 0.2374548 
    ## Run 248 stress 0.2350166 
    ## Run 249 stress 0.2398586 
    ## Run 250 stress 0.2370063 
    ## Run 251 stress 0.2386392 
    ## Run 252 stress 0.2367985 
    ## Run 253 stress 0.2368196 
    ## Run 254 stress 0.2341811 
    ## Run 255 stress 0.2375984 
    ## Run 256 stress 0.236558 
    ## Run 257 stress 0.235367 
    ## Run 258 stress 0.2342758 
    ## Run 259 stress 0.2377774 
    ## Run 260 stress 0.2404078 
    ## Run 261 stress 0.2477551 
    ## Run 262 stress 0.2409076 
    ## Run 263 stress 0.2378621 
    ## Run 264 stress 0.2383172 
    ## Run 265 stress 0.2310107 
    ## ... New best solution
    ## ... Procrustes: rmse 0.06663266  max resid 0.2443927 
    ## Run 266 stress 0.234445 
    ## Run 267 stress 0.2375192 
    ## Run 268 stress 0.2398231 
    ## Run 269 stress 0.2427041 
    ## Run 270 stress 0.2384114 
    ## Run 271 stress 0.2376482 
    ## Run 272 stress 0.234759 
    ## Run 273 stress 0.2436208 
    ## Run 274 stress 0.2374954 
    ## Run 275 stress 0.2423879 
    ## Run 276 stress 0.2374785 
    ## Run 277 stress 0.2386467 
    ## Run 278 stress 0.2370219 
    ## Run 279 stress 0.2367836 
    ## Run 280 stress 0.2328143 
    ## Run 281 stress 0.2350397 
    ## Run 282 stress 0.2367797 
    ## Run 283 stress 0.2386374 
    ## Run 284 stress 0.2373857 
    ## Run 285 stress 0.2375099 
    ## Run 286 stress 0.2384361 
    ## Run 287 stress 0.2421663 
    ## Run 288 stress 0.2332282 
    ## Run 289 stress 0.2380342 
    ## Run 290 stress 0.2372818 
    ## Run 291 stress 0.2364577 
    ## Run 292 stress 0.2403051 
    ## Run 293 stress 0.2372664 
    ## Run 294 stress 0.2319103 
    ## Run 295 stress 0.2387637 
    ## Run 296 stress 0.237508 
    ## Run 297 stress 0.2383274 
    ## Run 298 stress 0.2391097 
    ## Run 299 stress 0.2358838 
    ## Run 300 stress 0.2317076 
    ## Run 301 stress 0.2396182 
    ## Run 302 stress 0.2425557 
    ## Run 303 stress 0.2369573 
    ## Run 304 stress 0.2372257 
    ## Run 305 stress 0.2358497 
    ## Run 306 stress 0.239801 
    ## Run 307 stress 0.2373582 
    ## Run 308 stress 0.2407613 
    ## Run 309 stress 0.2344532 
    ## Run 310 stress 0.2369107 
    ## Run 311 stress 0.234225 
    ## Run 312 stress 0.2378409 
    ## Run 313 stress 0.2330545 
    ## Run 314 stress 0.242632 
    ## Run 315 stress 0.2372947 
    ## Run 316 stress 0.244848 
    ## Run 317 stress 0.2387162 
    ## Run 318 stress 0.232726 
    ## Run 319 stress 0.2373622 
    ## Run 320 stress 0.2362113 
    ## Run 321 stress 0.2396332 
    ## Run 322 stress 0.2413679 
    ## Run 323 stress 0.2383175 
    ## Run 324 stress 0.2380269 
    ## Run 325 stress 0.2351827 
    ## Run 326 stress 0.2404292 
    ## Run 327 stress 0.2407722 
    ## Run 328 stress 0.2302106 
    ## ... New best solution
    ## ... Procrustes: rmse 0.03240022  max resid 0.1161295 
    ## Run 329 stress 0.2404505 
    ## Run 330 stress 0.2380866 
    ## Run 331 stress 0.2396991 
    ## Run 332 stress 0.2386208 
    ## Run 333 stress 0.2392623 
    ## Run 334 stress 0.2336446 
    ## Run 335 stress 0.2375542 
    ## Run 336 stress 0.2425594 
    ## Run 337 stress 0.2335924 
    ## Run 338 stress 0.2328226 
    ## Run 339 stress 0.2389859 
    ## Run 340 stress 0.2371751 
    ## Run 341 stress 0.2371668 
    ## Run 342 stress 0.2387302 
    ## Run 343 stress 0.2379095 
    ## Run 344 stress 0.2364654 
    ## Run 345 stress 0.236966 
    ## Run 346 stress 0.2383532 
    ## Run 347 stress 0.2375555 
    ## Run 348 stress 0.2350945 
    ## Run 349 stress 0.239695 
    ## Run 350 stress 0.2370225 
    ## Run 351 stress 0.2405434 
    ## Run 352 stress 0.2369617 
    ## Run 353 stress 0.2380981 
    ## Run 354 stress 0.2368904 
    ## Run 355 stress 0.2422526 
    ## Run 356 stress 0.237262 
    ## Run 357 stress 0.2401652 
    ## Run 358 stress 0.243921 
    ## Run 359 stress 0.2325327 
    ## Run 360 stress 0.2393137 
    ## Run 361 stress 0.2369452 
    ## Run 362 stress 0.2361563 
    ## Run 363 stress 0.2328554 
    ## Run 364 stress 0.2336422 
    ## Run 365 stress 0.2378813 
    ## Run 366 stress 0.235963 
    ## Run 367 stress 0.2356567 
    ## Run 368 stress 0.2393929 
    ## Run 369 stress 0.2387782 
    ## Run 370 stress 0.2395064 
    ## Run 371 stress 0.2372191 
    ## Run 372 stress 0.237122 
    ## Run 373 stress 0.2391476 
    ## Run 374 stress 0.2398654 
    ## Run 375 stress 0.237975 
    ## Run 376 stress 0.2320097 
    ## Run 377 stress 0.232573 
    ## Run 378 stress 0.2379365 
    ## Run 379 stress 0.2369432 
    ## Run 380 stress 0.2348825 
    ## Run 381 stress 0.2336021 
    ## Run 382 stress 0.2376991 
    ## Run 383 stress 0.2368117 
    ## Run 384 stress 0.2416757 
    ## Run 385 stress 0.2344346 
    ## Run 386 stress 0.2377182 
    ## Run 387 stress 0.2385878 
    ## Run 388 stress 0.2368033 
    ## Run 389 stress 0.2387124 
    ## Run 390 stress 0.2332202 
    ## Run 391 stress 0.2404613 
    ## Run 392 stress 0.2417184 
    ## Run 393 stress 0.2357073 
    ## Run 394 stress 0.2430561 
    ## Run 395 stress 0.2381867 
    ## Run 396 stress 0.2327261 
    ## Run 397 stress 0.2375731 
    ## Run 398 stress 0.2364723 
    ## Run 399 stress 0.2425346 
    ## Run 400 stress 0.238698 
    ## Run 401 stress 0.2377828 
    ## Run 402 stress 0.2357314 
    ## Run 403 stress 0.2425869 
    ## Run 404 stress 0.2359495 
    ## Run 405 stress 0.2376058 
    ## Run 406 stress 0.2361988 
    ## Run 407 stress 0.2352968 
    ## Run 408 stress 0.2384272 
    ## Run 409 stress 0.2342138 
    ## Run 410 stress 0.2390252 
    ## Run 411 stress 0.2359703 
    ## Run 412 stress 0.2367492 
    ## Run 413 stress 0.2409673 
    ## Run 414 stress 0.2334393 
    ## Run 415 stress 0.2369825 
    ## Run 416 stress 0.2384647 
    ## Run 417 stress 0.2385943 
    ## Run 418 stress 0.2396511 
    ## Run 419 stress 0.2334394 
    ## Run 420 stress 0.2314091 
    ## Run 421 stress 0.2433825 
    ## Run 422 stress 0.2338505 
    ## Run 423 stress 0.2387272 
    ## Run 424 stress 0.2338772 
    ## Run 425 stress 0.2369498 
    ## Run 426 stress 0.2368214 
    ## Run 427 stress 0.2313424 
    ## Run 428 stress 0.2371451 
    ## Run 429 stress 0.2314698 
    ## Run 430 stress 0.2384673 
    ## Run 431 stress 0.2371123 
    ## Run 432 stress 0.2336002 
    ## Run 433 stress 0.2362233 
    ## Run 434 stress 0.2388123 
    ## Run 435 stress 0.2366271 
    ## Run 436 stress 0.2419146 
    ## Run 437 stress 0.2324751 
    ## Run 438 stress 0.2398718 
    ## Run 439 stress 0.2350743 
    ## Run 440 stress 0.2381396 
    ## Run 441 stress 0.2401936 
    ## Run 442 stress 0.237855 
    ## Run 443 stress 0.2376352 
    ## Run 444 stress 0.2327924 
    ## Run 445 stress 0.2365432 
    ## Run 446 stress 0.2377966 
    ## Run 447 stress 0.2389353 
    ## Run 448 stress 0.238413 
    ## Run 449 stress 0.2371314 
    ## Run 450 stress 0.2337565 
    ## Run 451 stress 0.2320839 
    ## Run 452 stress 0.2381213 
    ## Run 453 stress 0.2302105 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001315476  max resid 0.0003277053 
    ## ... Similar to previous best
    ## *** Solution reached

``` r
#make data frame of NMDS scores
nmds.silo3and9_log_bray.scores <- cbind(silo3and9$day, silo3and9$temp,data.frame(scores(nmds.silo3and9_log_bray)))
colnames(nmds.silo3and9_log_bray.scores)[1:2] <- c("day","temp")
ggplot(nmds.silo3and9_log_bray.scores, aes(NMDS1, NMDS2)) + geom_point(aes(col = day, shape = temp)) + theme_bw() + ggtitle("bray curtis NMDS of log ADJNSAF values with zeros replaced with 0.1")
```

![](ClusteringTechnicalReplicates_files/figure-markdown_github/unnamed-chunk-12-1.png)
