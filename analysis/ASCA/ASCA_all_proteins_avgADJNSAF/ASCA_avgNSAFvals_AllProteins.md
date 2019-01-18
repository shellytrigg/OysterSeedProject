ASCA on average NSAF values of all proteins
================
Shelly Trigg
1/17/2019

load libraries

    ## Warning: package 'dplyr' was built under R version 3.4.4

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## Warning: package 'tidyr' was built under R version 3.4.4

    ## Loading required package: MASS

    ## Warning: package 'MASS' was built under R version 3.4.4

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## Loading required package: abind

    ## Loading required package: pls

    ## Warning: package 'pls' was built under R version 3.4.4

    ## 
    ## Attaching package: 'pls'

    ## The following object is masked from 'package:stats':
    ## 
    ##     loadings

    ## Warning: package 'ggplot2' was built under R version 3.4.4

load data

``` r
#NSAF data from filtered proteins
data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/silo3and9_nozerovals_AVGs.csv", stringsAsFactors = FALSE)
```

Load data with annotations so protein IDs can be substituted with shorter names

subset the data for only CHOYP proteins; this removes contaminant proteins

perform ASCA

``` r
#create matrix to pass to ASCA command, excluding the silo and time info
ASCA_X <- as.matrix(choyp_data[,-c(1:2)])
#create matrix to pass to ASCA command with only the silo and time info
ASCA_F <- as.matrix(choyp_data[,c(1:2)])
#perform ASCA
ASCA <- ASCA.Calculate(ASCA_X, ASCA_F, equation.elements = "1,2,12", scaling = FALSE)
```

    ## Variance explained per principal component (if >1%):
    ## Whole data set   PC1: 26.97%   PC2: 21.08%   PC3: 10.76%   PC4: 9.43%    PC5: 7.07%    PC6: 6.64%    PC7: 4.77%    PC8: 3.83%    PC9: 3.17%    PC10: 2.54%   
    ## Factor 1         PC1: 39.62%   PC2: 25.87%   PC3: 13.38%   PC4: 9.34%    PC5: 8.24%    PC6: 3.55%    PC7:  NA%     PC8:  NA%     PC9:  NA%     PC10:  NA%    
    ## Factor 2         PC1: 70.07%   PC2: 29.93%   PC3:  NA%     PC4:  NA%     PC5:  NA%     PC6:  NA%     PC7:  NA%     PC8:  NA%     PC9:  NA%     PC10:  NA%    
    ## Interaction 12   PC1: 40.79%   PC2: 17.69%   PC3: 13.61%   PC4: 11.72%   PC5: 9.88%    PC6: 6.30%    PC7:  NA%     PC8:  NA%     PC9:  NA%     PC10:  NA%    
    ## 
    ## Percentage each effect contributes to the total sum of squares:
    ## Overall means    91.22%
    ## Factor 1         5.60%
    ## Factor 2         2.04%
    ## Interaction 12   3.96%
    ## Residuals        0.00%
    ## 
    ## Percentage each effect contributes to the sum of squares of the centered data:
    ## Factor 1         63.77%
    ## Factor 2         23.23%
    ## Interaction 12   45.08%
    ## Residuals        0.00%

Here is a summary of the ASCA results (e.g. variance explained by different factors; factor 1= time (days), factor 2 = temperature, interaction = interaction of time and temperature)

``` r
#print the ASCA summary
ASCA.GetSummary(ASCA)
```

    ## Variance explained per principal component (if >1%):
    ## Whole data set   PC1: 26.97%   PC2: 21.08%   PC3: 10.76%   PC4: 9.43%    PC5: 7.07%    PC6: 6.64%    PC7: 4.77%    PC8: 3.83%    PC9: 3.17%    PC10: 2.54%   
    ## Factor 1         PC1: 39.62%   PC2: 25.87%   PC3: 13.38%   PC4: 9.34%    PC5: 8.24%    PC6: 3.55%    PC7:  NA%     PC8:  NA%     PC9:  NA%     PC10:  NA%    
    ## Factor 2         PC1: 70.07%   PC2: 29.93%   PC3:  NA%     PC4:  NA%     PC5:  NA%     PC6:  NA%     PC7:  NA%     PC8:  NA%     PC9:  NA%     PC10:  NA%    
    ## Interaction 12   PC1: 40.79%   PC2: 17.69%   PC3: 13.61%   PC4: 11.72%   PC5: 9.88%    PC6: 6.30%    PC7:  NA%     PC8:  NA%     PC9:  NA%     PC10:  NA%    
    ## 
    ## Percentage each effect contributes to the total sum of squares:
    ## Overall means    91.22%
    ## Factor 1         5.60%
    ## Factor 2         2.04%
    ## Interaction 12   3.96%
    ## Residuals        0.00%
    ## 
    ## Percentage each effect contributes to the sum of squares of the centered data:
    ## Factor 1         63.77%
    ## Factor 2         23.23%
    ## Interaction 12   45.08%
    ## Residuals        0.00%

    ## $summary.pca
    ##            PC1       PC2       PC3        PC4        PC5        PC6
    ## data 0.2697456 0.2108149 0.1076162 0.09428784 0.07073832 0.06639430
    ## 1    0.3962349 0.2586869 0.1338009 0.09343167 0.08235720 0.03548848
    ## 2    0.7007309 0.2992691        NA         NA         NA         NA
    ## 12   0.4079499 0.1768866 0.1361433 0.11717182 0.09883591 0.06301248
    ##            PC7        PC8        PC9       PC10
    ## data 0.0476876 0.03833476 0.03173637 0.02537494
    ## 1           NA         NA         NA         NA
    ## 2           NA         NA         NA         NA
    ## 12          NA         NA         NA         NA
    ## 
    ## $summary.ssq
    ##                     Overall means          1          2         12
    ## Contribution to ssq     0.9122141 0.05598042 0.02039528 0.03957785
    ##                        Residuals
    ## Contribution to ssq 2.301613e-34

### Plot PCAs from ASCA

This first plot is the time (days) effect PCA

``` r
#plot PCA for factor 1, which is time in this case
ASCA.PlotScoresPerLevel(ASCA, ee = "1", pcs = "1,2")
```

![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_PCA_timeEffect_plot-1.png)

This next plot is the temperature effect PCA

``` r
#plot PCA for factor 2, which is temperature in this case
ASCA.PlotScoresPerLevel(ASCA, ee = "2", pcs = "1,2")
```

![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_PCA_tempEffect_plot-1.png)

This next plot is the time x temp interaction effect PCA

``` r
#plot PCA for factor interaction, which is time x temp in this case
timextemp_PC12 <- data.frame(ASCA$`12`$svd$t[,c(1,2)])
timextemp_PC12 <- cbind(data.frame(ASCA$`12`$level.combinations$row.patterns), timextemp_PC12)
colnames(timextemp_PC12)<- c("day","temp","PC1","PC2")
timextemp_PC12$day <- as.character(timextemp_PC12$day)
timextemp_PC12$temp <- as.character(timextemp_PC12$temp)
ggplot(timextemp_PC12, aes(PC1, PC2)) + geom_point(aes(col = temp, shape = day, size = 3)) + theme_bw() + ggtitle("PC1 vs PC2 for time x temperature interaction effect") + theme(plot.title = element_text(face = "bold")) + xlab(paste("PC1"," (",formatC(ASCA$`12`$svd$var.explained[1] * 100,digits=2,format="f"),"%)", sep = "")) + ylab(paste("PC2"," (",formatC(ASCA$`12`$svd$var.explained[2] * 100,digits=2,format="f"),"%)", sep = ""))
```

    ## Warning: The shape palette can deal with a maximum of 6 discrete values
    ## because more than 6 becomes difficult to discriminate; you have 7.
    ## Consider specifying shapes manually if you must have them.

    ## Warning: Removed 2 rows containing missing values (geom_point).

![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_PCA_timeXtempEffect_plot-1.png)

### Analysis of proteins affected by temperature

plot PC2 loadings for temperature effect

``` r
#extract protein names from ASCA data; these will be combined with loadings values; the order is maintained
protnames <- data.frame(colnames(ASCA$data))
#combine protein names with ASCA loadings for PC2 for temperature, since this component showed the greatest separation between temperatures
d <- cbind(protnames, ASCA$`2`$svd$v[,2])
#rename the columns
colnames(d) <- c("protein", "PC2loadings")
#make a dataframe of proteins with PC loadings greater than zero for the loadings plot
d_great <- d[which(d$PC2loadings > 0),]
#make a dataframe of proteins with PC loadings less than zero for the loadings plot
d_less <- d[which(d$PC2loadings < 0),]
#plot PC1 loadings
plot(d_great[order(d_great$PC2loadings, decreasing = TRUE),2], xlab = "protein", ylab = "PC2 Loadings")
```

![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_PCA_tempEffect_PC2loadings-1.png)

``` r
plot(d_less[order(d_less$PC2loadings, decreasing = TRUE),2],xlab = "protein", ylab = "PC2 Loadings")
```

![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_PCA_tempEffect_PC2loadings-2.png)

To pull out proteins affected by temperature based on their influence in seperating treatment groups on PC2 of the temperature PCA, I picked an absolute value loadings threshold of 0.025. This means any protein that had a loadings value \> 0.025 or \< -0.025 was selected.

``` r
#make list of proteins with temperature PC1 loadings values >= 0.025
cutd <- d[which(abs(d$PC2loadings) >= 0.025),]
#number of proteins affected by temp
nrow(cutd)
```

    ## [1] 138

``` r
#make a list of cutoff proteins with normalized abundance
cut_data <-data[,which(colnames(data) %in% cutd$protein)]
rownames(cut_data) <- paste(data$day, data$temp, sep="_")
```

Plot a heatmap of proteins affected by temperature based on temperature effect PC2 loadings value cutoff

``` r
#reorder the cut data so that samples are grouped by temperature
cut_data_ord <- data[,c(1:3,which(colnames(data) %in% cutd$protein))]
cut_data_ord <- cut_data_ord[order(cut_data_ord$temp,cut_data_ord$day),]
rownames(cut_data_ord) <- paste(cut_data_ord$day, cut_data_ord$temp, sep="_")
cut_data_ord <- cut_data_ord[,-c(1:3)]
#replace CHOYP names with "entry names"
#transposed to eventually get row names as columns to match entry names to. 
#t.data.frame returns a matrix so I had convert back to data frame for all downstream rearranging because data frames are much easier to work with
cut_data_ord_t <- data.frame(t.data.frame(cut_data_ord))
#make row names which contain the protein IDs a column
cut_data_ord_t <- cbind(rownames(cut_data_ord_t), cut_data_ord_t)
#add the column name
colnames(cut_data_ord_t)[1] <- "Protein.ID"
#convert the Protein ID column to character
cut_data_ord_t[,1] <- as.character(cut_data_ord_t[,1])
#add the Entry names
cut_data_ord_t <- merge(data_w_uniprot[,c("Protein.ID","Entry.name.new")],cut_data_ord_t, by = "Protein.ID")
#make the entry names the row names
rownames(cut_data_ord_t) <- cut_data_ord_t[,"Entry.name.new"]
#remove the protein ID and entry name columns since the row names are added
cut_data_ord_t <- cut_data_ord_t[,-c(1:2)]
#this step removes the x from the column names while preserving the order of colnames
colnames(cut_data_ord_t) <- rownames(cut_data_ord)
#plot the heat map
heatmap3(as.matrix(cut_data_ord_t), Colv = NA, cexRow = 0.5)
```

![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_tempEffectPC2_cutoff0.025_heatmap_OrderedByTemp_ShortNames-1.png)

### Analysis of proteins affected by time (development associated proteins)

plot PC1 loadings for the time effect PCA ![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_PCA_timeEffect_PC1loadings-1.png)![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_PCA_timeEffect_PC1loadings-2.png)

plot PC2 loadings for the time effect PCA ![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_PCA_timeEffect_PC2loadings-1.png)![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_PCA_timeEffect_PC2loadings-2.png)

    ## [1] 230

Plot a heatmap of proteins affected by time based on time effect PC1 and PC2 loadings value cutoffs ![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_timeEffectPC1and2_cutoff0.025_heatmap_OrderedByTime-1.png)

### Analysis of proteins affected by the interaction of time and temp

plot PC1 loadings for the time x temperature interaction effect PCA ![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_PCA_timeXtempEffect_PC1loadings-1.png)![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_PCA_timeXtempEffect_PC1loadings-2.png)

plot PC2 loadings for the time x temperature interaction effect PCA ![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_PCA_timeXtempEffect_PC2loadings-1.png)![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_PCA_timeXtempEffect_PC2loadings-2.png)

    ## [1] 219

Plot a heatmap of proteins affected by the time x temperature interaction based on time x temperature effect PC1 and PC2 loadings value cutoffs ![](ASCA_avgNSAFvals_AllProteins_files/figure-markdown_github/avgNSAF_timeXtempEffectPC1and2_cutoff0.025_heatmap_OrderedByTime-1.png)
