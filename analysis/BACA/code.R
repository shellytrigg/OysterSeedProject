
#Downloading Bioconductor software (contains RDAVIDWebService, etc.)
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("RDAVIDWebService")

#documentation of package version
browseVignettes("RDAVIDWebService")

#Must have JAVA 8
library(rJava)

R CMD javareconf
library(RDAVIDWebService)
library(ggplot2)

install.packages("BACA")
library(BACA)

#Tutorial
data(gene.lists.ex)
str(gene.lists.ex)

result.kegg <- DAVIDsearch(gene.lists.ex, david.user = "krmitch7@uw.edu", 
             idType="ENTREZ_GENE_ID", annotation="KEGG_PATHWAY")
