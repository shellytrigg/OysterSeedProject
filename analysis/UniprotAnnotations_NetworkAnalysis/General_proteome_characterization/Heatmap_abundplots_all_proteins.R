library(heatmap3)
library(colorRamps)
library(plotrix)
library(ggplot2)

#heatmap of all proteins

NSAF_for_big_heatmap <- read.csv("~/Documents/GitHub/paper-OysterSeed-TimeXTemp/Data/silo3and9_NSAF_AVGs.csv", stringsAsFactors = FALSE)
NSAF_for_big_heatmap[NSAF_for_big_heatmap == 0] <- 0.1
rownames(NSAF_for_big_heatmap) <- NSAF_for_big_heatmap$proteinID
NSAF_for_big_heatmap <- NSAF_for_big_heatmap[grep("CHOYP", NSAF_for_big_heatmap$proteinID),]
NSAF_for_big_heatmap$proteinID <- NULL


#make data frames to store heatmap color values
#columns are days/temp; color 23C blue, 29C red, each day a different shade of gray
ColSideColors<-cbind(Temp=c(rep(c("steelblue2","red"),6)), Day=c(rep("#D9D9D9",2),rep("#BDBDBD",2),rep("#969696",2),rep("#737373",2),rep("#525252",2),rep("#252525",2)))

heatmap3(NSAF_for_big_heatmap, cexRow = 0.1, cexCol = 0.5,ColSideColors = ColSideColors,ColAxisColors=1)
heatmap3(NSAF_for_big_heatmap[,grep("23C", colnames(NSAF_for_big_heatmap))])


hm <- as.dist(1-cor(t(NSAF_for_big_heatmap), use="pa"))
hm2 <- hclust(hm, method = 'average')
plot(hm2)

#Identify clades in heatmap
# define some clusters

#make a list of proteins with cluster IDs based on where the tree is cut
mycl <- cutree(hm2, h=0.8)

#assign colors to cluster IDs
clusterCols <- colorRamps::primary.colors(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
heatmap3(NSAF_for_big_heatmap,cexRow = 0.5, cexCol = 0.5, RowSideColors=myClusterSideBar,ColSideColors = ColSideColors,RowAxisColors=1,ColAxisColors=1,hclustfun=hclust, distfun = function(x) as.dist(1 - cor(t(x), use = "pa")), method = "average", scale = "row", Colv=NA)

#get color name
clusterColor <- lapply(myClusterSideBar, color.id)
clusterColor <- lapply(clusterColor, `[[`,1)
clusterColor <- data.frame(unlist(clusterColor), stringsAsFactors = FALSE)

foo <- cbind(data.frame(mycl), clusterColor)
#add protein column to data frame for merging
foo$protein_ID <- rownames(foo)

#add protein column back to NSAF data frame for merging
NSAF_for_big_heatmap_proteinID <- NSAF_for_big_heatmap
NSAF_for_big_heatmap_proteinID$protein_ID <- rownames(NSAF_for_big_heatmap_proteinID)

clade.prot <- merge(NSAF_for_big_heatmap_proteinID, foo, by = "protein_ID")

colnames(clade.prot)[15] <- "ClusterColor"
colnames(clade.prot)[14] <- "ClusterID"



#normalize within terms (autoscale by row)
#calculate row means and SD
clade.prot$mean <- apply(clade.prot[,2:13], 1, mean)
clade.prot$sd <- apply(clade.prot[,2:13],1, sd)

#sweep second argument is margin which corresponds to 1 (rows) or 2 (columns); this is for autoscaling by row to make colors relative to one another like the heatmap does
#subtract the row mean from each value
clade.prot_norm <- sweep(clade.prot[,2:13],1,clade.prot$mean)
#divide each value by the row standard deviation
clade.prot_norm <- sweep(clade.prot_norm[,1:12],1,clade.prot$sd, FUN = "/")

clade.prot_norm <- cbind(clade.prot[,c("protein_ID", "ClusterID", "ClusterColor")], clade.prot_norm)


STACKED_NSAF <- tidyr::gather(clade.prot_norm, daytemp, NSAF, 4:15)
STACKED_NSAF$day <- gsub(".*Day*|_.*","",STACKED_NSAF$daytemp)
STACKED_NSAF$temp <- gsub(".*_","",STACKED_NSAF$daytemp)
str(STACKED_NSAF)
STACKED_NSAF$day <- as.integer(STACKED_NSAF$day)

ggplot(STACKED_NSAF,aes(day,NSAF)) + geom_line(aes(group = protein_ID,alpha = 0.5)) + facet_wrap(~ClusterColor)  + scale_x_continuous(breaks = c(3,5,7,9,11,13), labels = c(3,5,7,9,11,13))


ggplot(STACKED_NSAF[which(STACKED_NSAF$ClusterColor=="gray50"),],aes(day,NSAF)) + geom_line(aes(group = protein_ID,color = factor(temp)), alpha = 0.02)  + geom_boxplot(aes(group = day)) + scale_x_continuous(breaks = c(3,5,7,9,11,13), labels = c(3,5,7,9,11,13))
