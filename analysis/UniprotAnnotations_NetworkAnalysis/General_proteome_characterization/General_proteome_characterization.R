library(tidyr)
library(GSEABase)
library(plyr)
library(RColorBrewer)

#read in data
NSAF <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/silo3and9_nozerovals_AVGs.csv", stringsAsFactors = FALSE)
sr_lab_goslim <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/background/GOSlim_terms.csv", stringsAsFactors = FALSE)

#I had already changed the zeros to 0.1, so I will sub that back with zero for this analysis
NSAF[NSAF == 0.1] <- 0

# #order NSAF data by temp and remove day 0
# NSAF <- NSAF[which(NSAF$day != 0),]
# NSAF <- NSAF[order(NSAF$temp),]
# 
# NSAF_t <- data.frame(t(NSAF[,-c(1,2)]))
# colnames(NSAF_t) <- paste(NSAF$day, NSAF$temp, sep = "_")
# 

#separate 3 and 9

silo3 <- NSAF[grep("23", NSAF$temp),]
silo9 <- NSAF[grep("29", NSAF$temp),]

#reformat to long list
silo3_t <- data.frame(t(silo3[,-c(1,2)]))
colnames(silo3_t) <- paste0("D",silo3$day)

#remove proteins that were only in silo9
View(silo3_t[which(apply(silo3_t, 1, var) == 0),])
no_val_proteins <- rownames(silo3_t[which(apply(silo3_t, 1, var) == 0),])

View(silo3_t[-c(which(rownames(silo3_t) %in% no_val_proteins)),])
silo3_t_nozero <- silo3_t[-c(which(rownames(silo3_t) %in% no_val_proteins)),]
#confirm it worked
nrow(silo3_t)-nrow(silo3_t_nozero)
length(no_val_proteins)

#add protein column
silo3_t_nozero$protein_ID <- rownames(silo3_t_nozero)
silo3_t_nozero <- silo3_t_nozero[grep("CHOYP", silo3_t_nozero$protein_ID),]

#merge with uniprot
#uniprot <- read.csv("/Volumes/web/metacarcinus/Cgigas/all_giga-uniprot-blastP-out.nopipe.annotations.tab", sep ="\t", header = FALSE, stringsAsFactors = FALSE)
#combine uniprot and NSAF data
silo3_t_nozero_unip<- merge(silo3_t_nozero, uniprot, by = "protein_ID", all.x = TRUE)
#exclude proteins that didn't map to uniprot DB
silo3_t_nozero_unip_mapped <- silo3_t_nozero_unip[-grep("unmapped", silo3_t_nozero_unip$Entry),]
nrow(silo3_t_nozero_unip_mapped)
silo3_t_nozero_unip_mapped <- silo3_t_nozero_unip_mapped[which(silo3_t_nozero_unip_mapped$evalue <= 10^-10),]
nrow(silo3_t_nozero_unip_mapped)


##get go terms for each day
silo3_all_srlab_terms_per_day <- data.frame()
silo3_all_srlab_BPterms_per_day <- data.frame()

days <- colnames(silo3_t_nozero_unip_mapped[,2:7])
for(i in days){
  silo3_pro_GO <- silo3_t_nozero_unip_mapped[which(silo3_t_nozero_unip_mapped[,i] != 0),c("protein_ID","GO_IDs")]

  silo3_pro_GOid_term <- data.frame()
#loop through each line of data frame containing one column with "protein_ID" and one column with list of "GO_IDs" separated by ";"
  for (j in 1:nrow(silo3_pro_GO)){
    #create a row that lists all GO IDs associated with one protein listed across rows, where each different GO ID gets put in it's own column
    silo3_pro_GOid_term_row <- data.frame(t(data.frame(strsplit(as.character(silo3_pro_GO$GO_IDs[j]),'; ', fixed = TRUE))), stringsAsFactors = FALSE)
    #add each row created in the line above to the empty data frame
    #this will add NAs in columns where less terms exist; e.g. if a protein only has two terms and another has 10, there will be 8 NAs added to the row for the protein with 2 terms
    silo3_pro_GOid_term <- rbind.fill(silo3_pro_GOid_term,silo3_pro_GOid_term_row)
  }

  #add protein IDs back to GO IDs
  silo3_pro_GOid_term <- cbind(silo3_pro_GO[,"protein_ID"], silo3_pro_GOid_term)
  #this results in a table with protein ID listed in one column and each next column contains a GO ID that was listed with the protein in Uniprot DB.
  
  #reshape data so that all GO ID columns are gathered in one column called "GO" 
  STACKED_silo3_pro_GOid_term <- tidyr::gather(silo3_pro_GOid_term,"protein_ID","GO", 2:ncol(silo3_pro_GOid_term))
  #exlude middle column which just contains the string "protein_ID" in each row
  STACKED_silo3_pro_GOid_term <- STACKED_silo3_pro_GOid_term[,c(1,3)]
  #remove duplicate rows
  STACKED_silo3_pro_GOid_term <- unique(STACKED_silo3_pro_GOid_term)
  colnames(STACKED_silo3_pro_GOid_term)[1] <- "protein_ID"
  #remove any rows where GO column has NA value. 
  STACKED_silo3_pro_GOid_term <- STACKED_silo3_pro_GOid_term[which(!is.na(STACKED_silo3_pro_GOid_term$GO)),]
  #this resulting data frame has two columns "protein_ID" and "GO"
  
  ##get go slim terms for each day
  STACKED_silo3_pro_GOid_term_sr_labGOslim <- merge(STACKED_silo3_pro_GOid_term, sr_lab_goslim, by = "GO", all.x = TRUE)
  
  STACKED_silo3_pro_GOid_term_sr_labGOslim$GOSlim_bin <- ifelse(is.na(STACKED_silo3_pro_GOid_term_sr_labGOslim$GOSlim_bin),"unmapped",STACKED_silo3_pro_GOid_term_sr_labGOslim$GOSlim_bin)
  
  silo3_all_srlab_terms <- data.frame(table(STACKED_silo3_pro_GOid_term_sr_labGOslim$GOSlim_bin))
  silo3_all_srlab_terms$day <- i
  View(silo3_all_srlab_terms)
  
  silo3_all_srlab_BPterms <- data.frame(table(STACKED_silo3_pro_GOid_term_sr_labGOslim[which(STACKED_silo3_pro_GOid_term_sr_labGOslim$aspect=="P"),"GOSlim_bin"]))
  silo3_all_srlab_BPterms$day <- i
    
  silo3_all_srlab_terms_per_day <- rbind(silo3_all_srlab_terms_per_day,silo3_all_srlab_terms)
  silo3_all_srlab_BPterms_per_day <- rbind(silo3_all_srlab_BPterms_per_day,silo3_all_srlab_BPterms)
}

#compare GO slim frequencies between temperatures
compare_terms <- merge(silo3_all_srlab_BPterms_per_day, silo9_all_srlab_BPterms_per_day, by = c("Var1", "day"))
colnames(compare_terms) <- c("GOslim", "day", "23C", "29C")
print(compare_terms)

##make pie chart for each day
#make pie charts of sr lab GO slim terms for all mapped proteins
#for color palette expansion: https://www.r-bloggers.com/how-to-expand-color-palette-with-ggplot-and-rcolorbrewer/
#for pie chart: http://www.sthda.com/english/wiki/ggplot2-pie-chart-quick-start-guide-r-software-and-data-visualization


getPalette = colorRampPalette(brewer.pal(11, "Spectral"))
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

for(i in days){
p <- ggplot(silo3_all_srlab_BPterms_per_day[which(silo3_all_srlab_BPterms_per_day$day== i),], aes(x = "", y = Freq, fill = Var1)) + geom_bar(stat = "identity") + coord_polar("y", start = 0) + scale_fill_manual(values = getPalette(nrow(silo3_all_srlab_BPterms_per_day[which(silo3_all_srlab_BPterms_per_day$day== i),]))) + blank_theme + theme(axis.text.x = element_blank())
ggsave(paste0("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/pieS3",i,".pdf"),p)
}

    
    

###########
###All for silo 9
########
#reformat to long list
silo9_t <- data.frame(t(silo9[,-c(1,2)]))
colnames(silo9_t) <- paste0("D",silo9$day)

#remove proteins that were only in silo9
View(silo9_t[which(apply(silo9_t, 1, var) == 0),])
no_val_proteins <- rownames(silo9_t[which(apply(silo9_t, 1, var) == 0),])

View(silo9_t[-c(which(rownames(silo9_t) %in% no_val_proteins)),])
silo9_t_nozero <- silo9_t[-c(which(rownames(silo9_t) %in% no_val_proteins)),]
#confirm it worked
nrow(silo9_t)-nrow(silo9_t_nozero)
length(no_val_proteins)

#add protein column
silo9_t_nozero$protein_ID <- rownames(silo9_t_nozero)
silo9_t_nozero <- silo9_t_nozero[grep("CHOYP", silo9_t_nozero$protein_ID),]

#merge with uniprot
#uniprot <- read.csv("/Volumes/web/metacarcinus/Cgigas/all_giga-uniprot-blastP-out.nopipe.annotations.tab", sep ="\t", header = FALSE, stringsAsFactors = FALSE)
#combine uniprot and NSAF data
silo9_t_nozero_unip<- merge(silo9_t_nozero, uniprot, by = "protein_ID", all.x = TRUE)
#exclude proteins that didn't map to uniprot DB
silo9_t_nozero_unip_mapped <- silo9_t_nozero_unip[-grep("unmapped", silo9_t_nozero_unip$Entry),]
nrow(silo9_t_nozero_unip_mapped)
silo9_t_nozero_unip_mapped <- silo9_t_nozero_unip_mapped[which(silo9_t_nozero_unip_mapped$evalue <= 10^-10),]
nrow(silo9_t_nozero_unip_mapped)


##get go terms for each day
silo9_all_srlab_terms_per_day <- data.frame()
silo9_all_srlab_BPterms_per_day <- data.frame()

days <- colnames(silo9_t_nozero_unip_mapped[,2:7])
for(i in days){
  silo9_pro_GO <- silo9_t_nozero_unip_mapped[which(silo9_t_nozero_unip_mapped[,i] != 0),c("protein_ID","GO_IDs")]
  
  silo9_pro_GOid_term <- data.frame()
  #loop through each line of data frame containing one column with "protein_ID" and one column with list of "GO_IDs" separated by ";"
  for (j in 1:nrow(silo9_pro_GO)){
    #create a row that lists all GO IDs associated with one protein listed across rows, where each different GO ID gets put in it's own column
    silo9_pro_GOid_term_row <- data.frame(t(data.frame(strsplit(as.character(silo9_pro_GO$GO_IDs[j]),'; ', fixed = TRUE))), stringsAsFactors = FALSE)
    #add each row created in the line above to the empty data frame
    #this will add NAs in columns where less terms exist; e.g. if a protein only has two terms and another has 10, there will be 8 NAs added to the row for the protein with 2 terms
    silo9_pro_GOid_term <- rbind.fill(silo9_pro_GOid_term,silo9_pro_GOid_term_row)
  }
  
  #add protein IDs back to GO IDs
  silo9_pro_GOid_term <- cbind(silo9_pro_GO[,"protein_ID"], silo9_pro_GOid_term)
  #this results in a table with protein ID listed in one column and each next column contains a GO ID that was listed with the protein in Uniprot DB.
  
  #reshape data so that all GO ID columns are gathered in one column called "GO" 
  STACKED_silo9_pro_GOid_term <- tidyr::gather(silo9_pro_GOid_term,"protein_ID","GO", 2:ncol(silo9_pro_GOid_term))
  #exlude middle column which just contains the string "protein_ID" in each row
  STACKED_silo9_pro_GOid_term <- STACKED_silo9_pro_GOid_term[,c(1,3)]
  #remove duplicate rows
  STACKED_silo9_pro_GOid_term <- unique(STACKED_silo9_pro_GOid_term)
  colnames(STACKED_silo9_pro_GOid_term)[1] <- "protein_ID"
  #remove any rows where GO column has NA value. 
  STACKED_silo9_pro_GOid_term <- STACKED_silo9_pro_GOid_term[which(!is.na(STACKED_silo9_pro_GOid_term$GO)),]
  #this resulting data frame has two columns "protein_ID" and "GO"
  
  ##get go slim terms for each day
  STACKED_silo9_pro_GOid_term_sr_labGOslim <- merge(STACKED_silo9_pro_GOid_term, sr_lab_goslim, by = "GO", all.x = TRUE)
  
  STACKED_silo9_pro_GOid_term_sr_labGOslim$GOSlim_bin <- ifelse(is.na(STACKED_silo9_pro_GOid_term_sr_labGOslim$GOSlim_bin),"unmapped",STACKED_silo9_pro_GOid_term_sr_labGOslim$GOSlim_bin)
  
  silo9_all_srlab_terms <- data.frame(table(STACKED_silo9_pro_GOid_term_sr_labGOslim$GOSlim_bin))
  silo9_all_srlab_terms$day <- i
  View(silo9_all_srlab_terms)
  
  silo9_all_srlab_BPterms <- data.frame(table(STACKED_silo9_pro_GOid_term_sr_labGOslim[which(STACKED_silo9_pro_GOid_term_sr_labGOslim$aspect=="P"),"GOSlim_bin"]))
  silo9_all_srlab_BPterms$day <- i
  
  silo9_all_srlab_terms_per_day <- rbind(silo9_all_srlab_terms_per_day,silo9_all_srlab_terms)
  silo9_all_srlab_BPterms_per_day <- rbind(silo9_all_srlab_BPterms_per_day,silo9_all_srlab_BPterms)
}
##make pie chart for each day
#make pie charts of sr lab GO slim terms for all mapped proteins
#for color palette expansion: https://www.r-bloggers.com/how-to-expand-color-palette-with-ggplot-and-rcolorbrewer/
#for pie chart: http://www.sthda.com/english/wiki/ggplot2-pie-chart-quick-start-guide-r-software-and-data-visualization


getPalette = colorRampPalette(brewer.pal(11, "Spectral"))
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

for(i in days){
  p <- ggplot(silo9_all_srlab_BPterms_per_day[which(silo9_all_srlab_BPterms_per_day$day== i),], aes(x = "", y = Freq, fill = Var1)) + geom_bar(stat = "identity") + coord_polar("y", start = 0) + scale_fill_manual(values = getPalette(nrow(silo9_all_srlab_BPterms_per_day[which(silo9_all_srlab_BPterms_per_day$day== i),]))) + blank_theme + theme(axis.text.x = element_blank())
  ggsave(paste0("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/pieS9",i,".pdf"),p)
}

    
    
    
    
    
    
#go terms
#Subset "protein_ID" column and "GO IDs" column containing list of GO_IDs separated by ";" from all_sig0.1_pro_logFC_pval 
silo3_pro_GO <- silo3_t_nozero_unip_mapped[,c("protein_ID","GO_IDs")]

silo3_pro_GOid_term <- data.frame()
#loop through each line of data frame containing one column with "protein_ID" and one column with list of "GO_IDs" separated by ";"
for (i in 1:nrow(silo3_pro_GO)){
  #create a row that lists all GO IDs associated with one protein listed across rows, where each different GO ID gets put in it's own column
  silo3_pro_GOid_term_row <- data.frame(t(data.frame(strsplit(as.character(silo3_pro_GO$GO_IDs[i]),'; ', fixed = TRUE))), stringsAsFactors = FALSE)
  #add each row created in the line above to the empty data frame
  #this will add NAs in columns where less terms exist; e.g. if a protein only has two terms and another has 10, there will be 8 NAs added to the row for the protein with 2 terms
  silo3_pro_GOid_term <- rbind.fill(silo3_pro_GOid_term,silo3_pro_GOid_term_row)
}

#add protein IDs back to GO IDs
silo3_pro_GOid_term <- cbind(silo3_t_nozero_unip_mapped[,"protein_ID"], silo3_pro_GOid_term)
#this results in a table with protein ID listed in one column and each next column contains a GO ID that was listed with the protein in Uniprot DB.
str(silo3_pro_GOid_term)

#reshape data so that all GO ID columns are gathered in one column called "GO" 
STACKED_silo3_pro_GOid_term <- tidyr::gather(silo3_pro_GOid_term,"protein_ID","GO", 2:ncol(silo3_pro_GOid_term))
#exlude middle column which just contains the string "protein_ID" in each row
STACKED_silo3_pro_GOid_term <- STACKED_silo3_pro_GOid_term[,c(1,3)]
#remove duplicate rows
STACKED_silo3_pro_GOid_term <- unique(STACKED_silo3_pro_GOid_term)
colnames(STACKED_silo3_pro_GOid_term)[1] <- "protein_ID"
#remove any rows where GO column has NA value. 
STACKED_silo3_pro_GOid_term <- STACKED_silo3_pro_GOid_term[which(!is.na(STACKED_silo3_pro_GOid_term$GO)),]
#this resulting data frame has two columns "protein_ID" and "GO"




###Next map all GO IDs to GO slim terms
#make list of unique GO terms without a protein ID column
GOids <- unique(STACKED_silo3_pro_GOid_term$GO)
#goslims with GSEA
myCollection <- GOCollection(GOids)
#I downloaded goslim_generic.obo from http://geneontology.org/docs/go-subset-guide/
#then i moved it to the R library for GSEABase in the extdata folder
fl <- system.file("extdata", "goslim_generic.obo", package="GSEABase")
slim <- getOBOCollection(fl)
slims <- data.frame(goSlim(myCollection, slim, "BP"))
slims$GOid <- rownames(slims)
slims$Term <- as.character(slims$Term)
rownames(slims) <- NULL
#this results in a data frame with columns: 
# "GOid" which is the GO slim GO ID
# "Term" which is the GO slim term corresponding to the GO ID
# "Percent" which is the percent of GO IDs in my list that mapped to each GO slim term
# "Count" which is the number of GO IDs in my list that mapped to each GO slim term

write.csv(slims, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/silo3_goslims.csv", row.names = FALSE, quote = FALSE)


#########
#For silo 9

#reformat to long list
silo9_t <- data.frame(t(silo9[,-c(1,2)]))
colnames(silo9_t) <- silo9$day

#remove proteins that were only in silo9
View(silo9_t[which(apply(silo9_t, 1, var) == 0),])
no_val_proteins <- rownames(silo9_t[which(apply(silo9_t, 1, var) == 0),])

View(silo9_t[-c(which(rownames(silo9_t) %in% no_val_proteins)),])
silo9_t_nozero <- silo9_t[-c(which(rownames(silo9_t) %in% no_val_proteins)),]
#confirm it worked
nrow(silo9_t)-nrow(silo9_t_nozero)
length(no_val_proteins)

#add protein column
silo9_t_nozero$protein_ID <- rownames(silo9_t_nozero)
silo9_t_nozero <- silo9_t_nozero[grep("CHOYP", silo9_t_nozero$protein_ID),]

#merge with uniprot
#uniprot <- read.csv("/Volumes/web/metacarcinus/Cgigas/all_giga-uniprot-blastP-out.nopipe.annotations.tab", sep ="\t", header = FALSE, stringsAsFactors = FALSE)
#combine uniprot and NSAF data
silo9_t_nozero_unip<- merge(silo9_t_nozero, uniprot, by = "protein_ID", all.x = TRUE)
#exclude proteins that didn't map to uniprot DB
silo9_t_nozero_unip_mapped <- silo9_t_nozero_unip[-grep("unmapped", silo9_t_nozero_unip$Entry),]
nrow(silo9_t_nozero_unip_mapped)
silo9_t_nozero_unip_mapped <- silo9_t_nozero_unip_mapped[which(silo9_t_nozero_unip_mapped$evalue <= 10^-10),]
nrow(silo9_t_nozero_unip_mapped)



#go terms
#Subset "protein_ID" column and "GO IDs" column containing list of GO_IDs separated by ";" from all_sig0.1_pro_logFC_pval 
silo9_pro_GO <- silo9_t_nozero_unip_mapped[,c("protein_ID","GO_IDs")]

silo9_pro_GOid_term <- data.frame()
#loop through each line of data frame containing one column with "protein_ID" and one column with list of "GO_IDs" separated by ";"
for (i in 1:nrow(silo9_pro_GO)){
  #create a row that lists all GO IDs associated with one protein listed across rows, where each different GO ID gets put in it's own column
  silo9_pro_GOid_term_row <- data.frame(t(data.frame(strsplit(as.character(silo9_pro_GO$GO_IDs[i]),'; ', fixed = TRUE))), stringsAsFactors = FALSE)
  #add each row created in the line above to the empty data frame
  #this will add NAs in columns where less terms exist; e.g. if a protein only has two terms and another has 10, there will be 8 NAs added to the row for the protein with 2 terms
  silo9_pro_GOid_term <- rbind.fill(silo9_pro_GOid_term,silo9_pro_GOid_term_row)
}

#add protein IDs back to GO IDs
silo9_pro_GOid_term <- cbind(silo9_t_nozero_unip_mapped[,"protein_ID"], silo9_pro_GOid_term)
#this results in a table with protein ID listed in one column and each next column contains a GO ID that was listed with the protein in Uniprot DB.
str(silo9_pro_GOid_term)

#reshape data so that all GO ID columns are gathered in one column called "GO" 
STACKED_silo9_pro_GOid_term <- tidyr::gather(silo9_pro_GOid_term,"protein_ID","GO", 2:ncol(silo9_pro_GOid_term))
#exlude middle column which just contains the string "protein_ID" in each row
STACKED_silo9_pro_GOid_term <- STACKED_silo9_pro_GOid_term[,c(1,3)]
#remove duplicate rows
STACKED_silo9_pro_GOid_term <- unique(STACKED_silo9_pro_GOid_term)
colnames(STACKED_silo9_pro_GOid_term)[1] <- "protein_ID"
#remove any rows where GO column has NA value. 
STACKED_silo9_pro_GOid_term <- STACKED_silo9_pro_GOid_term[which(!is.na(STACKED_silo9_pro_GOid_term$GO)),]
#this resulting data frame has two columns "protein_ID" and "GO"


STACKED_silo9_pro_GOid_term_sr_labGOslim <- merge(STACKED_silo9_pro_GOid_term, sr_lab_goslim, by = "GO", all.x = TRUE)

STACKED_silo9_pro_GOid_term_sr_labGOslim$GOSlim_bin <- ifelse(is.na(STACKED_silo9_pro_GOid_term_sr_labGOslim$GOSlim_bin),"unmapped",STACKED_silo9_pro_GOid_term_sr_labGOslim$GOSlim_bin)

silo9_all_srlab_terms <- data.frame(table(STACKED_silo9_pro_GOid_term_sr_labGOslim$GOSlim_bin))
silo9_all_srlab_BPterms <- data.frame(table(STACKED_silo9_pro_GOid_term_sr_labGOslim[which(STACKED_silo9_pro_GOid_term_sr_labGOslim$aspect=="P"),"GOSlim_bin"]))



###Next map all GO IDs to GO slim9 terms
#make list of unique GO terms without a protein ID column
GOids <- unique(STACKED_pro_GOid_term$GO)
#goslim9s with GSEA
myCollection <- GOCollection(GOids)
#I downloaded goslim9_generic.obo from http://geneontology.org/docs/go-subset-guide/
#then i moved it to the R library for GSEABase in the extdata folder
fl <- system.file("extdata", "goslim_generic.obo", package="GSEABase")
slim <- getOBOCollection(fl)
slim9s <- data.frame(goSlim(myCollection, slim, "BP"))
slim9s$GOid <- rownames(slim9s)
slim9s$Term <- as.character(slim9s$Term)
rownames(slim9s) <- NULL
#this results in a data frame with columns: 
# "GOid" which is the GO slim9 GO ID
# "Term" which is the GO slim9 term corresponding to the GO ID
# "Percent" which is the percent of GO IDs in my list that mapped to each GO slim9 term
# "Count" which is the number of GO IDs in my list that mapped to each GO slim9 term

write.csv(slim9s, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/silo9_goslims.csv", row.names = FALSE, quote = FALSE)


sr_lab_goslim <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/background/GOSlim_terms.csv", stringsAsFactors = FALSE)





####################
#try minimal set




#try get ancestors and remove redundant terms
###this creates a list of ancestor GO IDs for each GO ID in my list
##if the GO ID is not in the "go" object from the package, the entry is "NULL"
#https://jonlefcheck.net/2013/05/20/continuing-a-for-loop-in-r-after-an-error/

term_prop <- list() # create an empty list that will get filled in by loop
for(i in 1:length(sig_GOids)){ # for each line in my GO IDs list
  temp_term_prop <- try(go$ancestors[[sig_GOids[i]]], TRUE) #make a list of all ancester GO IDs for each GO ID in my list
  if(isTRUE(class(temp_term_prop)=="try-error")) {next} else {term_prop[[i]] = temp_term_prop} # if the "go" data doesn't contain my GO ID, go on to the next GO ID. 
}

#create an empty data frame the length the list of ancestor list made above and with two columns
ancestors <- data.frame(matrix(0,length(term_prop),2))
#names the two columns
colnames(ancestors) <- c("orig.GO","GOan")

#make an empty data frame to get filled in by the loop
ances_test <- data.frame()
for(i in 1:length(term_prop)){ #for each GO ID in the ancestors list (which is the same length and order as the sig_GOids list)
  ancestors$orig.GO[i] <- sig_GOids[i] # fill in the orig. GO ID column 
  ancestors$GOan[i] <- paste(term_prop[[i]], collapse = "_") # fill in the GO ancestor column with all GO ancestor IDs corresponding to original GO term separated by an underscore
  ancestors_row <- data.frame(t(data.frame(strsplit(as.character(ancestors$GOan[i]),'_', fixed = TRUE)))) #spread ancestor IDs out across multiple columns
  ances_test <- rbind.fill(ances_test,ancestors_row) #add each row to the data frame
}
#convert factors to characters
ances_test <- data.frame(lapply(ances_test, as.character), stringsAsFactors = FALSE)
#add original GO IDs back to ancestor GO IDs
ances_test <- cbind(data.frame(sig_GOids, stringsAsFactors = FALSE), ances_test)

###list GO IDs not in 'go' object (don't have ancestors)
length(ances_test[which(is.na(ances_test$X1)),1])
#[1] GO:0062023 GO:0103025 GO:0102102 GO:0102131 GO:0090736 GO:1905905 GO:0061844 GO:1905907 GO:0106036

#reshape data so that all ancestor GO ID columns are gathered in one column called "Ancterm" 
STACKED_anc_term <- tidyr::gather(ances_test,"Ancestor","Ancterm", 2:ncol(ances_test))
STACKED_anc_term <- STACKED_anc_term[,c(1,3)]
STACKED_anc_term <- unique(STACKED_anc_term)
STACKED_anc_term <- STACKED_anc_term[which(!is.na(STACKED_anc_term$Ancterm)),]
