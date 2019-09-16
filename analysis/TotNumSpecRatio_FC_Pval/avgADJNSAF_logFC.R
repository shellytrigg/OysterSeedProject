

#read in avg NSAF data
avgNSAF_data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/silo3and9_nozerovals_AVGs.csv", stringsAsFactors = FALSE)
#exclude day 0
avgNSAF_data <- avgNSAF_data[-1,]

#order by day and temp
avgNSAF_data <- avgNSAF_data[order(avgNSAF_data$day,avgNSAF_data$temp),]

#calculate log FC NSAF values
logFC_NSAF <- data.frame()
for (i in seq(1,nrow(avgNSAF_data),2)){
  for(j in 1:ncol(avgNSAF_data[,-c(1:2)])){
    logFC_NSAF[i,j] <- log(avgNSAF_data[i+1,j+2],2) - log(avgNSAF_data[i,j+2],2)
  }
}

#remove empty rows
logFC_NSAF <- logFC_NSAF[which(!is.na(logFC_NSAF$V1)),]

#bind with day column
logFC_NSAF <- cbind(unique(avgNSAF_data[,"day"]), logFC_NSAF)
#rename columns
colnames(logFC_NSAF) <- c("day",colnames(avgNSAF_data[,-c(1:2)]))

#transform data

logFC_NSAF_t <- data.frame(t(logFC_NSAF[,-1]),stringsAsFactors = FALSE)

#make protein column
logFC_NSAF_t$protein_ID <- rownames(logFC_NSAF_t)

#rename columns
colnames(logFC_NSAF_t) <- c("D_3_logFC_NSAF", "D_5_logFC_NSAF","D_7_logFC_NSAF","D_9_logFC_NSAF", "D_11_logFC_NSAF","D_13_logFC_NSAF", "protein_ID")

#write output to file
write.csv(logFC_NSAF_t, "~/Documents/GitHub/OysterSeedProject/analysis/TotNumSpecRatio_FC_Pval/avgADJNSAF_logFC_DAYSCOMPARED.csv", row.names = FALSE, quote = FALSE)
