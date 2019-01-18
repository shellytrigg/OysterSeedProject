#check if averaged NSAF values correlate with summed NSAF values

avg <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/silo3and9_nozerovals_AVGs.csv", stringsAsFactors = FALSE)
sum <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/silo3and9_nozerovals_SUMs.csv", stringsAsFactors = FALSE)

avg_sum <- cbind(avg,sum)
rownames(avg_sum) <- paste(avg$day,avg$temp, sep ="_")
avg_sum <- avg_sum[,order(colnames(avg_sum))]

#create an empty data frame
df_cor <- data.frame(matrix(0,nrow(avg_sum),0))
#loop through the data and calculate the standard deviation between replicates for each protein
#this calculates the SD for 
for (i in seq(1,ncol(avg_sum),2)){
  #this calculates the SD for each odd number row and the row following it
  df_cor_col <- cor(avg_sum[,i],avg_sum[,i+1])
  #this sequencially combines rows of data together after the SD is generated
  df_cor <- cbind(df_cor, df_cor_col)
}

#make a frequency table of correlation coefficients
data.frame(table(data.frame(unlist(df_cor))))

#They are all 1 so avg and sum NSAF values are 100% correlated