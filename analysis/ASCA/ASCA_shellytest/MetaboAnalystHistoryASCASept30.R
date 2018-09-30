mSet<-InitDataObjects("conc", "ts", FALSE)
mSet<-SetDesignType(mSet, "time")
mSet<-Read.TextData(mSet, "Replacing_with_your_file_path", "rowts", "disc");
mSet<-SanityCheckData(mSet)
#remove proteins if they aren't detected in >= 75% of samples
mSet<-RemoveMissingPercent(mSet, percent=0.75)
#replace NA values (where protein undetected) with column minimum
#this may not be the right thing to do...we need to figure out what this should be
mSet<-ImputeVar(mSet, method="colmin")
#chose IQR filter
mSet<-FilterVariable(mSet, "iqr", "F", 25)
#selected autoscaling because the non-normalized data looked very skewed to the left
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
mSet<-Normalization(mSet, "NULL", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_1_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_1_", "png", 72, width=NA)
mSet<-Perform.ASCA(mSet, 1, 1, 2, 2)
mSet<-PlotModelScree(mSet, "asca_scree_0_", "png", 72, width=NA)
mSet<-PlotASCAModel(mSet, "asca_fa_0_", "png", 72, width=NA, "a",FALSE)
mSet<-PlotASCAModel(mSet, "asca_fb_0_", "png", 72, width=NA, "b",FALSE)
mSet<-PlotInteraction(mSet, "asca_fab_0_", "png", 72,FALSE, width=NA)
mSet<-Perform.ASCA.permute(mSet, 20)
mSet<-PlotASCA.Permutation(mSet, "asca_perm_0_", "png", 72, width=NA)
mSet<-CalculateImpVarCutoff(mSet, 0.05, 0.9)
mSet<-PlotAscaImpVar(mSet, "asca_impa_0_", "png", 72, width=NA, "a")
mSet<-PlotAscaImpVar(mSet, "asca_impb_0_", "png", 72, width=NA, "b")
mSet<-PlotAscaImpVar(mSet, "asca_impab_0_", "png", 72, width=NA, "ab")
