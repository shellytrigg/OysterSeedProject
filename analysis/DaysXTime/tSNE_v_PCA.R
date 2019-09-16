###Run tSNE on day*time correct NSAF data

tsne_out <- Rtsne(NSAF_normTxD[,-c(1:2)],pca=FALSE,perplexity=3,theta=0.0)
ggplot(data.frame(tsne_out$Y), aes(data.frame(tsne_out$Y)[,1], data.frame(tsne_out$Y)[,2])) + geom_point(aes(col = pca_meta$day, shape = pca_meta$temp, size = 0.5)) + theme_bw() + ggtitle("tSNE of ADJNSAF values scaled by day*time")

### Run tSNE on uncorrected NSAF data
tsne_out <- Rtsne(NSAF_matrix,pca=FALSE,perplexity=3,theta=0.0)
ggplot(data.frame(tsne_out$Y), aes(data.frame(tsne_out$Y)[,1], data.frame(tsne_out$Y)[,2])) + geom_point(aes(col = pca_meta$day, shape = pca_meta$temp, size = 0.5)) + theme_bw() + ggtitle("tSNE of ADJNSAF values")

### Run tSNE on log2 transformed NSAF data
tsne_out <- Rtsne(log(NSAF_matrix,2),pca=FALSE,perplexity=3,theta=0.0)
ggplot(data.frame(tsne_out$Y), aes(data.frame(tsne_out$Y)[,1], data.frame(tsne_out$Y)[,2])) + geom_point(aes(col = pca_meta$day, shape = pca_meta$temp, size = 0.5)) + theme_bw() + ggtitle("tSNE of log2 ADJNSAF values")

### Run PCA on log2 transformed NSAF data
pca <- prcomp(log(NSAF_matrix,2))
ggplot(data.frame(pca$x), aes(PC1, PC2)) + geom_point(aes(col = pca_meta$day, shape = pca_meta$temp, size = 0.5)) + theme_bw() + ggtitle("PCA of log2 ADJNSAF values")

### Run PCA on uncorrected NSAF data
pca <- prcomp(NSAF_matrix,2)
ggplot(data.frame(pca$x), aes(PC1, PC2)) + geom_point(aes(col = pca_meta$day, shape = pca_meta$temp, size = 0.5)) + theme_bw() + ggtitle("PCA of ADJNSAF values")

