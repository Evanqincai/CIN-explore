files <- Sys.glob("*.correctedDepth.txt.merge_bin.txt")
library(tidyverse)
library(matrixStats)
library(robustHD)
mat <- read.table(files[1], header=T)
files[1] <- gsub(".correctedDepth.txt.merge_bin.txt","",files[1])
names(mat)[2] <- files[1]
mat <- subset(mat, select=c('id',files[1]))
for(f in files[-1]){
 d <- read.table(f, header=T)
#id      chrom_50bin_mean
 f <- gsub(".correctedDepth.txt.merge_bin.txt","",f)
 names(d)[2] <- f
 mat <- mat %>% inner_join(d, by="id")
}
mat <- as.data.frame(mat)
row.names(mat) <- mat$id
mat <- subset(mat ,select=-c(id))
mean <-  (as.data.frame(rowMeans(mat)))[,1]
sd <- rowSds(as.matrix(mat))
l <- length(colnames(mat))
mat$mean <- mean
mat$sd <- sd
zscore_mat <- mat[,1:l]
for(i in colnames(zscore_mat)){
 for(j in row.names(zscore_mat)){
  zscore_mat[j,i] <- (zscore_mat[j,i] - mat[j,'mean'])/mat[j,'sd']
 }
}
tmp <- zscore_mat^2
g_zscore <- colSums(tmp)
g_zscore <- c(g_zscore, mean(g_zscore), sd(g_zscore))
zscore_mat$mean <- mean
zscore_mat$sd <- sd
zscore_mat['g-wide',] <- g_zscore
write.table(zscore_mat, 'new_all.txt', sep="\t", quote=FALSE,row.names=T)
q()
