library(zoo)
args <- commandArgs(T)
in_file <- args[1]
out_file <- args[2]
d <- read.table(in_file , header=T)
genome_roll_mean <- data.frame()
for(i in 1:22){
 chrom <- d[which(d$chr==i),]
# chrom[is.na(chrom$log2_TNratio_corrected),]$log2_TNratio_corrected <- 0
 chrom <- na.omit(chrom)
 chrom$cn <- 2^chrom$log2_TNratio_corrected * 2
 chrom_zoo <- zoo(chrom$cn, chrom$start)
 chrom_roll_mean <- rollapply(chrom_zoo, 50, mean,  by=50, align = c("right"), partial = T)
 chrom_roll_mean <- as.data.frame(chrom_roll_mean)
 chrom_roll_mean$pos <- row.names(chrom_roll_mean)
 row.names(chrom_roll_mean) <- NULL
#chrom_roll_mean pos     chr
 chrom_roll_mean$chr <- i
 chrom_roll_mean$id <- paste(chrom_roll_mean$chr, chrom_roll_mean$pos, sep="_")
 chrom_roll_mean <- subset(chrom_roll_mean, select=c('id', 'chrom_roll_mean'))
 colnames(chrom_roll_mean) <- c('id', 'chrom_50bin_mean')
 genome_roll_mean <- rbind(genome_roll_mean, chrom_roll_mean)
}
write.table(genome_roll_mean, out_file , sep="\t", quote=FALSE, row.names=FALSE)
q()
