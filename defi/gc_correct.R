library(getopt)
command=matrix(c("coverage","c",1,"character",
                 "output","o",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if (!is.null(args$help) || is.null(args$coverage) ||  is.null(args$output)) {
    cat(paste(getopt(command, usage = T), "\n"))
    q()
}
coverage_file <- args$coverage
corrected_output <- args$output

library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Rsamtools)
class(Homo.sapiens)
library(devtools)
library(biovizBase)

load("/lustre/rde/user/lvf/project/CIN/delfi_scripts/filters.hg19.rda")
load("/lustre/rde/user/lvf/project/CIN/delfi_scripts/gaps.hg19.rda")
gc.correct <- function(coverage, bias) {
    bias_gc <- bias[coverage <  as.numeric(quantile(coverage, probs=c(0.1,0.9))[2]) & coverage > as.numeric(quantile(coverage ,probs=c(0.1,0.9))[1])]
    coverage_gc <- coverage[coverage <  as.numeric(quantile(coverage, probs=c(0.1,0.9))[2]) & coverage > as.numeric(quantile(coverage ,probs=c(0.1,0.9))[1])]
    i <- seq(min(bias_gc, na.rm=TRUE), max(bias_gc, na.rm=TRUE), by = 0.001)
#    coverage.trend <- loess(coverage_gc ~ bias_gc)
#    coverage.model <- loess(predict(coverage.trend, i) ~ i)
#    coverage.pred <- predict(coverage.model, bias)
#    coverage.corrected <- coverage  +  median(coverage) - coverage.pred
#    coverage.corrected
    coverage_loess <- loess(coverage_gc ~ bias_gc)
    coverage_pred <- predict(coverage_loess, bias)
    coverage.corrected <- coverage - coverage_pred + median(coverage)
}

gc.correct_delfi <- function(coverage, bias) {
#    bias_gc <- bias[coverage <  as.numeric(quantile(coverage, probs=c(0.1,0.9))[2]) & coverage > as.numeric(quantile(coverage ,probs=c(0.1,0.9))[1])]
#    coverage_gc <- coverage[coverage <  as.numeric(quantile(coverage, probs=c(0.1,0.9))[2]) & coverage > as.numeric(quantile(coverage ,probs=c(0.1,0.9))[1])]
    i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
    coverage.trend <- loess(coverage ~ bias)
    coverage.model <- loess(predict(coverage.trend, i) ~ i)
    coverage.pred <- predict(coverage.model, bias)
    coverage.corrected <- coverage +  median(coverage) - coverage.pred
    coverage.corrected
}



d <- as.data.frame(fread(coverage_file))
colnames(d) <- c("chr", "start", "end", "coverage", "gc_content")
#colnames(d) <- c("chr", "start", "end", "gc_content", "coverage")
log2read_count <- log2(d$coverage + 0.0001)
log2read_count_corrected <- gc.correct(log2read_count , d$gc_content*100)
d$coverage_corrected  <- round(2^log2read_count_corrected)
#d$coverage_corrected  <- gc.correct(d$coverage, d$gc_content * 100) 
d$coverage_corrected_delfi  <- gc.correct_delfi(d$coverage, d$gc_content * 100) 
read_count_corrected_range <- makeGRangesFromDataFrame(d, keep.extra.columns=T)
chromosomes <- GRanges(paste0("chr", 1:22),
                       IRanges(0, seqlengths(Hsapiens)[1:22]))
tcmeres <- gaps.hg19[grepl("centromere|telomere", gaps.hg19$type)]
arms <- GenomicRanges::setdiff(chromosomes, tcmeres)
arms <- arms[-c(25,27,29,41,43)]
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
arms$arm <- armlevels
read_count_corrected_range <- read_count_corrected_range[-queryHits(findOverlaps(read_count_corrected_range,gaps.hg19))]
read_count_corrected_range <- read_count_corrected_range[queryHits(findOverlaps(read_count_corrected_range, arms))]
read_count_corrected_range$arm <- armlevels[subjectHits(findOverlaps(read_count_corrected_range, arms))]
read_count_corrected_range <- data.frame(read_count_corrected_range)
read_count_corrected_range[is.na(read_count_corrected_range)] <- 0

loess_res <- loess(read_count_corrected_range$coverage ~ read_count_corrected_range$gc_content,control = loess.control(surface = "direct") ,degree = 2)
write.table(read_count_corrected_range,corrected_output,sep="\t",row.names=FALSE,quote=FALSE)
read_count_corrected_range <- read_count_corrected_range[read_count_corrected_range$coverage_corrected <  as.numeric(quantile(read_count_corrected_range$coverage_corrected , probs=c(0.1,0.9))[2]) & read_count_corrected_range$coverage_corrected > as.numeric(quantile(read_count_corrected_range$coverage_corrected ,probs=c(0.1,0.9))[1]),]
#loess_res <- loess(read_count_corrected_range$coverage ~ read_count_corrected_range$gc_content,control = loess.control(surface = "direct") ,degree = 2)
max_y <- as.numeric(max(read_count_corrected_range$coverage) * 1.1)
min_y <- as.numeric(min(read_count_corrected_range$coverage) * 0.8)
pred <- predict(loess_res, read_count_corrected_range$gc_content)
p1 <- paste0(corrected_output , ".png")
png(p1, type="cairo", height = 1000, width=1000)
par(mfrow=c(2,2))
plot(read_count_corrected_range$gc_content , read_count_corrected_range$coverage, xlab = "gc content", ylab = "read count",  cex.axis = 2,  cex.lab = 2, main="raw", cex.main = 3,ylim = c(min_y, max_y))
plot(read_count_corrected_range$gc_content , read_count_corrected_range$coverage_corrected_delfi, xlab = "gc content", ylab = "read count",  cex.axis = 2,  cex.lab = 2, main = "delfi" ,cex.main = 3 , ylim = c(min_y, max_y))
plot(read_count_corrected_range$gc_content , read_count_corrected_range$coverage, xlab = "gc content", ylab = "read count",  cex.axis = 2,  cex.lab = 2, main="raw", cex.main = 3, ylim = c(min_y, max_y))
lines(read_count_corrected_range$gc_content, pred, col = "red")
plot(read_count_corrected_range$gc_content , read_count_corrected_range$coverage_corrected, xlab = "gc content", ylab = "read count", cex.axis = 2,  cex.lab = 2, main = "delfi winsorize", cex.main = 3, ylim = c(min_y, max_y))

dev.off()
q()
