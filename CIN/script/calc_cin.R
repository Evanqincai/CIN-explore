library(tidyverse)
args <- commandArgs(T)
corrected_cn_file <- args[1]
corrected_cn <- read.table(corrected_cn_file, header=T)
baseline <- read.table('/lustre/rde/user/lvf/project/HRD/ULP/Trimmed/baseline/baseline_result/baseline/50M_baseline',header=T,row.names=1)
#id      chrom_50bin_mean
baseline$id <- row.names(baseline)
baseline <- subset(baseline , select=c('id','mean','sd'))
hit_region <- corrected_cn %>% inner_join(baseline, by='id')
instability_region <- hit_region[(hit_region$chrom_50bin_mean > hit_region$mean + 3*hit_region$sd ) | (hit_region$chrom_50bin_mean < hit_region$mean - 3*hit_region$sd ),]
instability_region_number <- length(instability_region$id)
hit_region$zscore <- (hit_region$chrom_50bin_mean - hit_region$mean)/hit_region$sd
genome_wide_zscore <- (sum(hit_region$zscore^2) - baseline['g-wide', 'mean'])/baseline['g-wide', 'sd']
genome_wide_zscore <- round(genome_wide_zscore, digits=2)
sample_id <- gsub(".correctedDepth.txt.merge_bin.txt", "", corrected_cn_file)
res <- paste(sample_id, instability_region_number, genome_wide_zscore)
print(res)
q()
