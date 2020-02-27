#! /bin/bash
path_in=$1
path_out=$2

for i in 1 2.5 5;
do
    mkdir ${path_out}/${i}/Defi
    mkdir ${path_out}/${i}/Defi/{coverage,coverage_shell}
    cat ${path_out}/${i}/sample.id|while read  ID; do
        echo "/lustre/rde/user/lvf/tools/anaconda2/bin/bedtools multicov -bams ${path_in}/${i}/${ID}.20200214/basic_analysis/align/${ID}.sorted.rmdup.realign.bam  -bed /lustre/rde/user/lvf/project/CIN/dilution/region_bed/1M.bed > ${path_out}/${i}/Defi/coverage/${ID}.sorted.rmdup.realign.bam.1M.coverage 
	/lustre/rde/user/lvf/tools/anaconda2/bin/bedtools nuc -fi /lustre/rde/user/lvf/database/hg19/hg19.fa -bed ${path_out}/${i}/Defi/coverage/${ID}.sorted.rmdup.realign.bam.1M.coverage | tail -n +2 | cut -f 1-4,6 >  ${path_out}/${i}/Defi/coverage/${ID}.sorted.rmdup.realign.bam.1M.coverage.gc
	/lustre/rde/user/lvf/tools/software/Rscript /lustre/rde/user/lvf/project/CIN/dilution/gc_correct.R -c ${path_out}/${i}/Defi/coverage/${ID}.sorted.rmdup.realign.bam.1M.coverage.gc -o ${path_out}/${i}/Defi/coverage/${ID}.sorted.rmdup.realign.bam.1M.coverage.gc.corrected " >${path_out}/${i}/Defi/coverage_shell/${ID}.coverage.sh
    done;
done
