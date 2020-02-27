if [ -z $3 ]; then
    echo "bash $(basename $0) WIG_FILES_DIR WIG_FILE_LIST OUT_FILE"
    exit 1
fi
set -e

WIG_FILES_DIR=$1
OUT_FILE=$3
WIG_FILE_LIST=$2
ls ${WIG_FILES_DIR}/*wig > ${WIG_FILE_LIST}
/lustre/rde/user/lvf/project/HRD/BRCA_classify/LP_WGS/until_1016/R-3.6.1/bin/Rscript /lustre/rde/user/lvf/project/HRD/BRCA_classify/LP_WGS/until_1016/ichorCNA/scripts/createPanelOfNormals.R \
   --filelist  ${WIG_FILE_LIST} \
   --gcWig  /lustre/rde/user/lvf/project/HRD/BRCA_classify/LP_WGS/until_1016/ichorCNA/inst/extdata/gc_hg19_1000kb.wig \
   --centromere /lustre/rde/user/lvf/project/HRD/BRCA_classify/LP_WGS/until_1016/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
   --outfile ${OUT_FILE}
BASELINE=`realpath ${OUT_FILE}`
echo "Your baseline is ${BASELINE}.rds"
