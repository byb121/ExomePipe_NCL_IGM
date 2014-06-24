#! /bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
# #$ -l mem_free=40G
# $ -j h="node1[7-9]|node20"

# Load Modules
#/******************** on lampredi2
#module load apps/R/3.0.0/gcc-4.4.6+lapack-3.4.1+blas-1
#module load apps/gatk/3.1.1/noarch
#REF_DIR="/users/a5907529/data/GATK/bundle2.8/hg19"
#********************/

#/******************** on fmscluster
module load apps/java
export _JAVA_OPTIONS="-XX:-UseLargePages -Xms4096m"
GATKDIR="/home/a5907529/WORKING_DATA/QC_GOOD/GATK3"
REF_DIR="/sharedlustre/IGM/bundle2.8/hg19"
#********************/

SAMPLE_PATH=`pwd`
OUTPUT_VCF_PREFIX="New_HaplotypeCaller"
RAW_VCF_INPUT="$SAMPLE_PATH/HaplotyperCaller.vcf"

REF_FILE="${REF_DIR}/ucsc.hg19.fasta"
HAPMAP="${REF_DIR}/hapmap_3.3.hg19.vcf"
DBSNP="${REF_DIR}/dbsnp_138.hg19.vcf"
OMNI="${REF_DIR}/1000G_omni2.5.hg19.vcf"
MILLS_INDEL="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg19.vcf"
SNP_1000G="${REF_DIR}/1000G_phase1.snps.high_confidence.hg19.vcf"

SNP_RECAL_FILE="${RAW_VCF_INPUT}.SNP.recali"
SNP_TRANCH_FILE="${RAW_VCF_INPUT}.SNP.tranches"
SNP_RECAL_R_SCRIPT="${RAW_VCF_INPUT}.SNP.rscript"
SNP_RECALI_OUTPUT="$SAMPLE_PATH/${OUTPUT_VCF_PREFIX}_recali_SNP.vcf"

INDEL_RECAL_FILE="${RAW_VCF_INPUT}.INDEL.recali"
INDEL_TRANCH_FILE="${RAW_VCF_INPUT}.INDEL.tranches"
INDEL_RECAL_R_SCRIPT="${RAW_VCF_INPUT}.INDEL.rscript"
INDEL_RECALI_OUTPUT="$SAMPLE_PATH/${OUTPUT_VCF_PREFIX}_recali_SNP_INDEL.vcf"

# SNP error model
java -Xmx4g -jar $GATKDIR/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_FILE \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP} \
-resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 ${SNP_1000G} \
-an QD \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-nt 4 \
-mode SNP \
-input ${RAW_VCF_INPUT} \
-recalFile $SNP_RECAL_FILE \
-tranchesFile $SNP_TRANCH_FILE \
-rscriptFile $SNP_RECAL_R_SCRIPT

java -Xmx4g -jar $GATKDIR/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_FILE \
-input ${RAW_VCF_INPUT} \
--ts_filter_level 99.5 \
-recalFile $SNP_RECAL_FILE \
-tranchesFile $SNP_TRANCH_FILE \
-mode SNP \
-o $SNP_RECALI_OUTPUT

#indel model
java -Xmx4g -jar $GATKDIR/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_FILE \
--maxGaussians 4 \
-resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS_INDEL \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP} \
-an QD \
-an ReadPosRankSum \
-an MQRankSum \
-an FS \
-mode INDEL \
-nt 4 \
-input $RAW_VCF_INPUT \
-recalFile $INDEL_RECAL_FILE \
-tranchesFile $INDEL_TRANCH_FILE \
-rscriptFile $INDEL_RECAL_R_SCRIPT

java -Xmx4g -jar $GATKDIR/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_FILE \
-input $SNP_RECALI_OUTPUT \
--ts_filter_level 99.0 \
-recalFile $INDEL_RECAL_FILE \
-tranchesFile $INDEL_TRANCH_FILE \
-mode INDEL \
-o $INDEL_RECALI_OUTPUT