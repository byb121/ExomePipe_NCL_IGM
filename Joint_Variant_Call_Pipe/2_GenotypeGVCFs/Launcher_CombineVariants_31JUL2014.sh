#! /bin/bash

#**********************************  important!!  *******************************************
#if use the script after trimming: fastq file name patterns need to be changed as following: 
# comment out - "map_recali_perLane_recali_perSample_covOnTargets_26MAY2014.sh" line 112, 113
# and use line 108, 109 
# the comment out - "detectSampleLanes.pl" line 13 use line 16
#********************************************************************************************

## Script to launch Exome Pipeline (Picard - GATK bam recalibration - GATK call variants & varscan & dindel & bedtools for coverage) for multiple samples
## INALL is an array of directory names, with one directory per sample
## SAMPLE_PATH is the pathway to the sample directories ${SAMPLE_PATH}/${INALL[$SAMPLE_NUMBER]}


#/************************* Paras need to be adjusted for diff samples
INPUT="/home/a5907529/WORKING_DATA/Sophie_A2463/scripts/GenotypeGVCFs/HC_out" #folder of GenotypeGVCFs output
SCRIPT_PATH="/home/a5907529/WORKING_DATA/Sophie_A2463/scripts/GenotypeGVCFs"
REF_FILE="/sharedlustre/IGM/GATK/bundle2.8/b37/human_g1k_v37_decoy.fasta"
OUTPUT_FILE="/home/a5907529/WORKING_DATA/Sophie_A2463/scripts/GenotypeGVCFs/HC_out_combined.vcf" #gatk tuned alignment file output directory name
JAVA_TMP_DIR_NAME="_hulala" #tail of temp java working directory under scratch dir
#**************************/

###############################
## Submitting jobs ##
JOB_ID1="Yaobo_Exome_CombineVariants"
echo $JOB_ID1
arr=("${REF_FILE}" "${INPUT}" "${OUTPUT_FILE}" "${JOB_ID1}" "$JAVA_TMP_DIR_NAME")
#echo qsub -N "${JOB_ID1}" ${SCRIPT_PATH}/CombineVariants.sh "${arr[@]}"
qsub -N "${JOB_ID1}" ${SCRIPT_PATH}/CombineVariants.sh "${arr[@]}"
#sh ${SCRIPT_PATH}/CombineVariants.sh "${arr[@]}"


