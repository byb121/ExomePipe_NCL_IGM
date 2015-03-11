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

SAMPLE_PATH="/home/a5907529/WORKING_DATA/Sophie_A2463"
SCRIPT_PATH="/home/a5907529/WORKING_DATA/Sophie_A2463/scripts/Align"
TARGETS="${SCRIPT_PATH}/refSeq.CDS.b37.2014Aug.bed"
Library_ID="Sophie.A2463.b1" #Better to use {PI name}.{Batch Number} to identify the library

REFDIR="/sharedlustre/IGM/GATK/bundle2.8/b37"
#SCRATCH_DIR="/users/a5907529/scratch"
SCRATCH_DIR="/users/a5907529/lustre/pipe.temp.sophie"
SEQ_PLATFORM="ILLUMINA" # Valid values are: ILLUMINA, SOLID, LS454, HELICOS and PACBIO
COV_OUT_DIR_NAME="coverage" #coverage file output directory name. Just a name, it will be placed under the sample directory
DUP_FREE_BAM_DIR_NAME="dup_free_bam" # folder to ouput bam file without any GATK process, this can be input of freebayes
GATK_OUT_DIR_NAME="gatk" #gatk tuned alignment file output directory name
WRKGDIR_NAME="_no_QC" # tail of temp working directory under scratch dir, in case same sample running at the same using the same temp directory
JAVA_TMP_DIR_NAME="_no_QC" #tail of temp java working directory under scratch dir
#**************************/

#rm -r ${SCRATCH_DIR}
#mkdir ${SCRATCH_DIR}

###############################

## Submitting jobs ##
INALL=( "${SAMPLE_PATH}"/Sample_* )
for SampleFolder in "${INALL[@]}"; do
	SAMPLE_ID="${SampleFolder##*/Sample_}"
	#if [[ $SAMPLE_ID == NG-7460_SA_D64614 ||  $SAMPLE_ID == NG-7460_VH_D197700 ]]; then
	if [[ $SAMPLE_ID == F99248 ]]; then
	       	JOB_ID1="Yaobo_Exome_${SAMPLE_ID}_1st_stage"
		arr=("${SAMPLE_ID}" "${SAMPLE_PATH}" "${SCRIPT_PATH}" "$REFDIR" "${SCRATCH_DIR}" "$TARGETS" "${SEQ_PLATFORM}" "${Library_ID}" "${COV_OUT_DIR_NAME}" "${GATK_OUT_DIR_NAME}" "${DUP_FREE_BAM_DIR_NAME}" "${WRKGDIR_NAME}" "${JAVA_TMP_DIR_NAME}")
	        echo "${SAMPLE_ID}"
		#echo qsub -N "${JOB_ID1}" ${SCRIPT_PATH}/redupFastuniq_31Jul2014.sh "${arr[@]}"
		qsub -N "${JOB_ID1}" ${SCRIPT_PATH}/redupFastuniq_31Jul2014.sh "${arr[@]}"
		#sh ${SCRIPT_PATH}/redupFastuniq_31Jul2014.sh "${arr[@]}"
	fi
done

