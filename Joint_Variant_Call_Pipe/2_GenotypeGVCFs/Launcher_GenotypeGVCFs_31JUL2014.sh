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
GVCF_LIST_FILE="/home/a5907529/WORKING_DATA/Linda_20141024/scripts/2_GenotypeGVCFs/GVCF.list"
SCRIPT_PATH="/home/a5907529/WORKING_DATA/Linda_20141024/scripts/2_GenotypeGVCFs"
REF_FILE="/sharedlustre/IGM/GATK/bundle2.8/b37/human_g1k_v37_decoy.fasta"
OUT_DIR="/home/a5907529/WORKING_DATA/Linda_20141024/scripts/2_GenotypeGVCFs/HC_out" #gatk tuned alignment file output directory name
JAVA_TMP_DIR_NAME="_no_QC" #tail of temp java working directory under scratch dir

if [ ! -d $OUT_DIR ]; then
	mkdir $OUT_DIR
else
	echo "Directory: $OUT_DIR exists, will not overwrite"
fi

#**************************/

###############################
REGIONS=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y" "MT")

#process input GVCFs to -V gvcf.1 -V gvdf.2 ......
readarray GVCFs < $GVCF_LIST_FILE
INPUT=""
for GVCF in "${GVCFs[@]}"; do
	GVCF=(`echo -n $GVCF`)
	INPUT="$INPUT -V $GVCF"
done

## Submitting jobs ##
for LOOKING_REGION in "${REGIONS[@]}"; do
	#if [[ $LOOKING_REGION == 1 ]]; then
	       	JOB_ID1="Yaobo_Exome_${LOOKING_REGION}_HC_GenotypeGVCFs"
		OUTPUT_FILE="$OUT_DIR/$JOB_ID1.vcf"
		echo $JOB_ID1
		arr=("${REF_FILE}" "${LOOKING_REGION}" "${INPUT}" "${OUTPUT_FILE}" "${JOB_ID1}" "$JAVA_TMP_DIR_NAME")
		#echo qsub -N "${JOB_ID1}" ${SCRIPT_PATH}/GenotypeGVCFs.sh "${arr[@]}"
		qsub -N "${JOB_ID1}" ${SCRIPT_PATH}/GenotypeGVCFs.sh "${arr[@]}"
		#sh ${SCRIPT_PATH}/GenotypeGVCFs.sh "${arr[@]}"
	#fi
done

