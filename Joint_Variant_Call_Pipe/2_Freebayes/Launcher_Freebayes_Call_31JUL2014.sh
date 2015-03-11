#! /bin/bash


## Script to launch Exome Pipeline (Picard - GATK bam recalibration - GATK call variants & varscan & dindel & bedtools for coverage) for multiple samples
## INALL is an array of directory names, with one directory per sample
## SAMPLE_PATH is the pathway to the sample directories ${SAMPLE_PATH}/${INALL[$SAMPLE_NUMBER]}


#/************************* Paras need to be adjusted for diff samples
BAM_LIST_FILE="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/Freebayes/Bam.list"
SCRIPT_PATH="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/Freebayes"
REF_FILE="/sharedlustre/IGM/GATK/bundle2.8/b37/human_g1k_v37_decoy.fasta"
OUT_DIR="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/Freebayes/Free_out" #gatk tuned alignment file output directory name

if [ ! -d $OUT_DIR ]; then
	mkdir $OUT_DIR
else
	echo "Directory: $OUT_DIR exists, will not overwrite"
fi

#**************************/

#REF_FILE=$1
#LOOKING_REGION=$2 #eg: chr1
#INPUT=$3 #Bam list file
#OUTPUT_FILE=$4
#JOB_ID=$5

###############################
REGIONS=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y' 'MT')

#REGIONS=('1:0-249250621' '2:0-243199373' '3:0-198022430' '4:0-191154276' '5:0-180915260' '6:0-171115067' '7:0-159138663' '8:0-146364022' '9:0-14121343' '10:0-135534747' '11:0-135006516' '12:0-133851895' '13:0-115169878' '14:0-107349540' '15:0-102531392' '16:0-90354753' '17:0-81195210' '18:0-78077248' '19:0-59128983' '20:0-63025520' '21:0-48129895' '22:0-51304566' 'X:0-155270560' 'Y:0-59373566' 'MT:0-16569')

#REGIONS=('chrM' 'chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX' 'chrY')

## Submitting jobs ##
for LOOKING_REGION in "${REGIONS[@]}"; do
	#if [[ $LOOKING_REGION == '1:0-249250621' ]]; then
	       	JOB_ID1="Yaobo_Exome_${LOOKING_REGION/:/_}_Freebayes"
		OUTPUT_FILE="$OUT_DIR/$JOB_ID1.vcf"
		echo $JOB_ID1
		arr=("${REF_FILE}" "${LOOKING_REGION}" "${BAM_LIST_FILE}" "${OUTPUT_FILE}" "${JOB_ID1}")
		#echo qsub -N "${JOB_ID1}" ${SCRIPT_PATH}/Freebayes_Call.sh "${arr[@]}"
		qsub -N "${JOB_ID1}" ${SCRIPT_PATH}/Freebayes_Call.sh "${arr[@]}"
		#sh ${SCRIPT_PATH}/GenotypeGVCFs.sh "${arr[@]}"
	#fi
done

