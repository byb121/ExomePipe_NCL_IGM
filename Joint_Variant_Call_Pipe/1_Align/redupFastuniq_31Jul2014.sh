#! /bin/bash
#$ -cwd -V
#$ -j y
#$ -l h_vmem=70G # fastuniq require a lot of RAM
#$ -M yaobo.xu@ncl.ac.uk
#$ -m e

set -e
res1=$(date +%s.%N) # use to calculate whole run time of the job
#cwd=`pwd`; export cwd

#/************************************
# note: fastq file patterns are used in the perl script detectSampleLanes.pl and READ_FILE1 & READ_FILE2
# they need to be adjusted accordingly
#************************************/

echo $'\n'"["`date`"]: Job stared."

SAMPLE_ID=$1
SAMPLE_PATH=$2
SCRIPTS_DIR=$3
REF_DIR=$4
SCRATCH_DIR=$5
TARGETS=$6
SEQ_PLATFORM=$7 # Valid values are: ILLUMINA, SOLID, LS454, HELICOS and PACBIO
Library_ID=$8 #Better to use {PI name}.{Batch Number} to identify the library
COV_DIR_NAME=$9
GATK_OUT_DIR_NAME=${10}
DUP_FREE_BAM_DIR_NAME=${11}
WRKGDIR_NAME=${12}
JAVA_TMP_DIR_NAME=${13}

#/************************ fmscluster add modules
# need to add line here to load fastuniq module
#GATK 3.1.1 /home/a5907529/WORKING_DATA/QC_GOOD/GATK3
#*************************/

#/************************ lampredi2 add modules
#module load apps/fastuniq/1.1/gcc-4.4.6
#*************************/

INDIR="${SAMPLE_PATH}/Sample_${SAMPLE_ID}"

LANES_STRING=`perl ${SCRIPTS_DIR}/detectSampleLanes.pl $SAMPLE_PATH $SAMPLE_ID`
LANES=($LANES_STRING)

if [ ${#LANES[@]} -lt 1 ]
then
	exit "Error: No lane info detected, check file name patterns and paths"
fi

BAM_FILE_LIST=""
echo $'\n'"started to loop through lanes $LANES_STRING.."
echo "----------------------------- hua li li de fen jie xian -----------------------------"
echo "----------------------------- hua li li de fen jie xian -----------------------------"
for LANE in "${LANES[@]}"
do
	echo $'\n'"---------------- hua li de fen jie xian ----------------"
	echo "-------------------------------------------------"
	echo "remove duplicated reads from lane $LANE for sample ${SAMPLE_ID}"
	echo "-------------------------------------------------"
	
	#/********************** if QC ed with the Trimming.sh
	#READ_FILE1="$INDIR/${SAMPLE_ID}_L00${LANE}*val_1.fq"
        #READ_FILE2="$INDIR/${SAMPLE_ID}_L00${LANE}*val_2.fq"
	#**********************/	

	READ_FILE1=( $INDIR/${SAMPLE_ID}_*_L00${LANE}_R1_001.fastq* )
	READ_FILE2=( $INDIR/${SAMPLE_ID}_*_L00${LANE}_R2_001.fastq* )
	
	READ_FILE1=${READ_FILE1[0]}
	READ_FILE1=${READ_FILE1%'.gz'}
	READ_FILE2=${READ_FILE2[0]}
	READ_FILE2=${READ_FILE2%'.gz'}
	
	echo $READ_FILE1; echo $READ_FILE2
	# can handle zipped file now, will do unzip automatically
	if [ -f "$READ_FILE1.gz" ]
	then
		echo "find gzipped fastq file $READ_FILE1.gz, unzip now.."
		gunzip -c "$READ_FILE1.gz" > "$READ_FILE1"
	else
		if [ -f $READ_FILE1 ]
		then
			echo "find file $READ_FILE1"
		else
			exit "fasq file name error, can not find file $READ_FILE1"
		fi
	fi

	if [ -f "$READ_FILE2.gz" ]
	then
		echo "find gzipped fastq file $READ_FILE2.gz, unzip now.."
                gunzip -c "$READ_FILE2.gz" > "$READ_FILE2"
	else
		if [ -f $READ_FILE2 ]
		then
                        echo "find file $READ_FILE2"
                else
                        exit "fasq file name error, can not find file $READ_FILE2"
                fi
        fi
	
	READ_FILE1_nodup="$INDIR/${SAMPLE_ID}_L00${LANE}_R1_001.nodup.fastq"
        READ_FILE2_nodup="$INDIR/${SAMPLE_ID}_L00${LANE}_R2_001.nodup.fastq"
	FASTQ_LIST_FILE="$INDIR/${SAMPLE_ID}_L00${LANE}.temp.qlist"
	# remvoe dups with FastUniq
	echo -e "$READ_FILE1\n$READ_FILE2" > $FASTQ_LIST_FILE
	echo fastuniq -i "${FASTQ_LIST_FILE}" -t q -o "${READ_FILE1_nodup}" -p "${READ_FILE2_nodup}"
	fastuniq -i "${FASTQ_LIST_FILE}" -t q -o "${READ_FILE1_nodup}" -p "${READ_FILE2_nodup}"
	
	echo rm $READ_FILE1
	rm ${READ_FILE1}
	echo rm $READ_FILE2
	rm ${READ_FILE2}
	
	echo "-------------------------------------------------"
	echo "reads from lane $LANE for sample ${SAMPLE_ID} is duplicated free now"
	echo "-------------------------------------------------"

done

echo $'\n'"["`date`"]: The fastuniq job is DONE!!"
echo "["`date`"]: Submitting following job"
JOB_ID="Yaobo_Exome_${SAMPLE_ID}_2nd_stage"
arr=("${SAMPLE_ID}" "${SAMPLE_PATH}" "${SCRIPTS_DIR}" "${REF_DIR}" "${SCRATCH_DIR}" "${TARGETS}" "${SEQ_PLATFORM}" "${Library_ID}" "${COV_DIR_NAME}" "${GATK_OUT_DIR_NAME}" "${DUP_FREE_BAM_DIR_NAME}" "${WRKGDIR_NAME}" "${JAVA_TMP_DIR_NAME}")
#qsub -N "${JOB_ID}" ${SCRIPTS_DIR}/map_recali_perLane_recali_perSample_covOnTargets_GVCF_31Jul2014.sh "${arr[@]}"
qsub -N "${JOB_ID}" ${SCRIPTS_DIR}/map_redup.sh "${arr[@]}"

# runtime calculation
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)
printf "fastuniq runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds
echo "exit status $?"

