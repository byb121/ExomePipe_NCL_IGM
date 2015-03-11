#! /bin/bash
#$ -cwd -V
#$ -j y
# #$ -l h_vmem=70G # fastuniq require a lot of RAM
# #$ -M yaobo.xu@ncl.ac.uk
# #$ -m e

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
JOB_ID="Yaobo_Exome_${SAMPLE_ID}_2nd_stage"
arr=("${SAMPLE_ID}" "${SAMPLE_PATH}" "${SCRIPTS_DIR}" "${REF_DIR}" "${SCRATCH_DIR}" "${TARGETS}" "${SEQ_PLATFORM}" "${Library_ID}" "${COV_DIR_NAME}" "${GATK_OUT_DIR_NAME}" "${DUP_FREE_BAM_DIR_NAME}" "${WRKGDIR_NAME}" "${JAVA_TMP_DIR_NAME}")

qsub -N "${JOB_ID}" ${SCRIPTS_DIR}/redo_coverage.sh "${arr[@]}"

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

