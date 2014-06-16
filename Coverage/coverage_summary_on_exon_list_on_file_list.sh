#! /bin/bash
#$ -cwd
#$ -j y

#/******************* on Lampredi2
module load apps/perl/5.16.1/gcc-4.4.6
#*******************/

#/******************* on fmscluster
#module load apps/perl/5.18.2
#*******************/

SAMPLE_PATH=$1
SAMPLE_ID_LIST=$2
OUTPUT_prefix=$3
COV_FILE_NAME=$4
SCRIPT_PATH=$5

perl ${SCRIPT_PATH}/coverage_summary_on_exon_list_on_file_list.pl \
--Output_prefix $OUTPUT_prefix \
--SamplePath $SAMPLE_PATH \
--SampleNames $SAMPLE_ID_LIST \
--CovFileName $COV_FILE_NAME
