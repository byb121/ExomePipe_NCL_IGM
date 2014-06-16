#! /bin/bash
#$ -cwd
#$ -j y
#$ -V

#/***************** on fmsclustergw
#module load apps/python27/2.7.4
#module load apps/perl/5.18.2
#*****************/

#/***************** on lampredi2
module load apps/python/2.7.3/gcc-4.4.6
module load apps/perl/5.16.1/gcc-4.4.6
#*****************/

FASTQ_FILE=$1
OUTPUT_DIR=$2
fastqc --nogroup --noextract -o ${OUTPUT_DIR} ${FASTQ_FILE}
echo "done."
