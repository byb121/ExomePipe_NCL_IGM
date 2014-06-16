#! /bin/bash
#$ -cwd
#$ -j y

#/******************* on Lampredi2
module load apps/samtools/0.1.18/gcc-4.4.6
module load apps/perl/5.16.1/gcc-4.4.6
module load apps/bedtools/2.17.0/gcc-4.4.6
#*******************/

#/******************* on fmsclustergw
#module load apps/samtools/0.1.19
#module load apps/bedtools/2.19.0
#module load apps/perl/5.18.2
#*******************/

BAM=$1
TARGETS=$2
OUTPUT=$3

coverageBed -abam $BAM -b $TARGETS -hist -split  > $OUTPUT

