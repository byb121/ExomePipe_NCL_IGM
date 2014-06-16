#! /bin/bash
#$ -cwd
#$ -j y
#$ -V
#$ -l h_vmem=12g

#/********************** on fmsclustergw
module add apps/perl/5.18.2
module add apps/python27/2.7.4
module add apps/java
export _JAVA_OPTIONS="-XX:-UseLargePages -Xms4096m"
#**********************/

#/********************** on lampredi2
#module load apps/python/2.7.3/gcc-4.4.6
#module load apps/perl/5.16.1/gcc-4.4.6
#**********************/

READ_1_FILE=$1
READ_2_FILE=$2
SAMPLE_DIR=$3
SAMPLE_NAME=$4
FASTQC_OUTPUT=$5
SCRIPT_DIR_PERL=$6

TRIMMED_FASTQ_1="${SAMPLE_DIR}/${SAMPLE_NAME}_1.trimmed.fastq"
Ns_RMOVED_FASTQ_1="${SAMPLE_DIR}/${SAMPLE_NAME}_1.trimmed.fastq.polyN_removed.txt"

TRIMMED_FASTQ_2="${SAMPLE_DIR}/${SAMPLE_NAME}_2.trimmed.fastq"
Ns_RMOVED_FASTQ_2="${SAMPLE_DIR}/${SAMPLE_NAME}_2.trimmed.fastq.polyN_removed.txt"

#seqtk to trim off begaining 15bp and ending 3bp
echo $'\n'"seqtk trimfq -e 3 ${READ_1_FILE} > ${TRIMMED_FASTQ_1}"
seqtk trimfq -e 3 ${READ_1_FILE} > ${TRIMMED_FASTQ_1}
echo $'\n'"seqtk trimfq -e 3 ${READ_2_FILE} > ${TRIMMED_FASTQ_2}"
seqtk trimfq -e 3 ${READ_2_FILE} > ${TRIMMED_FASTQ_2}


#perl remove Ns tails but also remove read pairs that is shorter than 20 after trimming Ns
echo $'\n'"perl ${SCRIPT_DIR_PERL}/removeNtailsFromFastQ.pl --1 ${TRIMMED_FASTQ_1} --2 ${TRIMMED_FASTQ_2}"
perl ${SCRIPT_DIR_PERL}/removeNtailsFromFastQ.pl --1 ${TRIMMED_FASTQ_1} --2 ${TRIMMED_FASTQ_2}

#rm ${TRIMMED_FASTQ_1}
#rm ${TRIMMED_FASTQ_2}

#trim_galore to remove adaptors and low quality tails and perform fastqc on the output
echo $'\n'"trim_galore -q 20 --paired --length 20 --stringency 5 -o ${SAMPLE_DIR} --fastqc_args \"--nogroup --noextract --outdir ${FASTQC_OUTPUT}\" ${Ns_RMOVED_FASTQ_1} ${Ns_RMOVED_FASTQ_2}"
trim_galore -q 20 --paired --length 20 --stringency 5 -o ${SAMPLE_DIR} --fastqc_args "--nogroup --noextract --outdir ${FASTQC_OUTPUT}" ${Ns_RMOVED_FASTQ_1} ${Ns_RMOVED_FASTQ_2}

#rm ${Ns_RMOVED_FASTQ_1}
#rm ${Ns_RMOVED_FASTQ_2}
