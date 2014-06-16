#! /bin/bash
#$ -cwd -V
#$ -l h_vmem=16g

#/********************** on lampredi2
module load apps/gatk/3.1.1/noarch
module load apps/gatkqueue/3.1.1/noarch
#**********************/


#/********************** on fmscluster
#module add apps/java
#GATKQUEUEDIR="/home/a5907529/WORKING_DATA/Project_Fulcrum/GATK3"
#export _JAVA_OPTIONS="-XX:-UseLargePages -Xms4096m"
#**********************/

WKDIR=`pwd`
BAM_LIST_FILE="$WKDIR/BamFile.list"
SCATTER_NUMBER="100"
OUTPUT_FILE="$WKDIR/HaplotyperCaller.vcf"
ACTIVE_REGIONS="$WKDIR/activeRegions.txt"
REF_FILE="/sharedlustre/IGM/bundle2.8/hg19/ucsc.hg19.fasta"
TEMP_DIR="$WKDIR/HaplotypeCaller_tmp"
mkdir $TEMP_DIR

java -Djava.io.tmpdir=$TEMP_DIR -Xmx12g -jar $GATKQUEUEDIR/Queue-ben.jar -S $WKDIR/HaplotypeCaller.scala -I $BAM_LIST_FILE -R $REF_FILE -sg $SCATTER_NUMBER -V_out $OUTPUT_FILE -jobRunner GridEngine -retry 2 -run -l DEBUG -activeRegionsOut $ACTIVE_REGIONS

#rm -r $TEMP_DIR
