#! /bin/bash
#$ -cwd
#$ -j y
#$ -l h_vmem=8g

# Load Modules
#/********************on Lampredi
#module load apps/gatk/3.1.1/noarch
#*********************/


#/********************on fmscluster
GATKDIR="/home/a5907529/WORKING_DATA/Project_Fulcrum/GATK3"
module add apps/java
export _JAVA_OPTIONS="-XX:-UseLargePages -Xms4096m"
#********************/

INPUT_VCF=$1
OUTPUT_PATH=$2
SAMPLE=$3
REF=$4

OUTPUT_VCF="$OUTPUT_PATH/${SAMPLE}.selected.vcf"

java -Xmx4g -jar $GATKDIR/GenomeAnalysisTK.jar -R $REF \
-T SelectVariants \
--variant $INPUT_VCF \
-env \
-ef \
-o $OUTPUT_VCF \
-sn $SAMPLE

