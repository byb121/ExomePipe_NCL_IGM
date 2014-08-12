#! /bin/bash
#$ -cwd -V -j y
#$ -l h_vmem=16g
# #$ -M yaobo.xu@ncl.ac.uk
# #$ -m e

set -e
res1=$(date +%s.%N) # use to calculate whole run time of the job

#/********************** on lampredi2
#module load apps/gatk/3.1.1/noarch
#module load apps/gatkqueue/3.1.1/noarch
#**********************/


#/********************** on fmscluster
module add apps/java/jre-1.7.0_55
module add apps/gatk/3.2-protected
SCRATCH_DIR=$TMPDIR
export _JAVA_OPTIONS="-XX:-UseLargePages"
GATKDIR="$GATK_ROOT"
#**********************/

REF_FILE=$1
echo "Ref file: $REF_FILE"
LOOKING_REGION=$2 #eg: chr1
echo "Region: $LOOKING_REGION"
INPUT=$3
echo "Input GVCF para: $INPUT"
OUTPUT_FILE=$4
echo "Ouptut: $OUTPUT_FILE"
JOB_ID=$5
echo "Job id: $JOB_ID"
JAVA_TMP_DIR_NAME=$6
echo "Java temp dir name: $JAVA_TMP_DIR_NAME"
JAVA_TEMP_DIR="$TMPDIR/GenotypeGVCFs_${JOB_ID}${JAVA_TMP_DIR_NAME}"

if [ ! -d $JAVA_TEMP_DIR ]; then
	echo mkdir $JAVA_TEMP_DIR
	mkdir $JAVA_TEMP_DIR
else
	echo "Warning: Temprorary directory $JAVA_TEMP_DIR exists!"
fi

java -Djava.io.tmpdir=$JAVA_TEMP_DIR -Xmx12g -jar $GATKDIR/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $REF_FILE -o $OUTPUT_FILE -A InbreedingCoeff -A FisherStrand -A QualByDepth -A ChromosomeCounts -A GenotypeSummaries -L $LOOKING_REGION \
$INPUT

#rm -r $TEMP_DIR

echo $'\n'"["`date`"]: The whole job is DONE!!"

# runtime calculation
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)
echo "----------------------------- hua li li de fen jie xian -----------------------------"
printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds
echo "exit status $?"
