#! /bin/bash
#$ -cwd -V -j y
#$ -l h_vmem=16g
# #$ -M yaobo.xu@ncl.ac.uk
# #$ -m e

set -e
res1=$(date +%s.%N) # use to calculate whole run time of the job

#/********************** on lampredi2
#**********************/


#/********************** on fmscluster
SCRATCH_DIR=$TMPDIR
#**********************/

REF_FILE=$1
LOOKING_REGION=$2 #eg: chr1
INPUT=$3 #Bam list file
OUTPUT_FILE=$4
JOB_ID=$5

echo "Freebayes is started to call variants on region: $LOOKING_REGION"

freebayes -v ${OUTPUT_FILE} -f $REF_FILE -L ${INPUT} --use-mapping-quality -C 5 --ploidy 2 -r ${LOOKING_REGION} \
--min-alternate-count 5 \
--min-mapping-quality 30 \
--min-base-quality 30 \
--min-coverage 5 \
--genotype-qualities

echo "Done."

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
