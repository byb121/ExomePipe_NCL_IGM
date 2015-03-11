#! /bin/bash
#$ -cwd -V -j y
# #$ -l h_vmem=30g
#$ -M yaobo.xu@ncl.ac.uk
#$ -m e

set -e
res1=$(date +%s.%N) # use to calculate whole run time of the job

#/********************** on lampredi2
#module load apps/gatk/3.1.1/noarch
#module load apps/gatkqueue/3.1.1/noarch
#**********************/


#/********************** on fmscluster
module add apps/perl/5.18.2
SCRATCH_DIR=$TMPDIR
#**********************/

SCRIPT="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/MergeV/VCF_combiner_HC_FREE.pl"
REF_FAI="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/MergeV/human_g1k_v37_decoy.fasta.fai"
FREE="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/Freebayes/Free_out_combined.vcf"
HC="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/GenotypeGVCFs/New_HC_out_recali_SNP_INDEL.vcf"
OUTPUT="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/MergeV/merged_Sophie_A1969_Rachel_new_20141023.vcf"
BATCH_MAF="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/MergeV/merged_Sophie_A1969_Rachel_new_20141023_batch_MAF.txt"

perl $SCRIPT --fai $REF_FAI --HC $HC --FREE $FREE --output $OUTPUT --MAF ${BATCH_MAF}

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
