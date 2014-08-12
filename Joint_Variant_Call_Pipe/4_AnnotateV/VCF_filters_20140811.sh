#! /bin/bash
#$ -cwd
#$ -l h_vmem=30g
#$ -j y -V
#$ -M yaobo.xu@ncl.ac.uk
#$ -m e

#/***********************
module load apps/perl
#***********************/

SCRIPT=$1
SAMPLE_SELECTED_VCF=$2
OUTPUT=$3
OUTPUT_Everthing=$4
INPUT_VCF=$5
SAMPLE=$6
CNV=$7
InterestedGENES=$8
AnnovarDIR=$9
PERL_SCRIPT_DIR=${10}

if [[ $SCRIPT != *After_Annovared* ]]
then
	perl $PERL_SCRIPT_DIR/FilterVCF4samples.pl \
	--vcf ${INPUT_VCF} \
	--samples $SAMPLE \
	--output ${SAMPLE_SELECTED_VCF}
fi

perl $SCRIPT --vcf ${SAMPLE_SELECTED_VCF} \
--out ${OUTPUT} \
--outAll ${OUTPUT_Everthing} \
--add_genotypeCall_flags Yes \
--CNV $CNV \
--InterestedGenes $InterestedGENES \
--AnnovarDIR $AnnovarDIR

