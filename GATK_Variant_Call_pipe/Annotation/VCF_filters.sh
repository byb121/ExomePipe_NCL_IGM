#! /bin/bash
#$ -cwd
#$ -j y
#$ -V

#/*********************** on Lampredi2
module load apps/gatk/3.1.1/noarch
module load apps/perl/5.16.1/gcc-4.4.6
#***********************/

SCRIPT=$1
SAMPLE_SELECTED_VCF=$2
OUTPUT=$3
OUTPUT_Everthing=$4
INPUT_VCF=$5
EXP=$6
REF_FILE=$7
CNV=$8
InterestedGENES=$9
AnnovarDIR=${10}

if [[ $SCRIPT != *After_Annovared* ]]
then
	java -Xmx4g -jar $GATKDIR/GenomeAnalysisTK.jar -T SelectVariants -R $REF_FILE \
	--variant $INPUT_VCF \
	-env -ef \
	-o ${SAMPLE_SELECTED_VCF} \
	$EXP
fi

perl $SCRIPT --vcf ${SAMPLE_SELECTED_VCF} \
--out ${OUTPUT} \
--outAll ${OUTPUT_Everthing} \
--add_genotypeCall_flags Yes \
--CNV $CNV \
--InterestedGenes $InterestedGENES \
--AnnovarDIR $AnnovarDIR
