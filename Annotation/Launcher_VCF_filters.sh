#! /bin/bash

################## Paras need to be adjusted for diff samples ##############
INALL=(D05685 D08064 D10020 D100563 D11016 D114041 D115022 D122641 D122836 D13415 D146989 D163417 D164942 D169918 D170484 D171748 D171987 D174381 D36421 D62925 D83732 D86968 D88892 D9646 D99945 D16618_D050227_D62018 D122836_D99945_D62925)
INPUT_VCF="$SAMPLEDIR/New_HaplotypeCaller_recali_SNP_INDEL.vcf"
REF_FILE="/users/data/GATK/bundle2.8/hg19/ucsc.hg19.YAOBO.fasta"
AnnovarDIR="/users/a5907529/data/Files_HG/vcf_annotation_november2013"

ANOVAR="N" #use to indicate if Anovar output has beed produced, Y means yes, then After_Anovar script will be used
PERL_SCRIPT_DIR=$(`pwd`)
PERL_SCRIPT="${PERL_SCRIPT_DIR}/VCF_2_annotated_xlsx_20140501.pl"
PERL_SCRIPT_AF_ANOVAR="${PERL_SCRIPT_DIR}/VCF_2_annotated_xlsx_20140501_After_Annovared.pl" #note this name has to contain 'After_Annovared', otherwise GATK will start to select samples

GENO_CALL_FLAG="Yes"

CNV="/users/a5907529/lustre/Sophie_unsovled_cases/CNVs/batch_20130528/batch_20130528_annotated_CNVs.txt"
InterestedGENES="/users/a5907529/lustre/Sophie_unsovled_cases/Genes_PID_01MAY2014_YAOBO.txt"

VCF_TAIL="_selected.vcf"
FILTERED_TAIL="_filtered.xlsx" # output excel file name tail of filtered variants, has to be ended with xlsx
EVERYTHING_TAIL="_Everything.xlsx" # output excel file name tail of all variants, has to be ended with xlsx


###############################

for SAMPLE in "${INALL[@]}"
do
	IFS='_' read -a EXPs <<< "$SAMPLE"  #### split string
	STRING=""
	for EXP in "${EXPs[@]}"
	do
		STRING="${STRING}-se "
    		STRING="${STRING}$EXP "
	done
	echo $STRING
	SAMPLE_SELECTED_VCF="${SAMPLEDIR}/${SAMPLE}${VCF_TAIL}"
	OUTPUT="${SAMPLEDIR}/${SAMPLE}${FILTERED_TAIL}"
	OUTPUT_ALL="${SAMPLEDIR}/${SAMPLE}${EVERYTHING_TAIL}"
	if [ $ANOVAR == "Y" ]
	then
		SCRIPT="${PERL_SCRIPT_AF_ANOVAR}"
	else
		SCRIPT="${PERL_SCRIPT}"
	fi
	#echo "${SAMPLEDIR}/VCF_filters.sh $SCRIPT ${SAMPLE_SELECTED_VCF} $OUTPUT $OUTPUT_ALL $INPUT_VCF $STRING $REF_FILE $CNV $InterestedGENES"
	arr=("$SCRIPT" "${SAMPLE_SELECTED_VCF}" "$OUTPUT" "$OUTPUT_ALL" "$INPUT_VCF" "$STRING" "$REF_FILE" "$CNV" "$InterestedGENES" "$AnnovarDIR")
	qsub ${SAMPLEDIR}/VCF_filters.sh "${arr[@]}"
done
