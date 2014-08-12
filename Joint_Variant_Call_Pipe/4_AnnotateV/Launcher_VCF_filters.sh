#! /bin/bash

################## Paras need to be adjusted for diff samples ##############
#INALL=(D05685 D170484 F05163 D08064 D171748 F05460 D10020 D171987 F05506 D100563 D174381 F06341 D11016 D180850 F07191 D114041 D36421 F08196 D115022 D62925 F0852 D122641 D83732 F08560 D122836,D99945,D62925 D86968 F08744 D122836 D88892 F09507,F09508 D129266 D9646 F09591 D134130 D99945 F10149 D13415 F00603 F11151 D146989 F02209 F13349 D163417 F03450 F98212 D164942 F03513 F98407 D16618,D050227,D62018 F03566 F98468 D169918 F03665 F99129,F99359 F07230 NG-7460_HT NG-7460_RB)

#INALL=(F09507,F09508)
#INALL=(D16618,D050227,D62018)
INALL=(F07230)
SAMPLEDIR="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/MergeV" # the folder contain the VCF file which contain all of the samples
INPUT_VCF="${SAMPLEDIR}/merged_Sophie_A1969_20140808.vcf"
AnnovarDIR="/sharedlustre/IGM/annovar_2014jul14"

ANOVAR="N" #use to indicate if Anovar output has beed produced, Y means yes, then After_Anovar script will be used
PERL_SCRIPT_DIR=$(pwd) # change if the annotation perl script is stored somewhere else
PERL_SCRIPT="${PERL_SCRIPT_DIR}/VCF_2_annotated_excel_20140811.pl"
PERL_SCRIPT_AF_ANOVAR="${PERL_SCRIPT_DIR}/VCF_2_annotated_excel_20140811_After_Annovared.pl" #note this name has to contain 'After_Annovared', otherwise GATK will start to select samples

GENO_CALL_FLAG="Yes"

CNV="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/AnnotateV/CNV_all_A1969.txt"
InterestedGENES="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/AnnotateV/Genes_PID_01MAY2014_YAOBO.txt"

VCF_TAIL="_selected.vcf"
FILTERED_TAIL="_filtered.xlsx" # output excel file name tail of filtered variants, has to be ended with xlsx
EVERYTHING_TAIL="_Everything.xlsx" # output excel file name tail of all variants, has to be ended with xlsx

###############################

for SAMPLE in "${INALL[@]}"
do
	STRING=${SAMPLE//,/_}
	SAMPLE_SELECTED_VCF="${SAMPLEDIR}/${STRING}${VCF_TAIL}"
	OUTPUT="${SAMPLEDIR}/${STRING}${FILTERED_TAIL}"
	OUTPUT_ALL="${SAMPLEDIR}/${STRING}${EVERYTHING_TAIL}"
	if [ $ANOVAR == "Y" ]
	then
		SCRIPT="${PERL_SCRIPT_AF_ANOVAR}"
	else
		SCRIPT="${PERL_SCRIPT}"
	fi
	#echo "${SAMPLEDIR}/VCF_filters.sh $SCRIPT ${SAMPLE_SELECTED_VCF} $OUTPUT $OUTPUT_ALL $INPUT_VCF $SAMPLE $REF_FILE $CNV $InterestedGENES"
	arr=("$SCRIPT" "${SAMPLE_SELECTED_VCF}" "$OUTPUT" "$OUTPUT_ALL" "$INPUT_VCF" "$SAMPLE" "$CNV" "$InterestedGENES" "$AnnovarDIR" "${PERL_SCRIPT_DIR}")
	#qsub ${PERL_SCRIPT_DIR}/VCF_filters_20140811.sh "${arr[@]}"
	sh ${PERL_SCRIPT_DIR}/VCF_filters_20140811.sh "${arr[@]}"
done

