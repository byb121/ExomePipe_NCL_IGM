#! /bin/bash

################## Paras need to be adjusted for diff samples ##############
INALL=(D05685 D08064 D10020 D100563 D11016 D114041 D115022 D122641 D122836,D99945,D62925 D129266 D134130 D13415 D146989 D163417 D164942 D16618,D050227,D62018 D169918 D170484 D171748 D171987 D174381 D180850 D36421 D83732 D86968 D88892 D9646 F00603 F02209 F03450 F03513 F03566 F03665 F05163 F05460 F05506 F06341 F07191 F07230 F08196 F0852 F08560 F08744 F09507,F09508 F09591 F10149 F11151 F13349 F98212 F98407 F98468 F99129,F99359)

#INALL=(NG-7460_AM_D173655 NG-7460_FMC_D189345 NG-7460_HT NG-7460_NMC_D173654 NG-7460_PM_D173656 NG-7460_RB NG-7460_RKB_D131490 NG-7460_SA_D64614 NG-7460_TP_D66204 NG-7460_VH_D197700)

#INALL=(F07230 D163417 F13349 D115022 D134130 D174381 D146989)

#INALL=(F98407 D129266 F02209 F03513 F98468 F03450)
#INALL=(D180850)

SAMPLEDIR="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/MergeV" # the folder contain the VCF file which contain all of the samples
INPUT_VCF="${SAMPLEDIR}/merged_Sophie_A1969_20140808.vcf"
AnnovarDIR="/sharedlustre/IGM/annovar_2014jul14"

ANOVAR="N" #use to indicate if Anovar output has beed produced, Y means yes, then After_Anovar script will be used
PERL_SCRIPT_DIR=$(pwd) # change if the annotation perl script is stored somewhere else
PERL_SCRIPT="${PERL_SCRIPT_DIR}/VCF_2_annotated_excel_20141010.pl"
PERL_SCRIPT_AF_ANOVAR="${PERL_SCRIPT_DIR}/VCF_2_annotated_excel_20141010_After_Annovared.pl" #note this name has to contain 'After_Annovared', otherwise GATK will start to select samples

GENO_CALL_FLAG="Yes"

BATCH_MAF="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/MergeV/merged_Sophie_A1969_20140808_batch_MAF.txt"
CNV="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/AnnotateV/CNV_all_A1969.txt"
InterestedGENES="/home/a5907529/WORKING_DATA/Sophie_A1969/scripts/AnnotateV/Genes_PID_01DEC2014_YAOBO.txt"

VCF_TAIL="_selected.vcf"
FILTERED_TAIL="_filtered.xlsx" # output excel file name tail of filtered variants, has to be ended with xlsx
EVERYTHING_TAIL="_Everything.xlsx" # output excel file name tail of all variants, has to be ended with xlsx
CUSTOM_VCF_TAIL="_selected.custom.vcf"

###############################

for SAMPLE in "${INALL[@]}"
do
	STRING=${SAMPLE//,/_}
	SAMPLE_SELECTED_VCF="${SAMPLEDIR}/${STRING}${VCF_TAIL}"
	OUTPUT="${SAMPLEDIR}/${STRING}${FILTERED_TAIL}"
	OUTPUT_ALL="${SAMPLEDIR}/${STRING}${EVERYTHING_TAIL}"
	OUTPUT_VCF="${SAMPLEDIR}/${STRING}${CUSTOM_VCF_TAIL}"

	if [ $ANOVAR == "Y" ]
	then
		SCRIPT="${PERL_SCRIPT_AF_ANOVAR}"
	else
		SCRIPT="${PERL_SCRIPT}"
	fi
	arr=("$SCRIPT" "${SAMPLE_SELECTED_VCF}" "$OUTPUT" "$OUTPUT_ALL" "$INPUT_VCF" "$SAMPLE" "$CNV" "$InterestedGENES" "$AnnovarDIR" "${PERL_SCRIPT_DIR}" "${BATCH_MAF}" "${OUTPUT_VCF}")
	qsub ${PERL_SCRIPT_DIR}/VCF_filters_20140811.sh "${arr[@]}"

 	#sh ${PERL_SCRIPT_DIR}/VCF_filters_20140811.sh "${arr[@]}"
done

