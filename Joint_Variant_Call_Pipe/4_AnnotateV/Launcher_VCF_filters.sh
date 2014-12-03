#! /bin/bash

################## Paras need to be adjusted for diff samples ##############
#INALL=(D05685 D170484 F05163 D08064 D171748 F05460 D10020 D171987 F05506 D100563 D174381 F06341 D11016 D180850 F07191 D114041 D36421 F08196 D115022 D62925 F0852 D122641 D83732 F08560 D122836,D99945,D62925 D86968 F08744 D122836 D88892 F09507,F09508 D129266 D9646 F09591 D134130 D99945 F10149 D13415 F00603 F11151 D146989 F02209 F13349 D163417 F03450 F98212 D164942 F03513 F98407 D16618,D050227,D62018 F03566 F98468 D169918 F03665 F99129,F99359 F07230 NG-7460_HT NG-7460_RB)

#INALL=(241013P_D186893)

#INALL=(D197292,D197408,D197290,D168414 D72902 D43046 F01426 D40211 F0193 L1405621_004,D132669,D129374 L1405621_004,F1159,F09660 D185139,F13546 F09207 241013F_D186908,241013P_D186893 D190098 D175429,D185425,D185472,D185428 F99248 F13447 D14_18979,D14_18978,GNB8513 D134555,D134427 D197723,D197722 D118569 D198023,D198022)

INALL=(F12_437_Fibro_p5 F12_94_Fibro_pX_3 SB-Ad3_Fibro_p7 F12_296_iPSC_Cl1_p26 F12_445_Fibro_p4 HLHS_Fibro_p11 SB-Neo1_Fibro_p15 F12_437_Fibro_p5,F12_94_Fibro_pX_3,F12_296_iPSC_Cl1_p26,F12_445_Fibro_p4,SB-Ad3_Fibro_p7,SB-Neo1_Fibro_p15)

#INALL=(A2463_19,A2463_21,A2463_17,A2463_16)
SAMPLEDIR="/home/a5907529/WORKING_DATA/Linda_20141024/scripts/3_MergeV" # the folder contain the VCF file which contain all of the samples
INPUT_VCF="${SAMPLEDIR}/merged_Linda_20141029.vcf"
AnnovarDIR="/sharedlustre/IGM/annovar_2014jul14"

ANOVAR_DONE="N" #use to indicate if Anovar output has been produced, Y means yes, then After_Anovar script will be used
PERL_SCRIPT_DIR=$(pwd) # change if the annotation perl script is stored somewhere else
PERL_SCRIPT="${PERL_SCRIPT_DIR}/VCF_2_annotated_excel_20141010.pl"
PERL_SCRIPT_AF_ANOVAR="${PERL_SCRIPT_DIR}/VCF_2_annotated_excel_20141010_After_Annovared.pl" #note this name has to contain 'After_Annovared', otherwise GATK will start to select samples

GENO_CALL_FLAG="Yes"

BATCH_MAF="${SAMPLEDIR}/merged_Linda_20141029_batch_MAF.txt"
CNV="/home/a5907529/WORKING_DATA/Sophie_A2463/scripts/AnnotateV/CNV_all_A2463_merged_ID_changed.txt"
InterestedGENES="${PERL_SCRIPT_DIR}/Majlinda_of_genes_20141029.txt"

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
	if [ ${ANOVAR_DONE} == "Y" ]
	then
		SCRIPT="${PERL_SCRIPT_AF_ANOVAR}"
	else
		SCRIPT="${PERL_SCRIPT}"
	fi
	#echo "${SAMPLEDIR}/VCF_filters.sh $SCRIPT ${SAMPLE_SELECTED_VCF} $OUTPUT $OUTPUT_ALL $INPUT_VCF $SAMPLE $REF_FILE $CNV $InterestedGENES"
	arr=("$SCRIPT" "${SAMPLE_SELECTED_VCF}" "$OUTPUT" "$OUTPUT_ALL" "$INPUT_VCF" "$SAMPLE" "$CNV" "$InterestedGENES" "$AnnovarDIR" "${PERL_SCRIPT_DIR}" "${BATCH_MAF}" "${OUTPUT_VCF}")
	qsub ${PERL_SCRIPT_DIR}/VCF_filters_20140811.sh "${arr[@]}"
# 	sh ${PERL_SCRIPT_DIR}/VCF_filters_20140811.sh "${arr[@]}"
done

