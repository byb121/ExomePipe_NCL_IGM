#! /bin/bash

################## Paras need to be adjusted for diff samples ##############
INALL=(F12_296_iPSC_Cl1_p26 F12_437_Fibro_p5 F12_445_Fibro_p4 F12_94_Fibro_pX_3 HLHS-1 HLHS-2 HLHS_Fibro_p11 J-HLHS4 J-HLHS5 SB-Ad3_Fibro_p7 SB-Neo1_Fibro_p15)

SAMPLEDIR="/home/a5907529/WORKING_DATA/Linda_20141024/scripts/3_MergeV" # the folder contain the VCF file which contain all of the samples
INPUT_VCF="${SAMPLEDIR}/merged_Linda_Rachel_2015feb26.vcf"
AnnovarDIR="/sharedlustre/IGM/annovar_2014jul14"

ANOVAR="Y" #use to indicate if Anovar output has beed produced, Y means yes, then After_Anovar script will be used
PERL_SCRIPT_DIR=$(pwd) # change if the annotation perl script is stored somewhere else
PERL_SCRIPT="${PERL_SCRIPT_DIR}/VCF_2_annotated_excel_20150317.pl"
PERL_SCRIPT_AF_ANOVAR="${PERL_SCRIPT_DIR}/VCF_2_annotated_excel_20150317_After_Annovared.pl" #note this name has to contain 'After_Annovared', otherwise GATK will start to select samples

GENO_CALL_FLAG="Yes"

BATCH_MAF="/home/a5907529/WORKING_DATA/Linda_20141024/scripts/3_MergeV/merged_Linda_Rachel_2015feb26_batch_MAF.txt"
CNV="/home/a5907529/WORKING_DATA/Linda_20141024/ExomeDepth/Linda_all_CNVs.txt"
InterestedGENES="/home/a5907529/WORKING_DATA/Linda_20141024/scripts/4_AnnotateV/Majlinda_of_genes_20141029.txt"

VCF_TAIL="_selected.vcf"
FILTERED_TAIL="_filtered.xlsx" # output excel file name tail of filtered variants, has to be ended with xlsx
EVERYTHING_TAIL="_Everything.xlsx" # output excel file name tail of all variants, has to be ended with xlsx
CUSTOM_VCF_TAIL="_selectet.custom.vcf"

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
	arr=("$SCRIPT" "${SAMPLE_SELECTED_VCF}" "$OUTPUT" "$OUTPUT_ALL" "$INPUT_VCF" "$SAMPLE" "$CNV" "$InterestedGENES" "$AnnovarDIR" "${PERL_SCRIPT_DIR}" "${BATCH_MAF}" "${OUTPUT_VCF}" )
	qsub ${PERL_SCRIPT_DIR}/VCF_filters_20150317.sh "${arr[@]}"

 	#sh ${PERL_SCRIPT_DIR}/VCF_filters_20140811.sh "${arr[@]}"
done

