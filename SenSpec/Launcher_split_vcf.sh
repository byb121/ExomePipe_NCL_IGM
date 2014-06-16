#! /bin/bash

## Script to launch Exome Pipeline (Picard - GATK bam recalibration - GATK call variants & varscan & dindel & bedtools for coverage) for multiple samples
## INALL is an array of directory names, with one directory per sample
## SAMPLE_PATH is the pathway to the sample directories ${SAMPLE_PATH}/${INALL[$SAMPLE_NUMBER]}


################## Paras need to be adjusted for diff samples #########################
#INALL=(PFC-304-OW2_S1 PFC-305-EW2_S2 PFC-309-LS_S3 PFC-334-HC_S4 PFC-335-EY_S5 PFC-336-TF_S6 PFC-337-MV_S7 PFC-338-BJ_S8 PFC-339-TM_S9 PFC-340-MH_S10 PFC-341-JS_S13 PFC-342-JF_S14 PFC-344-DOH_S15 PFC-345-WC_S16 PFC-346-CM_S17 PFC-347-DS_S18 PFC-348-GP_S19 PFC-349-ES_S20 PFC-350-LE_S21 PFC-351-FSA_S22 PFC-352-EAL_S23 PFC-353-AA_S24 PFC-356-AF_S42 PFC-357-JW_S11 PFC-358-BC_S12 PFC-359-TC_S25 PFC-360-LW_S26 PFC-361-KD_S27 PFC-362-NR_S28 PFC-363-LV_S29 PFC-364-YB_S30 PFC-365-NH_S31 PFC-366-SS_S32 PFC-367-SC_S33 PFC-368-PT_S34 PFC-369-CC_S35 PFC-370-HLP_S48 PFC-371-BP_S36 PFC-372-AE_S37 PFC-373-DD_S38 PFC-374-WG_S39 PFC-375-PAN_S40 PFC-376-HA_S41 PFC-378-CAB_S43 PFC-379-JAW_S44 PFC-381-BA_S45 PFC-382-DW_S46 PFC-383-ALC_S47)

INALL=(PFC-304-OW2_S1)

INPUT_VCF="/home/a5907529/WORKING_DATA/QC_GOOD/scripts_FreeBayes/merged.vcf"
OUTPUT_PATH="/home/a5907529/WORKING_DATA/QC_GOOD/scripts_FreeBayes/SenSpec"
SCRIPT_PATH=${OUTPUT_PATH}
REF="/sharedlustre/IGM/bundle2.8/hg19/ucsc.hg19.fasta"
###############################

for SAMPLE_ID in "${INALL[@]}"
do
	#echo "INPUT_VCF:${INPUT_VCF}, OUTPUT_PATH:${OUTPUT_PATH}, SAMPLE_ID:${SAMPLE_ID}"
	qsub ${SCRIPT_PATH}/GATK_select_sample.sh ${INPUT_VCF} ${OUTPUT_PATH} ${SAMPLE_ID} $REF
done

