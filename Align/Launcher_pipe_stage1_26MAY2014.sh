#! /bin/bash

#**********************************  important!!  *******************************************
#if use the script after trimming: fastq file name patterns need to be changed as following: 
# comment out - "map_recali_perLane_recali_perSample_covOnTargets_26MAY2014.sh" line 112, 113
# and use line 108, 109 
# the comment out - "detectSampleLanes.pl" line 13 use line 16
#********************************************************************************************

## Script to launch Exome Pipeline (Picard - GATK bam recalibration - GATK call variants & varscan & dindel & bedtools for coverage) for multiple samples
## INALL is an array of directory names, with one directory per sample
## SAMPLE_PATH is the pathway to the sample directories ${SAMPLE_PATH}/${INALL[$SAMPLE_NUMBER]}


#/************************* Paras need to be adjusted for diff samples
INALL=(PFC-304-OW2_S1 PFC-305-EW2_S2 PFC-309-LS_S3 PFC-334-HC_S4 PFC-335-EY_S5 PFC-336-TF_S6 PFC-337-MV_S7 PFC-338-BJ_S8 PFC-339-TM_S9 PFC-340-MH_S10 PFC-341-JS_S13 PFC-342-JF_S14 PFC-344-DOH_S15 PFC-345-WC_S16 PFC-346-CM_S17 PFC-347-DS_S18 PFC-348-GP_S19 PFC-349-ES_S20 PFC-350-LE_S21 PFC-351-FSA_S22 PFC-352-EAL_S23 PFC-353-AA_S24 PFC-356-AF_S42 PFC-357-JW_S11 PFC-358-BC_S12 PFC-359-TC_S25 PFC-360-LW_S26 PFC-361-KD_S27 PFC-362-NR_S28 PFC-363-LV_S29 PFC-364-YB_S30 PFC-365-NH_S31 PFC-366-SS_S32 PFC-367-SC_S33 PFC-368-PT_S34 PFC-369-CC_S35 PFC-370-HLP_S48 PFC-371-BP_S36 PFC-372-AE_S37 PFC-373-DD_S38 PFC-374-WG_S39 PFC-375-PAN_S40 PFC-376-HA_S41 PFC-378-CAB_S43 PFC-379-JAW_S44 PFC-381-BA_S45 PFC-382-DW_S46 PFC-383-ALC_S47) 

#INALL=(test1 test2)

SAMPLE_PATH="/home/a5907529/WORKING_DATA/QC_GOOD"
SCRIPT_PATH="${SAMPLE_PATH}/scripts_partial_QC"
TARGETS="/home/a5907529/WORKING_DATA/Project_Fulcrum/nexterarapidcapture_exome_targetedregions03072013.bed"
Library_ID="Kostas.Patrick" #Better to use {PI name}.{Batch Number} to identify the library

REFDIR="/sharedlustre/IGM/bundle2.8/hg19"
#SCRATCH_DIR="/users/a5907529/scratch"
SCRATCH_DIR="/home/a5907529/WORKING_DATA/pipe.temp.partial_QC"
SEQ_PLATFORM="ILLUMINA" # Valid values are: ILLUMINA, SOLID, LS454, HELICOS and PACBIO
COV_OUT_DIR_NAME="coverage_partial_QC" #coverage file output directory name. Just a name, it will be placed under the sample directory
GATK_OUT_DIR_NAME="gatk_partial_QC" #gatk tuned alignment file output directory name
WRKGDIR_NAME="_partial_QC" # tail of temp working directory under scratch dir, in case same sample running at the same using the same temp directory
JAVA_TMP_DIR_NAME="_partial_QC" #tail of temp java working directory under scratch dir
#**************************/

rm -r ${SCRATCH_DIR}
mkdir ${SCRATCH_DIR}

###############################

## Submitting jobs ##
for SAMPLE_ID in "${INALL[@]}"
do
       	JOB_ID1="pipe_stage1_${SAMPLE_ID}"
	arr=("${SAMPLE_ID}" "${SAMPLE_PATH}" "${SCRIPT_PATH}" "$REFDIR" "${SCRATCH_DIR}" "$TARGETS" "${SEQ_PLATFORM}" "${Library_ID}" "${COV_OUT_DIR_NAME}" "${GATK_OUT_DIR_NAME}" "${WRKGDIR_NAME}" "${JAVA_TMP_DIR_NAME}")
        qsub ${SCRIPT_PATH}/map_recali_perLane_recali_perSample_covOnTargets_26MAY2014.sh "${arr[@]}"
done

