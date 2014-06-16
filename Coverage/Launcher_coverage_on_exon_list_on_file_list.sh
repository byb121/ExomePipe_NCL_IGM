#! /bin/bash

################## Paras need to be adjusted for diff samples #########################

#INALL=(D050227 D05685 D08064) #Sample names, corespond to each bam file in the BamFile.list, can be any string as long as not duplicated
INALL=(PFC-348-GP_S19 PFC-353-AA_S24)
SAMPLE_PATH=`pwd` #Path to ouput bedtools and summary files
INPUT_LIST="${SAMPLE_PATH}/BamFile.Padraig" #Bamfile list, full path to the bam files, one file each line
SCRIPT_PATH=$SAMPLE_PATH #Path to the betools and coverage summary scripts
TARGETS="${SAMPLE_PATH}/POLRMT.bed" #Regions to check coverage on in bed format.It can only have 4 column with the unique IDs on the 4th column
COV_tail="_POLRMT.txt"  #suffix of the bedtools output file
OUTPUT_prefix="POLRMT_353_348" # prefix of summary files
################## Paras need to be adjusted for diff samples #########################



#############################################################################################################
readarray INBAM < ${INPUT_LIST}
COUNT=0
for SAMPLE_ID in "${INALL[@]}"
do	
	BAM=${INBAM[$COUNT]}
	if [ $COUNT -eq 0 ]
        then
                SAMPLE_ID_LIST="${SAMPLE_ID}"
        else
                SAMPLE_ID_LIST="${SAMPLE_ID_LIST},${SAMPLE_ID}"
        fi
	#echo $BAM ${SAMPLE_ID}
	JOB_ID2="Bed_Cov_exon_list"
	OUTPUT_COV="${SAMPLE_PATH}/${SAMPLE_ID}${COV_tail}"
	arr=("$BAM" "$TARGETS" "$OUTPUT_COV")
	qsub -N $JOB_ID2 ${SCRIPT_PATH}/bedtools.sh "${arr[@]}"
	COUNT=$(($COUNT+1))
done

JOB_ID3="Cal_Cov_Exon"
arr2=("$SAMPLE_PATH" "$SAMPLE_ID_LIST" "$OUTPUT_prefix" "$COV_tail" "$SCRIPT_PATH")
qsub -hold_jid $JOB_ID2 -N $JOB_ID3 ${SCRIPT_PATH}/coverage_summary_on_exon_list_on_file_list.sh "${arr2[@]}"

