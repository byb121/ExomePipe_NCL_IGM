#! /bin/bash

INALL=(F12_296_iPSC_Cl1_p26 F12_445_Fibro_p4 F12_437_Fibro_p5 F12_94_Fibro_pX_3 SB-Ad3_Fibro_p7 SB-Neo1_Fibro_p15)

OUT_DIR="/users/a5907529/lustre/Prof_Linda_First6"
SAMPLE_PATH="/users/a5907529/data/ForMauro_2014-06-12/F14FTSAPHT0227_HUMjbhX/result"
SCRIPTDIR="$OUT_DIR/scripts/ExomePipe_2014JUNE/QC"

OUTPUT_DIR_FastQC="${OUT_DIR}/FastQC_raw"
mkdir ${OUTPUT_DIR_FastQC}

for SAMPLE_ID in "${INALL[@]}"
do      
	SAMPLE_DIR="${SAMPLE_PATH}/${SAMPLE_ID}/clean_data"
	OUTPUT_DIR_Unzipped="${OUT_DIR}/${SAMPLE_ID}"
	mkdir ${OUTPUT_DIR_Unzipped}
	ZIPPED_FASTQS="${SAMPLE_DIR}/${SAMPLE_ID}*fq.gz"
	for f in ${ZIPPED_FASTQS}; do
		echo "$f"
    		qsub $SCRIPTDIR/FastQC.sh $f ${OUTPUT_DIR_FastQC}
		qsub $SCRIPTDIR/unzip.sh $f ${OUTPUT_DIR_Unzipped}
	done
done
