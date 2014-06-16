#! /bin/bash
#$ -cwd
#$ -j y
#$ -V

FILE=$1
OUT_DIR=$2
SAMPLE_NAME=$3

filename="${FILE##*/}"                      # Strip longest match of */ from start
#dir="${fullpath:0:${#fullpath} - ${#filename}}" # Substring from 0 thru pos of filename
base="${filename%.[^.]*}"                  # Strip shortest match of . plus at least one non-dot char from end

#/***************** BGI data specific
base2=`echo $base | grep -o -P "\_L\d+\_"`
base3=`echo $base | grep -o -P "\d\..*\.fq$"`
#*****************/

OUT_FILE="${OUT_DIR}/${SAMPLE_NAME}$base2$base3"
echo "$FILE ${OUT_FILE}"
#gunzip -c $FILE > ${OUT_FILE}
echo "done."
