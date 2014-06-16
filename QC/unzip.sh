#! /bin/bash
#$ -cwd
#$ -j y
#$ -V

FILE=$1
OUT_DIR=$2

filename="${FILE##*/}"                      # Strip longest match of */ from start
#dir="${fullpath:0:${#fullpath} - ${#filename}}" # Substring from 0 thru pos of filename
base="${filename%.[^.]*}"                  # Strip shortest match of . plus at least one non-dot char from end
OUT_FILE="${OUT_DIR}/$base"
gunzip -c $FILE > ${OUT_FILE}
echo "done."
