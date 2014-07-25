#! /bin/bash
#$ -cwd -V
#$ -j y
#$ -pe smp 4
#$ -l mem_free=10G,h_vmem=15G
#$ -M yaobo.xu@ncl.ac.uk
#$ -m e

set -e
res1=$(date +%s.%N) # use to calculate whole run time of the job
#cwd=`pwd`; export cwd

#/************************************
# note: fastq file patterns are used in the perl script detectSampleLanes.pl and READ_FILE1 & READ_FILE2
# they need to be adjusted accordingly
#************************************/

echo $'\n'"["`date`"]: Job stared."

SAMPLE_ID=$1
SAMPLE_PATH=$2
SCRIPTS_DIR=$3
REF_DIR=$4
SCRATCH_DIR=$5
TARGETS=$6
SEQ_PLATFORM=$7 # Valid values are: ILLUMINA, SOLID, LS454, HELICOS and PACBIO
Library_ID=$8 #Better to use {PI name}.{Batch Number} to identify the library
COV_DIR_NAME=$9
GATK_OUT_DIR_NAME=${10}
DUP_FREE_BAM_DIR_NAME=${11}
WRKGDIR_NAME=${12}
JAVA_TMP_DIR_NAME=${13}

#/************ less frequently changed options
REF_FILE="${REF_DIR}/human_g1k_v37_decoy.fasta"
DBSNP_VCF="${REF_DIR}/dbsnp_138.b37.vcf"
# GATK base quality recalibration covariates
COVARIATES="-cov ReadGroupCovariate -cov QualityScoreCovariate -cov ContextCovariate -cov CycleCovariate"
#************/


#/************************ fmscluster add modules
# need to add line here to load fastuniq module
#GATK 3.1.1 /home/a5907529/WORKING_DATA/QC_GOOD/GATK3
#module load apps/samtools/0.1.19
#module load apps/bwa/0.7.6.a
#module load apps/picard/1.107
#module load apps/bedtools/2.19.0
#module load apps/perl/5.18.2
#module add apps/java
#SCRATCH_DIR=$TMPDIR
#export _JAVA_OPTIONS="-XX:-UseLargePages -Xms4096m"
#PICARDDIR="$PICARD_PATH"
#GATKDIR="/home/a5907529/WORKING_DATA/QC_GOOD/GATK3"
#INDIR="${SAMPLE_PATH}/${SAMPLE_ID}"
#WRKGDIR="${SCRATCH_DIR}/${SAMPLE_ID}${WRKGDIR_NAME}"
#JAVA_TMP_DIR="${SCRATCH_DIR}/${SAMPLE_ID}.Java.tmp.dir.${JAVA_TMP_DIR_NAME}"
#echo $'\n'"mkdir $JAVA_TMP_DIR"
#mkdir $JAVA_TMP_DIR
#Picard_sort="java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx8g -jar $PICARDDIR/SortSam.jar VALIDATION_STRINGENCY=LENIENT"
#Picard_nodups="java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx8g -jar $PICARDDIR/MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT"
#Picard_CleanSam="java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx8g -jar $PICARDDIR/CleanSam.jar VALIDATION_STRINGENCY=LENIENT"
#*************************/

#/************************ lampredi2 add modules
module load apps/fastuniq/1.1/gcc-4.4.6
module load apps/samtools/0.1.18/gcc-4.4.6
module load apps/bwa/0.7.4/gcc-4.4.6
module load apps/picard/1.85/noarch
module load apps/bedtools/2.17.0/gcc-4.4.6
module load apps/perl/5.16.1/gcc-4.4.6
module load apps/gatk/3.1.1/noarch
INDIR="${SAMPLE_PATH}/${SAMPLE_ID}"
WRKGDIR="${SCRATCH_DIR}/${SAMPLE_ID}${WRKGDIR_NAME}"
JAVA_TMP_DIR="${SCRATCH_DIR}/${SAMPLE_ID}.Java.tmp.dir.${JAVA_TMP_DIR_NAME}"
echo $'\n'"mkdir $JAVA_TMP_DIR"
mkdir $JAVA_TMP_DIR
Picard_sort="java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx8g -jar $PICARDDIR/java/SortSam.jar VALIDATION_STRINGENCY=LENIENT"
Picard_nodups="java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx8g -jar $PICARDDIR/java/MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT"
Picard_CleanSam="java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx8g -jar $PICARDDIR/java/CleanSam.jar VALIDATION_STRINGENCY=LENIENT"
#*************************/

#Output folders' names
COV_DIR="${INDIR}/${COV_DIR_NAME}"
GATK_OUT_DIR="${INDIR}/${GATK_OUT_DIR_NAME}"
DUP_FREE_BAM_DIR="${INDIR}/${DUP_FREE_BAM_DIR_NAME}"

# prepare folders
echo $'\n'mkdir $WRKGDIR
if [ ! -d $WRKGDIR ]; then
	mkdir $WRKGDIR
else 
	echo "$WRKGDIR exists"	
fi
echo $'\n'mkdir $COV_DIR
if [ ! -d $COV_DIR ]; then
	mkdir $COV_DIR
else
        echo "$COV_DIR exists"   
fi
echo $'\n'mkdir $GATK_OUT_DIR
if [ ! -d $GATK_OUT_DIR ]; then
	mkdir $GATK_OUT_DIR
else
        echo "$GATK_OUT_DIR exists"   
fi
echo $'\n'mkdir $DUP_FREE_BAM_DIR
if [ ! -d $DUP_FREE_BAM_DIR ]; then
	mkdir ${DUP_FREE_BAM_DIR}
else
        echo "${DUP_FREE_BAM_DIR} exists"   
fi

PICARD_TEMP="$WRKGDIR/Picard_Temp"
echo $'\n'"mkdir $PICARD_TEMP"
if [ ! -d $WRKGDIR ]; then
	mkdir $PICARD_TEMP
else
        echo "$PICARD_TEMP exists"   
fi
PICARD_LOG="$WRKGDIR/${SAMPLE_ID}_picard.log"

################## BWA Bit ####################
LANES_STRING=`perl ${SCRIPTS_DIR}/detectSampleLanes.pl $SAMPLE_PATH $SAMPLE_ID`
LANES=($LANES_STRING)

if [ ${#LANES[@]} -lt 1 ]; then
	exit "Error: No lane info detected, check file name patterns and paths"
fi

BAM_FILE_LIST=""
echo $'\n'"started to loop through lanes"
echo "----------------------------- hua li li de fen jie xian -----------------------------"
echo "----------------------------- hua li li de fen jie xian -----------------------------"
for LANE in "${LANES[@]}"
do
	echo "---------------- hua li de fen jie xian ----------------"
	echo "process lane $LANE for sample ${SAMPLE_ID}"
	SAM_FILE1="$WRKGDIR/${SAMPLE_ID}_${LANE}.sam"
	BAM_FILE1="$WRKGDIR/${SAMPLE_ID}_${LANE}.bam" #final lane bam after GATK base recalibration
	BAM_FILE_LIST="$BAM_FILE_LIST-I ${BAM_FILE1} "
	
	#/********************** if QC ed with the Trimming.sh
	#READ_FILE1="$INDIR/${SAMPLE_ID}_L00${LANE}*val_1.fq"
        #READ_FILE2="$INDIR/${SAMPLE_ID}_L00${LANE}*val_2.fq"
	#**********************/	

	READ_FILE1="$INDIR/${SAMPLE_ID}_L00${LANE}_R1_001.fastq"
	READ_FILE2="$INDIR/${SAMPLE_ID}_L00${LANE}_R2_001.fastq"
	
	# can handle zipped file now, will do unzip automatically
	if [ -f "$READ_FILE1.gz" ]; then
		echo "find gzipped fastq file $READ_FILE1.gz, unzip now.."
		gunzip -c "$READ_FILE1.gz" > "$READ_FILE1"
	else
		if [ -f $READ_FILE1 ]; then
			echo "find file $READ_FILE1"
		else
			exit "fasq file name error, can not find file $READ_FILE1"
		fi
	fi

	if [ -f "$READ_FILE2.gz" ]; then
		echo "find gzipped fastq file $READ_FILE2.gz, unzip now.."
                gunzip -c "$READ_FILE2.gz" > "$READ_FILE2"
	else
		if [ -f $READ_FILE2 ]; then
                        echo "find file $READ_FILE2"
                else
                        exit "fasq file name error, can not find file $READ_FILE2"
                fi
        fi
	
	READ_FILE1_nodup="$INDIR/${SAMPLE_ID}_L00${LANE}_R1_001.nodup.fastq"
        READ_FILE2_nodup="$INDIR/${SAMPLE_ID}_L00${LANE}_R2_001.nodup.fastq"
	FASTQ_LIST_FILE="$INDIR/${SAMPLE_ID}_L00${LANE}.temp.qlist"
	# remvoe dups with FastUniq
	echo -e "$READ_FILE1\n$READ_FILE2" > $FASTQ_LIST_FILE
	echo fastuniq -i "${FASTQ_LIST_FILE}" -t q -o "${READ_FILE1_nodup}" -p "${READ_FILE2_nodup}"
	fastuniq -i "${FASTQ_LIST_FILE}" -t q -o "${READ_FILE1_nodup}" -p "${READ_FILE2_nodup}"
	
	RAW_BAM="$SAM_FILE1.bam"
	DUP_FREE_BAM_FILE_LIST="$DUP_FREE_BAM_FILE_LIST-I ${RAW_BAM} "
	
	#bwa mapping
	echo $'\n'"["`date`"]: bwa starts to align the lane ${LANE}."
	echo bwa mem -R "@RG\tID:FlowCell.${LANE}.${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:${SEQ_PLATFORM}\tLB:${Library_ID}.${SAMPLE_ID}" -t 4 -M $REF_FILE $READ_FILE1_nodup $READ_FILE2_nodup "> $SAM_FILE1"
	bwa mem -R "@RG\tID:FlowCell.${LANE}.${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:${SEQ_PLATFORM}\tLB:${Library_ID}.${SAMPLE_ID}" -t 4 -M $REF_FILE $READ_FILE1_nodup $READ_FILE2_nodup > $SAM_FILE1 # "-M" is for Picard compatibility

	rm ${READ_FILE1}
	rm ${READ_FILE2}
	rm ${READ_FILE1_nodup}
	rm ${READ_FILE2_nodup}
	
	echo $'\n'"["`date`"]: PICARD to sort the sam file ${SAM_FILE1}"
	echo "$Picard_sort INPUT=$SAM_FILE1 OUTPUT=$RAW_BAM SORT_ORDER=coordinate TMP_DIR=$PICARD_TEMP"
	$Picard_sort INPUT=$SAM_FILE1 OUTPUT=$RAW_BAM SORT_ORDER=coordinate
	echo $'\n'"samtools index $RAW_BAM"
        samtools index $RAW_BAM
	
	#GATK to refine the alignment
	INTERVALS_LANE="$WRKGDIR/${SAMPLE_ID}_${LANE}_nodups.sorted.bam.intervals"
	REALIGNED_BAM_LANE="$WRKGDIR/${SAMPLE_ID}_${LANE}_nodups.sorted.realigned.bam"
	RECAL_TABLE_GRP_LANE="$WRKGDIR/${SAMPLE_ID}_${LANE}_nodups.sorted.realigned.bam.grp"

	echo $'\n'"["`date`"]:GATK: Creating realignment intervals for $RAW_BAM"
	echo java -Xmx10g -jar $GATKDIR/GenomeAnalysisTK.jar -rf BadCigar -T RealignerTargetCreator -nt 4 -I $RAW_BAM -R $REF_FILE -o $INTERVALS_LANE
	java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx12g -jar $GATKDIR/GenomeAnalysisTK.jar -rf BadCigar -T RealignerTargetCreator -nt 4 -I $RAW_BAM -R $REF_FILE -o $INTERVALS_LANE
	echo $'\n'"["`date`"]:GATK: Realigning reads..."
	echo java -Xmx10g -jar $GATKDIR/GenomeAnalysisTK.jar -rf BadCigar -T IndelRealigner -I $RAW_BAM -R $REF_FILE -targetIntervals $INTERVALS_LANE -o $REALIGNED_BAM_LANE
	java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx12g -jar $GATKDIR/GenomeAnalysisTK.jar -rf BadCigar -T IndelRealigner -I $RAW_BAM -R $REF_FILE -targetIntervals $INTERVALS_LANE -o $REALIGNED_BAM_LANE
	echo "rm $INTERVALS_LANE"
	rm $INTERVALS_LANE
	
	echo $'\n'"["`date`"]:GATK: Calculating recalibration tables..."
	echo java -Xmx10g -jar $GATKDIR/GenomeAnalysisTK.jar -rf BadCigar -T BaseRecalibrator -nct 4 -I $REALIGNED_BAM_LANE -R $REF_FILE $COVARIATES -knownSites $DBSNP_VCF -o $RECAL_TABLE_GRP_LANE
	java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx12g -jar $GATKDIR/GenomeAnalysisTK.jar -rf BadCigar -T BaseRecalibrator -nct 4 -I $REALIGNED_BAM_LANE -R $REF_FILE $COVARIATES -knownSites $DBSNP_VCF -o $RECAL_TABLE_GRP_LANE
	echo $'\n'"["`date`"]:GATK: Creating Recalibrated alignment file..."
	echo java -Xmx6g -Djava.io.tmpdir=${JAVA_TMP_DIR} -jar $GATKDIR/GenomeAnalysisTK.jar -rf BadCigar -T PrintReads -nct 4 -R $REF_FILE -I $REALIGNED_BAM_LANE -BQSR $RECAL_TABLE_GRP_LANE -o $BAM_FILE1
	java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx12g -jar $GATKDIR/GenomeAnalysisTK.jar -rf BadCigar -T PrintReads -nct 4 -R $REF_FILE -I $REALIGNED_BAM_LANE -BQSR $RECAL_TABLE_GRP_LANE -o $BAM_FILE1
	echo rm $REALIGNED_BAM_LANE
	rm $REALIGNED_BAM_LANE
	echo rm $RECAL_TABLE_GRP_LANE
	rm $RECAL_TABLE_GRP_LANE

	echo $'\n'"["`date`"]:GATK: Indexing bam files..."
	echo "samtools index $BAM_FILE1"
	echo "---------------- hua li de fen jie xian ----------------"
	samtools index $BAM_FILE1
done
echo $'\n'"------------------------------"
echo "all lanes have been processed!"
echo $'\n'"------------------------------"

echo $'\n'"Merge alignment and move files"
# Merge and clean
FINAL_BAM="${GATK_OUT_DIR}/${SAMPLE_ID}_nodups.realigned.recalibrated.bam"
DUP_FREE_BAM="${DUP_FREE_BAM_DIR}/${SAMPLE_ID}_nodups.bam"

if [ ${#LANES[@]} -eq 1 ]
then
	echo "mv $BAM_FILE1 $FINAL_BAM"
	mv $BAM_FILE1 $FINAL_BAM
	mv $RAW_BAM $DUP_FREE_BAM
else # merge with GATK printReads
	MERGED_BAM="$WRKGDIR/${SAMPLE_ID}.bam"
	echo $'\n'java -Xmx10g -Djava.io.tmpdir=${JAVA_TMP_DIR} -jar $GATKDIR/GenomeAnalysisTK.jar -T PrintReads -nct 4 -R $REF_FILE $BAM_FILE_LIST -o $MERGED_BAM
	java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx12g -jar $GATKDIR/GenomeAnalysisTK.jar -T PrintReads -nct 4 -R $REF_FILE $BAM_FILE_LIST -o $MERGED_BAM
	echo "samtools index $MERGED_BAM"	
	samtools index $MERGED_BAM
	
	echo $'\n'java -Xmx10g -Djava.io.tmpdir=${JAVA_TMP_DIR} -jar $GATKDIR/GenomeAnalysisTK.jar -T PrintReads -nct 4 -R $REF_FILE $DUP_FREE_BAM_FILE_LIST -o $DUP_FREE_BAM
        java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx12g -jar $GATKDIR/GenomeAnalysisTK.jar -T PrintReads -nct 4 -R $REF_FILE $DUP_FREE_BAM_FILE_LIST -o $DUP_FREE_BAM
        echo samtools index $DUP_FREE_BAM
        samtools index $DUP_FREE_BAM
	

	#clean
	for LANE in "${LANES[@]}"
	do
		BAM_FILE1="$WRKGDIR/${SAMPLE_ID}_${LANE}.bam"
        	RAW_BAM="$WRKGDIR/${SAMPLE_ID}_${LANE}.sam.bam"
		echo rm $BAM_FILE1
		rm $BAM_FILE1
		echo rm $RAW_BAM
                rm $RAW_BAM
	done
		
	# GATK recalibration at sample level #
	INTERVALS="$WRKGDIR/${SAMPLE_ID}_nodups.bam.intervals"
	echo $'\n'"["`date`"]: GATK: Creating realignment intervals..."
	echo java -Xmx10g -jar $GATKDIR/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 4 -I $MERGED_BAM -R $REF_FILE -o $INTERVALS
	java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx12g -jar $GATKDIR/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 4 -I $MERGED_BAM -R $REF_FILE -o $INTERVALS
	echo $'\n'"["`date`"]: GATK: Realigning reads..."
	echo java -Xmx10g -jar $GATKDIR/GenomeAnalysisTK.jar -T IndelRealigner -I $MERGED_BAM -R $REF_FILE -targetIntervals $INTERVALS -o $FINAL_BAM
	java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx12g -jar $GATKDIR/GenomeAnalysisTK.jar -T IndelRealigner $MERGED_BAM -R $REF_FILE -targetIntervals $INTERVALS -o $FINAL_BAM
	echo rm $MERGED_BAM
	rm $MERGED_BAM
	echo rm $INTERVALS
	rm $INTERVALS
	echo rm $SORTED_BAM_FILE_NODUPS
	rm $SORTED_BAM_FILE_NODUPS
fi
echo samtools index $FINAL_BAM
samtools index $FINAL_BAM

# cleaning bit #
echo rm -r $PICARD_TEMP
rm -r $PICARD_TEMP
echo rm -r $WRKGDIR
rm -r $WRKGDIR
echo rm -r $JAVA_TMP_DIR
rm -r $JAVA_TMP_DIR

#Coverage on targets file
BED_OUTPUT="${COV_DIR}/${SAMPLE_ID}_onTargets.txt"
echo $'\n'"coverageBed -abam $FINAL_BAM -b $TARGETS -hist -split > $BED_OUTPUT"
coverageBed -abam $FINAL_BAM -b $TARGETS -hist -split > $BED_OUTPUT

echo $'\n'"["`date`"]: The whole job is DONE!!"

# runtime calculation
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)
echo "----------------------------- hua li li de fen jie xian -----------------------------"
echo "----------------------------- hua li li de fen jie xian -----------------------------"
printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds
echo "exit status $?"
