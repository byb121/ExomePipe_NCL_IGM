ExomePipe_NCL_IGM
=================

Intro
-----------------

This is a package of scripts for exome-seq data analysis. It utilizes a HPC managed by Sun Grid Engine and focuses on germ line variants detection exclusively for now. The scripts used in each stage of the analysis are grouped in separated folders.

There are two types of scripts right now:
* Scripts doing an specific job, and
* Scripts controlling the workflow and submit ordered scripts to the queue engine of a HPC.

The 1st kind of scripts requires no change or minor change from project to project, but the 2nd type of scripts is required to change in order to fit a new project setup.

There are two main pipelines grouped into two folders:

* 1. GATK_Variant_Call_Pipe
Mainly to achieve the best practice suggested in GATK forum. The pipeline does mapping with BWA -> samtools -> Indel Realignment -> Base Quality Recalibration -> HaplotypeCaller -> Variants Quality Recalibration -> Annovar for variant filteration.

update (Nov 2014):
The pipeline has changed according to the new best practice (GATK V3.1 and above). GVCF files are produced after alignment, scalar scripts for HaplotypeCaller is no longer needed. The pipeline is now using FastUniq to remove duplicated reads instead of PICARD. PICARD is still used to convert SAM files into sorted BAMs. The pipeline does the following jobs:
FastUniq -> BWA -> PICARD -> Indel Realignment -> Base Quality Recalibration -> HaplotypeCaller(GVCFs) -> GenotypeGVCFs from GVCFs -> Variant Quality Recalibration -> Annovar for variant annotation and filtering


* 2. Joint_Variant_Call_Pipe
We constantly find GATK HaplotypeCaller misses true variants, while meantime Freebayes become more popular. In August 2014, I added pipeline to do a simple joint variant calls using GATK HaplotypeCaller and Freebayes. GVCF format is used as suggest in the new GATK Best Practice for GATK 3.2. 
The pipeline does:
Duplication removal with FastUniq -> mapping with BWA mem -> PICARD to sort sam and convert to bam -> Indel Realignment -> Base Quality Recalibration -> HaplotypeCaller to produce GVCFs, then call variants with Freebayes from PICARD output bam files and GATK GenotypeGVCFs from GVCFs -> Quality Recalibration of GATK output -> joint variants calls -> Annovar for variant annotation and filtering.


Sample Folder Structure
-----------------

In order to use the script with minimize changes on a new project, the fastq files should be organized in a specific structure.

Ideally, each biological sample should have its own folder under a project folder; fastq files for each sample should be under a same folder. Scripts can be under the same project folder. Fastq files should have same pattern of file names. This need to be reflected in the following scripts in both pipelines so that fastq files can be located:

* detectSampleLanes.pl
* map_recali_perLane_recali_perSample_covOnTargets_GVCF_31Jul2014.sh
* redupFastuniq_31Jul2014.sh


Other things
-----------------
* FastUniq, autoadapt may need to be installed
* Annovar installation and annotations may be very different, so the annotation scripts may need to change a lot

