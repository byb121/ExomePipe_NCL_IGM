ExomePipe_NCL_IGM
=================

Intro
-----------------

This is a package of scripts for exome-seq data analysis. It utilizes a HPC managed by Sun Grid Engine and focuses on germ line variants detection exclusively for now. The scripts used in each stage of the analysis are grouped in separated folders.

There are two types of scripts right now:
* Scripts doing an specific job, and
* Scripts controlling the workflow and submit ordered scripts to the queue engine of a HPC.

The 1st kind of scripts requires no change or minor change from project to project, but the 2nd type of scripts is required to change in a new project setup.

There are two main pipelines grouped in two different folders:
GATK_Variant_Call_Pipe: Mainly to achieve the best practice suggested in GATK forum. The pipeline does mapping with BWA -> samtools -> PICARD to remove duplicated reads -> Indel Realignment -> Base Quality Recalibration -> HaplotypeCaller -> Variants Quality Recalibration -> Annovar for variant filteration.

Joint_Variant_Call_Pipe: We constantly find GATK HaplotypeCaller misses true variants, while meantime Freebayes become more popular. In August 2014, I added pipeline to do a simple joint variant calls between GATK HaplotypeCaller and Freebayes, GVCF format is used as well as suggest by the GATK Best Practice for GATK 3.2.  The pipeline does duplication removal with FastUniq -> mapping with BWA mem -> PICARD to sort sam and convert to bam -> Indel Realignment -> Base Quality Recalibration -> HaplotypeCaller to produce GVCFs, then call variants with Freebayes from PICARD output bam files and GATK GenotypeGVCFs from GVCFs -> Quality Recalibration of GATK output -> joint variants calls -> Annovar for variant annotation and filtering.

Sample Folder Structure
-----------------

In order to use the script with minimize changes on a new project, the fastq files should be organized in a specific structure.

Ideally, each biological sample should have its own folder under a project folder; fastq files for each sample should be under a same folder. Scripts can be under the same project folder. Fastq files should have same pattern of file names. This need to be reflected in the following scripts in both pipelines so that fastq files can be located:

* detectSampleLanes.pl
* map_recali_perLane_recali_perSample_covOnTargets_GVCF_31Jul2014.sh
* map_recali_perLane_recali_perSample_covOnTargets_26MAY2014.sh
* redupFastuniq_31Jul2014.sh


Other things
-----------------

* Annovar installation and annotations may be very different, so the annotation scripts may need to change a lot

To be continued…
