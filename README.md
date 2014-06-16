ExomePipe_NCL_IGM
=================

Intro
-----------------

This is a tool kit from exome-seq data analysis. It utilizes a HPC managed by Sun Grid Engine and foucus on germline variants detection exclusively for now. The scripts used in each stage of the analysis are grouped in an folder. 

There are two types of scripts right now:
*The scripts doing the work, and
*The scripts control workflow and sumbit ordered scripts to the queue engine of an HPC.

Most of the time, the 1st kind of scripts requires no change from project to project, but the 2nd type of scrpts are required to change in a new project.

Ideally, each sample should have it is own folder under a project folder, fastq files for each sample should be under their sample folder separatedly. 


