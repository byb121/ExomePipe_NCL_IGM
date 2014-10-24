#!/bin/bash

Argus="$@"

module load apps/samtools

for BAM in $Argus; do
	echo samtools index $BAM
	samtools index $BAM
done
