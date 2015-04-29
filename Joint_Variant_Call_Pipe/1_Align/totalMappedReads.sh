#!/bin/bash

Argus="$@"

module load apps/samtools

for BAM in $Argus; do
	sizeeee=$(samtools view -c -F 256 $BAM)
	echo "$BAM #total reads: $sizeeee"
done
