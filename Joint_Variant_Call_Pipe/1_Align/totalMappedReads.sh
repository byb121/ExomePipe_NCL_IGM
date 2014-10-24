#!/bin/bash

Argus="$@"

module load apps/samtools

for BAM in $Argus; do
	sizeeee=$(samtools idxstats $BAM | awk '{ sum += $3 } END { print sum }')
	echo "$BAM #total reads: $sizeeee"
done
