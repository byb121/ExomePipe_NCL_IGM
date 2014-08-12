#!/bin/bash

Argus="$@"

module load apps/samtools/0.1.18/gcc-4.4.6

for BAM in $Argus; do
	sizeeee=$(samtools idxstats $BAM | awk '{ sum += $3 } END { print sum }')
	echo "$BAM #total reads: $sizeeee"
done
