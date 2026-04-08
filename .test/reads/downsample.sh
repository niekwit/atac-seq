#!/bin/bash

# Get first 1M reads from each fastq.gz file
FASTQ=$(ls large/*.fastq.gz)

for FILE in $FASTQ; do
	zcat "$FILE" | head -n 4000000 | gzip > "$(basename $FILE .fastq.gz).fastq.gz"
done
