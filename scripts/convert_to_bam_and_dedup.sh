#!/bin/bash
# This scripts converts SAM files to BAM, sorts the BAM, marks and removes duplicate reads.

#Convert SAM file to .bam | sort .bam with Samtools and Mark Duplicates with Picard Tools:

for i in *.sam;
do
	newfile="$(basename $i .sam)"
	samtools view -bS $i | samtools sort - $i.sorted
	java -Xmx2g -jar ../scripts/MarkDuplicates.jar INPUT=$i.sorted.bam OUTPUT=${newfile}.dedup.bam  METRICS_FILE=${newfile}.metricsfile MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=True
done;
