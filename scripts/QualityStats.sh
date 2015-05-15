#!/bin/bash

for i in *.clipped.fastq;
do
	echo working with "$i"
	newfile="$(basename $i .trimmed.clipped.fastq)"
	fastx_quality_stats -Q33 -i $i -o ${newfile}.qualstats.txt
done;

