#!/bin/bash
for i in *.clipped.fastq;
do
	echo working with "$i"
	newfile="$(basename $i .trimmed.clipped.fastq)"
	fastx_collapser -Q33 -v -i $i -o temp.collapsed.txt
	python ../scripts/fastqduplicatecounter.py temp.collapsed.txt temp.collapsed.headers.txt>> ${newfile}.duplicateCount.txt
	rm temp.collapsed.headers.txt
	rm temp.collapsed.txt
done;

