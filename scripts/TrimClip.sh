#!/bin/bash
for i in *.fastq;
do
	echo working with "$i"
	newfile="$(basename $i .fastq)"
	fastq_quality_trimmer -Q33 -v -t 20 - l 20 -i $i |fastx_clipper -Q33 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG - l 20 -n -v > ${newfile}.trimmed.clipped.fastq 
done;
