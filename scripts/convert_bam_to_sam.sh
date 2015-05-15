#!/bin/bash

# This scripts converts BAM files with duplicate reads removed to SAM files for differential gene expression analysis.

for i in *.dedup.bam;
do
        newfile="$(basename $i .bam)"
        samtools view -h $i > ${newfile}.sam
done;
