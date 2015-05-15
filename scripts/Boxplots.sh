#!/bin/bash

for i in *.qualstats.txt;
do
	echo working with "$i"
	newfile="$(basename $i .qualstats.txt)"
	bash fastq_quality_boxplot_graph.sh -i $i -t $i -o ${newfile}.Box.png
	bash fastx_nucleotide_distribution_graph.sh -i $i -t $i -o ${newfile}.Nuc.png
done;
convert *.png Boxplots.pdf
