#!/bin/bash

#this script is for bwa mapping fastq files to a reference library.

#Indexing reference library (the dictionary and the .fai index files will be needed tomorrow):  

bwa index -p assembly -a is ./assembly/Trinity.fasta

Dictionary="./assembly/Trinity.dict"

if [ -f  "$Dictionary" ]; 
then 
	echo "Trinity.dict already created"
else
	echo "file does not exist"
 	java -jar ../scripts/CreateSequenceDictionary.jar R= ./assembly/Trinity.fasta O= ./assembly/Trinity.dict
fi

samtools faidx ./assembly/Trinity.fasta

#Mapping and converting to SAM file format:

#For Single End lanes:
#Add your personal alignment parameters to each line before executing the script. We will leave -t (number of threads) to 2.


for i in *.clipped.fastq;
do
	echo working with "$i"
	newfile="$(basename $i .trimmed.clipped.fastq)"
	bwa aln -n .01 -k 5 -l 30 -t 2 assembly $i |bwa samse assembly - $i | ../scripts/AddReadGroup.pl -i - -o ${newfile}.sam -r ${newfile} -s ${newfile} -l ${newfile} -p Illumina


done;

#bwa samse assembly FR32.sai FR32_ATCACG_before.clipped.fastq >  FR32.sam

#
#For Paired-End lanes:
# for paired-end lanes there are four lines per set of two paired-end and one remaining singles .fastq files
# so you'll need to copy the following four lines for each set of paired-end files:

#Then replace ABBREVIATEDREFERENCEFILENAME with your abbreviated reference name from above,
#and SAMPLENAME with your sample name, e.g. "HMS_Cold_1".

#bwa aln -n 0.005 -k 5 -t 4 ABBREVIATEDREFERENCEFILENAME SAMPLENAME_trimmed_clipped_singles.fastq | bwa samse -r '@RG\tID:SAMPLENAME\tSM:SAMPLENAME\tPL:Illumina' ABBR$
