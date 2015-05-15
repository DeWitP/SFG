#!/usr/bin/env python

#This scripts takes a collapsed fasta file (from fastx_collapser) and calculates the fractions of Unique reads, duplicate reads and singletons.
#It then prints out the info to a tab-delimited text file that can be opened in Excel.
import sys

#Opens a fasta file with the database sequences as specified by user. Must be non-interleaved.
SEQS = open(sys.argv[1], 'r')

#Opens an output text file as specified by user
OUT = open(sys.argv[2], 'w')
duplicates=0
dupreadcount=0
singletons=0
totalreadcount=0
OUT.write('UniqueNumber'+'\t'+'Numseqscollapsed')
for line in SEQS:
	line=line.rstrip()
	if line[0:1] == '>':
		line=line.strip('>')
		cols=line.split('-')
		if int(cols[1]) > 1:
			duplicates+=1
			dupreadcount+=int(cols[1])
			totalreadcount+=int(cols[1])
			OUT.write('\n'+'\t'.join(cols))
		elif int(cols[1]) == 1:	
			totalreadcount+=int(cols[1])
			singletons+=1

percentsingletons=float(singletons)/float(totalreadcount)
percentduplicates=float(dupreadcount)/float(totalreadcount)
percentuniquedups=float(duplicates)/float(dupreadcount)
print 'Number of unique reads with duplicates'+'\t'+'Number of collapsed duplicates'+'\t'+ 'Number of singletons' + '\t' + 'Total number of reads' + '\t' + 'Percent singletons' + '\t' + 'Percent duplicate reads' + '\t' + 'Percent unique reads with duplicates'
print '\n' + str(duplicates) + '\t' + str(dupreadcount) + '\t' + str(singletons) + '\t' + str(totalreadcount) + '\t' + str(percentsingletons) + '\t' + str(percentduplicates) + '\t' + str(percentuniquedups)
OUT.close()