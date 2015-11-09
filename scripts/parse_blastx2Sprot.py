#!/usr/bin/env python

import sys

#parse_blastx2Sprot by DJ Barshis
#this script takes an -7 text outformat blast output from a blastx against sprot
#and extracts the top uniprot ID for each hit

#usage is parse_blastx2Sprot.py sprotblastoutput.txt outfilename.txt

INFILE = open(sys.argv[1], 'r') #Input file 1 is your blastresults.xml file

evalue = float(sys.argv[2]) #threshold evalue

OUT = open(sys.argv[3], 'w')

OUT.write("QueryName\tUniprotID")

queryname=''

for line in INFILE:
	line=line.rstrip()
	if line[0] == '#':
		continue
	else:
		cols=line.split('\t')
		if cols[0] == queryname:
			continue
		else:
			if float(cols[10]) <= evalue: #for parsing based on threshold value
				ID=cols[1].split('|')
				OUT.write('\n'+cols[0]+'\t'+ID[1])
				queryname=cols[0]

