#!/usr/bin/env python

import sys
from Bio.Blast import NCBIXML
#Usage, opens an outfile and then parses any number of .xml files into that outfile, printing all hits
#parse_blastn.py outfile.txt anynumberofinfiles.xml
OUT = open(sys.argv[1], 'w')

OUT.write("Query Name\tQuery Length\tSubject Name\tSubject Length\tAlignment Length\tQuery Start\tQuery End\tSubject Start\tSubject End\tQuery Sequence\tSubject Sequence\tHsp Score\tHsp Expect\tHsp Identities\tPercent Match\tNumber_of_gaps")
for xml_file in sys.argv[2:]:
	result_handle = open(xml_file)
	blast_records = NCBIXML.parse(result_handle)
	for rec in blast_records:
		for alignment in rec.alignments:
				for hsp in alignment.hsps:
					OUT.write('\n'+ str(rec.query) + '\t' + str(rec.query_length) + '\t' + str(alignment.title) + '\t' + str(alignment.length) + '\t' + str(hsp.align_length) + '\t' + str(hsp.query_start) + '\t' + str(hsp.query_end) + '\t' + str(hsp.sbjct_start) + '\t' + str(hsp.sbjct_end) + '\t' + str(hsp.query) + '\t' + str(hsp.sbjct) + '\t' + str(hsp.score) + '\t' + str(hsp.expect) + '\t' + str(hsp.identities) + '\t' + str(float(hsp.identities)/int(hsp.align_length)) + '\t' + str(hsp.gaps))

