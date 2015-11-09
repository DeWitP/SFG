#!/bin/bash

parse_blastn.py Blastx2nroutput.xml parsedBlastx2nr.txt
#if parse_blast.py does not work, and alternative parser is BlastParse.pl, though you'll have to change the evalue column in the ReParseBlastbycutoffs.py below

#the evalue is in column 12 if you're using the parse_blast.py or column 8 if you're using BlastParse.pl in the first step above
ReParseBlastbycutoffs.py parsedBlastx2nr.txt columntoparse evaluecutoff tophitreParsedBlastx2nr.txt