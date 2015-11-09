#!/usr/bin/env python

import urllib
import sys
import os
##### totalannotation.py by DJ Barshis
##### This script takes an input fasta file of sequence names and sequences, and blast results files of blasts against 
##### nr (parsed .txt with 1 hit per line) and swissprot and tremble (in -outfmt 7) uniprot databases 
##### and downloads the corresponding uniprot flat files from the www.uniprot.org web server,
##### extracts particular annotation information from the nr blast and each uniprot flat file and combines it into a meta-annotation table.
##### you will need to create a 2-line .txt file that has the names of the particular columns you would like to extract from the
##### nr parsed blast file separated by tabs (these files can be large so I suggest extracting the header using head or less in terminal
##### the second line consists of the "bad words" you want to skip over in you nr results separated by tabs.
##### I usually use "predicted	PREDICTED	hypothetical	unknown" or some combination thereof.

# usage is totalannotation.py YOUR_contigs.fasta BLASTx2nr.txt nrcolumnheadersandbadwords.txt BLASTx2Sprot.txt BLASTx2TrEMBL.txt evaluethreshold directoryforflatfiles(no slashes) outtablename.txt

#this is for setting how the script sorts your contigs into order
#change the word to 'text' for a text-based sorting or 'coral' for a
#palumbi-lab coral-specific numerical sorting
textorcoralsort = 'coral'

#innames, inseqs read_fasta_lists(sys.argv[1])
#sys.argv[2] = BLASTx2nr.txt
#sys.argv[3] = thingsfornr.txt
#uniprotIDs read_uniprot(sys.argv[4], sys.argv[5]) 
evalue=float(sys.argv[6])
directory=sys.argv[7] #name only, no /'s
#o=open(str(sys.argv[8]), 'w') # New data table file name

#####This reads in a fasta file and extracts the sequence names into a dictionary as the keys

def read_fasta_dict(file):
	fin = open(file, 'r')
	filelines=fin.readlines()
	filelines.append('EOF')
	count=0
	names={}
	seqs=[]
	numseqs=0
	for line in filelines:
		if line=='EOF':
			names[cols[0]]='%i' %(len(seq))
		line=line.strip()
		if line and line[0] == '>':                #indicates the name of the sequence
			if count>=1:
				names[cols[0]]='%i' %(len(seq))
			count+=1
			line=line[1:]
			cols=line.split(' ')
			seq=''
		else: seq +=line
	fin.close()
	return names

innames=read_fasta_dict(sys.argv[1])
print 'Read in fasta of %i sequences: ...' %(len(innames.keys()))

####This function reads in a parsed (every hit on one line) nr blast file and extracts certain columns and returns a dictionary

def nr_dict(file, colstoextract):
	fin = open(file, 'r')       # open input file
	cols2extract = open(colstoextract, 'r')
	d={}
	headers=[]
	contig=''
	linenum=0
	goodnrhits=0
	for line in fin:
		linenum+=1
		line=line.rstrip()
		cols=line.split('\t')
		if linenum == 1:
			headers=line	#Used to copy header to new files
# this loop is for extracting the column indexes for the column names specified on the first line of the stufffornr.txt file
			extractlinecount=0
			for aline in cols2extract:
				extractlinecount+=1
				if extractlinecount==1:
					aline=aline.rstrip()
					words=aline.split('\t')
					hitdescription=cols.index(words[0])
					nrEval=cols.index(words[1])
		if linenum >1:
			cols[0]=cols[0].split(' ')[0]
			if cols[0] == contig:
				d[cols[0]].append('%s\t%s' %(cols[hitdescription],cols[nrEval]))
			else: 
				if float(cols[8]) <= evalue:
					goodnrhits+=1
					contig = cols[0]
					numhit = 1
					d[cols[0]]=d.get(cols[0],[])
					d[cols[0]].append('%s\t%s' %(cols[hitdescription],cols[nrEval]))

	fin.close()
	cols2extract.close()
	return headers, d, goodnrhits

headers, d, goodnrhits=nr_dict(sys.argv[2], sys.argv[3])
print "Read in nr blast..."
print '%s%i' %('Number of good nr matches: ',goodnrhits)
print '%s%i' %('Number not matched in nr: ',len(innames.keys())-goodnrhits)
print "Searching for badwords..."

######This function parses the nr dictionary for hits that do not contain badwords (e.g. 'Predicted', 'hypothetical', etc.)

def parse_badwords(value, badwords):
	onlybad=0
	madegood=0
	badhits=[]
	goodhits=[]
	tophit=value[0]
	for each in value:
		numbadhits=0
		for item in badwords:
			if item in each:
				numbadhits+=1
		if numbadhits >=1:
			badhits.append(each)
		if numbadhits == 0:
			goodhits.append(each)
	if len(goodhits)==0:
		onlybad +=1
	if len(goodhits)>=1:
		madegood +=1
	goodhits+=badhits
	return tophit, goodhits, onlybad, madegood

badwordlist=[]
#reading in a list of badwords from stufffornr.txt
badwordfile=open(sys.argv[3],'r')
badwordline=0

for line in badwordfile:
	badwordline+=1
	if badwordline==2:
		line=line.rstrip()
		badwordlist=line.split('\t')

onlybadnrs=0
madegoodnrs=0

####this step loops through the entrys in your contig dictionary
####and calls the badword parser for each entry that has a match in the nr dictionary and returns the top hit and the top non-badword hit (if there is one)

for key,value in innames.items():
	if d.has_key(key):
		tophit, goodhits, onlybad, madegood= parse_badwords(d[key], badwordlist)
		innames[key]='%s\t%s\t%s' %(innames[key],tophit, goodhits[0])
		onlybadnrs+=onlybad
		madegoodnrs+=madegood
	else:
		innames[key]+='\t%s\t%s\t%s\t%s' %('No_sig_nr_hit','No_sig_nr_hit','No_sig_nr_hit','No_sig_nr_hit')

print '%s%i' %('Number of nr hits with only a bad word hit: ', onlybadnrs)
print '%s%i' %('Number of nr hits with a good word hit: ', madegoodnrs)

#######This function reads in the swissprot and trembl outputs and returns
#######a dictionary that contains the top uniprot ID from swissprot (if available) or trembl (if no swissprot match was found)

def read_uniprot(sprotfile,tremblfile):
	queryname=''
	uniprotIDs={}
	uniqueprotIDs={}
	sprotmatch=0
	tremblpeats=0
	tremblmatch=0
	sprot = open(sprotfile,'r')
	trembl = open(tremblfile,'r')
	for line in sprot:
		line=line.rstrip()
		if line[0] == '#':
			continue
		else:
			cols=line.split('\t')
			if cols[0] == queryname:
				continue
			else:
				if float(cols[10]) <= evalue: #for parsing based on threshold value
					sprotmatch+=1
					ID=cols[1].split('|')
					uniprotIDs[cols[0]]=uniprotIDs.get(cols[0],[])
					uniprotIDs[cols[0]].append(ID[1])
					innames[cols[0]]+='\t%s\t%s\t%s' %(ID[1],cols[2],cols[10])
					queryname=cols[0]
					if uniqueprotIDs.has_key(ID[1]):
						continue
					else:
						uniqueprotIDs[uniprotIDs[cols[0]][0]]=''
	print 'Read in swissprot blast ...'
	print '%s%i' %('Number of good swissprot matches: ', sprotmatch)
	for line in trembl:
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
					if uniprotIDs.has_key(cols[0]):
						uniprotIDs[cols[0]].append(ID[1])
						queryname=cols[0]
						tremblpeats+=1
					else:
						uniprotIDs[cols[0]]=uniprotIDs.get(cols[0],[])
						uniprotIDs[cols[0]].append(ID[1])
						innames[cols[0]]+='\t%s\t%s\t%s' %(ID[1],cols[2],cols[10])
						queryname=cols[0]
						tremblmatch+=1
						if uniqueprotIDs.has_key(uniprotIDs[cols[0]][0]):
							continue
						else:
							uniqueprotIDs[uniprotIDs[cols[0]][0]]=''

	print 'Read in TrEMBL blast ...'
	print '%s%i'%('Number of repeat matches from TrEMBL: ', tremblpeats)
	print '%s%i'%('Number of additional good matches from TrEMBL: ', tremblmatch)
	print '%s%i' %('flatfilesneeded: ',len(uniqueprotIDs.keys()))

	return uniprotIDs, uniqueprotIDs

#this line calls the uniprot reading function
uniprotIDs, uniquesforflats=read_uniprot(sys.argv[4], sys.argv[5])
print 'downloading flat files ...'

#this loop downloads all the uniprot flat files for the list of unique uniprotIDs that was parsed from the blast results
for key, value in uniquesforflats.items():
	if os.path.exists('./'+directory+'/'+key+'.txt'): #thanks JTL for this addition!
		continue
	else:
		urllib.urlretrieve('http://www.uniprot.org/uniprot/'+key+'.txt', './'+directory+'/'+key+'.txt')

print 'extracting relevant info from flat files ...'
print 'don\'t worry this takes awhile ...'

########this function extracts the relevant information from each individual flat file

def extractGO(contigname):
	if uniprotIDs.has_key(contigname):
		flatfile = open('./'+directory+'/'+uniprotIDs[contigname][0]+'.txt','r')
		ID='No_ID'
		DE='No_description'
		description=0
		KEGG='No_KEGG'
		flatfiledict={}
		GOcodes=''
		GOBP=''
		GOMF=''
		GOCC=''
		keywords=''
		for line in flatfile:
			line=line.rstrip()
			if line[0:2] == 'ID':
				line=line.split(' ')
				ID=line[3]
			if line[0:2] == 'DE' and description == 0:
				line=line.split('=')
				DE=line[1][:-1]
				description +=1
			if line[0:2] == 'DR':
				if line[5:9] == 'KEGG':
					line=line.split(';')
					KEGG=line[1].strip()
				if line[5:7] == 'GO':
					line=line.split(';')
					GOcodes+='%s \\\\ ' %(line[1].strip())
					if line[2].strip().split(':')[0] == 'C':
						GOCC+='%s (%s);' %(line[2].strip().split(':')[1], line[1].strip())
					if line[2].strip().split(':')[0] == 'P':
						GOBP+='%s (%s);' %(line[2].strip().split(':')[1], line[1].strip())
					if line[2].strip().split(':')[0] == 'F':
						GOMF+='%s (%s);' %(line[2].strip().split(':')[1], line[1].strip())
			if line[0:2] == 'KW':
				line=line[2:].split(';')
				for item in line:
					if item == '':
						continue
					else:
						keywords+='%s;' %(item.strip())
		if GOcodes=='':
			GOcodes='No_GOcodes'
		if GOBP=='':
			GOBP='No_GOBP'
		if GOMF=='':
			GOMF='No_GOMF'
		if GOCC=='':
			GOCC='No_GOCC'
		if keywords=='':
			keywords='No_keywords'
		outstring='\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(ID, DE, KEGG, GOcodes, GOBP, GOMF, GOCC, keywords)
		nomatch=0
	else:
		nomatch=1
		outstring='\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' %('No_Uniprotmatch','No_%_identity','No_evalue','No_ID','No_Description','No_KEGG','No_GO','No_GOCC','No_GOBP','No_GOMF','No_keywords')
	return outstring, nomatch

notmatched=0
extractingcounter=0
#####This loop calls the extraction function for each contig that has a uniprot match

for key, value in innames.items():
	extractingcounter+=1
	outstring, nomatch = extractGO(key)	
	innames[key]+=outstring
	notmatched+=nomatch

o=open(str(sys.argv[8]), 'w') # New data table file name

o.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %('ContigName', 'ContigLength', 'topnrMatch','topnrEvalue', 'nobadwordnrMatch', 'nobadwordnrEvalue','Uniprotmatch','%_identity','evalue','ID','Description','KEGG','GO','GOCC','GOBP','GOMF','Keywords')) #used if you want a specific header and filename for each file

print '%s%i' %('Hits not matched in sprot: ', notmatched)
print 'compiling extracted information ...'

############this if for sorting your contigs based on text order#############
if textorcoralsort == 'text':
	l=[]
	for key,value in innames.items():
		l.append((key,value))
	l.sort()
	for item in l:
		o.write('%s\t%s\n' % (item[0], item[1])) #writes each line of the tuple as separate tab delimited text
		
	o.close()

#############this is for sorting your contigs based on our coral specific contig names##############
if textorcoralsort == 'coral':
	l=[]
	joinedcontigcounter=600247
	
	for key,value in innames.items():
		name=key.split(' ')
		if name[0][0:6]=='contig':
			newname=name[0].split('_')
			if len(newname)==1:
				num=int(newname[0][6:])
			if len(newname)>1:
				joinedcontigcounter+=1
				num=joinedcontigcounter
		if name[0][0:6]=='c_sym_':
			newname=name[0].split('_')
			num=700000+int(newname[2])
		if name[0][0:6]=='d_sym_':
			newname=name[0].split('_')
			num=900000+int(newname[2])
				
		l.append((num,key,value))
	l.sort()
	for item in l:
	#	print item
		o.write('%s\t%s\n' % (item[1], item[2])) #writes each line of the tuple as separate tab delimited text
		
	o.close()
	
