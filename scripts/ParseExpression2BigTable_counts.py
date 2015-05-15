#!/usr/bin/env python

# Written by Daniel Barshis, June 2011

import sys
#sys.argv[1] Input file is a list of contig names in a single column with ContigName as the column header
#sys.argv[2] Output file name
#sys.argv[3:] Any number of files to add columns from


#makes a dictionary of a file where the first column is a list of contignames and subsequent columns with any associated information (i.e. counts, blast hits, etc.)
def make_dict1(file):
	fin = open(file, 'r')
	dict={}
	headers=[]
	count=0
	for line in fin:
		count+=1
		line=line.rstrip()
		cols=line.split('\t') #for tab-delimited text files
		if count==1:
			headers=cols[0:]
			#count+=1
		if count > 1:
			dict[cols[0]]=cols[1:]
	
	return dict, headers
		
dictbase, dictheaders=make_dict1(sys.argv[1]) #Input data table
xpressionfiles=sys.argv[3:] #Any number of expression files
	
#Loops through a list of secondary files that you would like to add information from
#the purpose would be to combine information from multiple files into one big "meta" table

xpressfiles=[]
xtrahits=[]
for file in xpressionfiles:
	Xvalues, Xheaders=make_dict1(file)
	xpressfiles.append(file+'_'+str(Xheaders[1])) # used to denote the column that you're extracting
	xtracount=0
#	For adding info when Xvalues is missing values that are in your dictbase
#  	for item in dictbase.keys():
#  		if Xvalues.has_key(item):
# #			dictbase[item]+=Xvalues[item] #appends a list of strings
#  			dictbase[item].append(Xvalues[item][2])
#  		else:
#  			xtracount+=1
#  			dictbase[item].append('No Sig Blast Hit')
#			dictbase[item].append('\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % ('NA','NA','NA','NA','NA','NA','NA',))
# 	For adding info when Xvalues has all the values in your dictbase or more
 	for item in Xvalues.keys():
		if dictbase.has_key(item):
#			dictbase[item]+=Xvalues[item] #appends a list of strings
			dictbase[item].append(Xvalues[item][0]) #appends a single column of data as a string appends quality, single counts
 		else:
			xtracount+=1
 	xtrahits.append(file+'='+str(xtracount))
#	print dictbase
o=open(str(sys.argv[2]), 'w') # New data table file name
#o.write('\t'.join(dictheaders)+'\t'+'\t'.join(Xheaders)+'\n') #used if you want all new headers from one single file
o.write('\t'.join(dictheaders)+'\t'+'\t'.join(xpressfiles)+'\n') #used if you want a specific header and filename for each file
print 'Hits not matched' + '\t'.join(xtrahits)
l=[]
for key,value in dictbase.items():
	l.append((key,value))   #translates dictbase into a tuple and adds just the contig # (minus"contig") as the first item
l.sort() #sorts tuple into numeric order
for item in l:
#	print item
	o.write(str(item[0])+'\t'+'\t'.join(item[1])+'\n') #writes each line of the tuple as separate tab delimited text

o.close()
