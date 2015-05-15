#!/usr/bin/env python

#########################################################################################################################################
#                                                                                                                                       #
#   Written by Jason Ladner                                                                                                             #
#   questions?: jtladner@gmail.com                                                                                                      #
#   Usage:                                                                                                                              #
#   python getgenosfromvcf_clean.py input_vcf_file  output_file_name  output_type  genotype_quality_cutoff  subset_contig_name_file     #
#                                                                                                                                       #
#   input_vcf_file = the .vcf file that you want to extract genotype information from                                                   #
#   output_file_name = the name of the file that the genotype info will be written to                                                   #
#                      it will be created if it doesn't exist, or written over if it does (*so be careful*)                             #
#   output_type = either 'cols' or 'rows' (enter on the command line without the quotes)                                                #
#                       'cols' specifies that SNPs should be in columns with 2 columns per SNP (one for each allele)                    #
#                              This is very similar to the input required by smartpca                                                   #
#                       'rows' specifies that SNPs should be in rows (one for each SNP),2 columns for each indiv (one for each allele)  #
#                              This is easier to import into excel to visualize the data                                                #
#   genotype_quality_cutoff = a floating point number that specifies the minimum quality allowed by genotype information to be used     #
#                             This parameter is optional (default = 0), but must be specified if subset_contig_name_file is specified   #
#   subset_contig_name_file = An optional file that if used should be a list of contig names, with one contig name per line             #
#                             This parameter can be used to specifiy a subset of contigs to extract genotype information from           #
#                                                                                                                                       #
#########################################################################################################################################


'''---------------------------> VCF Format <-----------------------------
    cols[0] = contig name
    cols[1] = SNP position
    cols[2] = ID???? - empty in my vcf files
    cols[3] = Reference base
    cols[4] = Alternative base
    cols[5] = SNP quality (Phred score)
    cols[6] = filter flags
    cols[7] = SNP Info
    cols[8] = genotype format (e.g., GT:AD:DP:GQ:PL)             ### This is different for some!!!!!!!!   It can even be different for different SNPs in the same .vcf
    cols[9:] = individual genotypes'''

import sys

#-------------------------------------------------------------#
#  This function takes a file and creates a dictionary with:  #
#    -one entry for each line                                 #
#    -keys will be first column of each line                  #
#    -values will be the complete line                        #
#-------------------------------------------------------------#
def file_dict(file):
    fin = open(file, 'r')       
    d={}
    
    for line in fin:
        line=line.rstrip()
        cols=line.split()
        d[cols[0]]=line
    fin.close()
    return d


#-------------------------------------------------------------#
#  This function takes the genotype description for a SNP     #
#   (e.g., GT:AD:DP:GQ:PL)                                    #
#  and returns the position of the genotype quality (GQ)      #
#  with the first position being 0                            #
#-------------------------------------------------------------#
def get_GQ_pos(genotype_description):            
    return genotype_description.split(':').index('GQ')

#-------------------------------------------------------------#
#  This function writes genotype data to an outfile           #
#   with 2 columns per SNP (one for each allele)              #
#                                                             #
#  This is very similar to the input required by smartpca     #
#-------------------------------------------------------------#
def write_snps_as_cols(outfile_name, headers, contigs, pos, ref, alt, all_gens):
    #Opens an output text file as specified by user
    OUT = open(outfile_name, 'w')
    
    OUT.write(headers[0] + '\t')
    for index in range(len(contigs)):                                        
        OUT.write("%s\t\t" % (contigs[index]))                                      #Can be easily changed so that the contig name is a combination of the contig_pos

    OUT.write('\n%s\t' % (headers[1]))         #Positions
    for index in range(len(contigs)):
        OUT.write(pos[index] + '\t\t')        
    OUT.write('\n%s\t' % (headers[2]))         #Reference Alleles
    for index in range(len(contigs)):
        OUT.write("%s\t\t" % (ref[index]))
    OUT.write('\n%s\t' % (headers[3]))         #Alternate alleles
    for index in range(len(contigs)):
        OUT.write('%s\t\t' % (alt[index]))
            
    indivcount=0
    for indiv in headers[4:]:
        OUT.write('\n%s' % (indiv))
        for index in range(len(contigs)):
            for each in all_gens[indivcount][index]:
                OUT.write('\t%s' % (each))
        indivcount+=1
    
    OUT.close()

#-------------------------------------------------------------#
#  This function writes genotype data to an outfile           #
#   with SNPs in rows (one for each SNP)                      #
#  2 columns for each individual (one for each allele         #
#  This is easier to import into excel to visualize the data  #
#-------------------------------------------------------------#
def write_snps_as_rows(outfile_name, headers, contigs, pos, ref, alt, all_gens):
    #Opens an output text file as specified by user
    OUT = open(outfile_name, 'w')

    for item in headers[:4]:
        OUT.write('%s\t' % (item))
    for item in headers[4:]:
        OUT.write('%s\t' % (item))
        
    for index in range(len(contigs)):
        OUT.write('\n%s\t%s\t%s\t%s' % (contigs[index], pos[index], ref[index], alt[index]))
        for ind in all_gens:
            OUT.write('\t%s%s' % (ind[index][0], ind[index][1]))
        
def extract_genos_from_vcf(vcf_file, outfile_name, output_type, quality_cutoff = 0, subset_contig_names = {"no_refs_given": ""}):
    #Opens an infile specified by the user. Should be a .vcf file 
    fin = open(sys.argv[1], 'r')
    if type(subset_contig_names) == type(""):                                                         #Checks to see if a list of names was provided
        subset_contig_names = file_dict(subset_contig_names)                      #Should be a list of the contigs of interest, with one name on each line

    linecount=0
    contigs=[]
    pos=[]
    ref=[]
    alt=[]
    all_gens=[]
    locus=[]
    headers=[]
    
    #Reads through the vcf file line by line
    for line in fin:
        cols=line.rstrip().split('\t')
        gens=[]                           #Will hold individual level genotype data, resets list form previous SNP
        gencount=0                        #Will be used below to keep track of the individuals when appending new SNPs
        
        if cols[0][0] == '#' and cols[0][1] != '#':       #To pull out headers
    
            headers=cols[0:2] + cols[3:5] + cols[9:]      #Creates list of headers of interest
            print "Indivs with genotypes in vcf file: %s" % ("\t".join(headers[4:]))                             #Prints the IDs of the individuals with genotypes in the vcf
    
        if cols[0][0] != '#':                             #Specifies only SNP lines
            if subset_contig_names==None :
                linecount+=1                                  #Keeps track of the numbers of SNPs
                
                contigs.append(cols[0])   #Creates a list of all the contig names
                pos.append(cols[1])       #Creates a list of all the SNP base positions
                ref.append(cols[3])       #Creates a list of all the reference bases
                alt.append(cols[4])       #Creates a list of all the alternative bases

                genotype_description=cols[8]
                GQ_position = get_GQ_pos(genotype_description)
                
                gens=cols[9:]                  #Creates a list of all the genotype info
                
                if linecount==1:                      #will only execute these commands once, when reading the first line of the file
                    for gen in gens:
                        all_gens.append([])           #Creates seperate list entry for each individual, will hold genotypes at all SNPs
                        locus.append([])              #Creates seperate list entry for each individual, will hold genoypes at a specific locus 
                
                for gen in gens:                    #Steps through the genotype info for each individual
                    parts=gen.split(':')            #Breaks apart the different sections of the genotype info
                    genotype=parts[0]               #Most probable genotype for the individual, with the two alleles seperated by '/'
                    genotypes=genotype.split('/')
                    if genotypes != ['.','.']:         #Had to be added because in new version of GATK, .vcf files are written so that if there is no genotype info for an individual, their genotype is './.' as opposed to './.:0:0:0:0', or whatever
                        gen_quality=float(parts[GQ_position])    #May need to change based on the structure of your genotype data
                        
                        if gen_quality<float(quality_cutoff):
                            genotypes=['.','.']
        
                    locus[gencount]=genotypes           #Replace genotype data from last locus with genotype form this locus
                    gencount+=1                         #Keeping track of the individual whose data is being examined
            
                locuscount=0
                for each in locus:
                    all_gens[locuscount].append(each)                  #Appends genotype data to the all_gens list
                    locuscount+=1                                      #keeps track of the number of loci with no missing data

    #--------------->'cols-style output is what you need to go into smartpca, but it still needs to be *****tweaked to include population*****
    if output_type=='cols':    #This will have the contig name, pos, ref and alt alleles listed as rows at the top and then each individual and all of its genotyoes listed in the rows that follow. There will be two columns of data for each SNP
    
        write_snps_as_cols(outfile_name, headers, contigs, pos,  ref, alt, all_gens)
    
    if output_type=='rows':               #This format has a row for each snp, it is easier to view and put into excel
    
        write_snps_as_rows(outfile_name, headers, contigs, pos,  ref, alt, all_gens)
    

###------->>>

if len(sys.argv)==4:             #No quality cutoff or cubset of contig names provided
    extract_genos_from_vcf(sys.argv[1], sys.argv[2], sys.argv[3])

elif len(sys.argv)==5:
    extract_genos_from_vcf(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], subset_contig_names = None)

elif len(sys.argv)==6:
    extract_genos_from_vcf(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

else:
    if len(sys.argv) <4: print "You did not enter enough parameters"
    elif len(sys.argv) >6: print "You entered too many parameters"
