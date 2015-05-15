#!/usr/bin/env python
import sys

#################################################################################################
#                                                                                               #
#   This program creates all necc. input files for smartpca analysis                            #
#   Usage: vcf2smartpca.py in.vcf prefix min_gen_quality_threshold                              #
#   min_gen_quality_threshold is for the individual genotypes, not the SNP as a whole           #
#   Author: Jason Ladner                                                                        #
#                                                                                               #
#################################################################################################

'''---------------------------> VCF Format <-----------------------------
    cols[0] = contig name
    cols[1] = SNP position
    cols[2] = ID???? - empty in my vcf files
    cols[3] = Reference base
    cols[4] = Alternative base
    cols[5] = SNP quality (Phred score)
    cols[6] = filter flags
    cols[7] = SNP Info
    cols[8] = genotype format (e.g., GT:AD:DP:GQ:PL)
    cols[9:] = individual genotypes'''

#-------------------------------------------------------------#
#  This function takes a file and creates a dictionary with:  #
#    -one entry for each line                                 #
#    -keys will be first column of each line                  #
#    -values will be the complete line                        #
#-------------------------------------------------------------#
def file_dict(file):
    fin = open(file, 'r')       # open input file
    d={}
    
    for line in fin:
        line=line.rstrip()
        cols=line.split()
        d[cols[0]]=line
    fin.close()
    return d

#---------------------------------------------------------------------------------#
#  This function creates two lists, which will be input to smartpca_input()    :  #
#    -'headers' should be a list, headers[4:] should be individual names          #
#    -'snp_genos' should be a list with one list per indiviual                    #
#        -Each indiv list will contain genotpe data for each SNP                  #
#    -'pops' is a dictionary with indiv names as the keys and their population    #
#      of origin as the values                                                    #
#---------------------------------------------------------------------------------#
def genos4smart(headers, snp_gens, pops):
    loci_names=[]
    for index, contig in enumerate(contigs):
        loci_names.append(contig + '_' + pos[index])

    lines=[]
    indivcount=0
    for indiv in headers[4:]:
        thisline=[indiv, pops[indiv]]
        for each in snp_gens[indivcount]:
            for hap in each:
                thisline.append(hap)
        indivcount+=1
        lines.append(thisline)
    return loci_names, lines

#---------------------------------------------------------------------------------#
#  This function creates all the necc. input files smartpca:                      #
#     -'loci_names' should be a list that contains the names of all the SNPs      #
#     -'lines' is a list that contains an entry for each indivdual                #
#      -Each entry is also a list, with an entry for the indiv's ID & population  #
#      as well as two entries for each locus saying the # of each allele          #
#                                                                                 #
#---------------------------------------------------------------------------------#
def smartpca_input(loci_names, lines, prefix, omit=[]):
#Opens an output text files with prefix as specified by user
    GENO = open(prefix + '_Geno', 'w')
    INDIV = open(prefix + '_Indiv', 'w')
    SNP = open(prefix + '_SNP', 'w')
    print "omitting: " + str(omit)
    pops=[]
    names=[]
    
    print '# of loci = %d' % (len(loci_names))
    
    loci=[]
    for locus in loci_names:
        loci.append([])
    
    for indiv in lines:
        
        if indiv[0] in omit: continue
        
        INDIV.write('%s\tU\t%s\n' % (indiv[0], indiv[1]))               #Writes an entry in the INDIV file for each individual: 'ID\tU\tPopulation'
        
        names.append(indiv[0])              #Adds each individual's name to the list names
        pops.append(indiv[1])               #Adds each individuals population to the list pops
        
        i=0
        k=0
        while i<len(indiv)-2:                              #the number of columns of actual data
            loci[k].append([indiv[i+2], indiv[i+3]])
            i+=2
            k+=1

#!!!!!!!!!!Can not include data for individuals with only info from one allele, because each snp is only coded with # of ref alleles!!!!!!!!!!!!
    for index, locus in enumerate(loci):
        indcount=0
        SNP.write(loci_names[index] + '\t1\t0.0\t1\n')                             #All SNPs specified to be in the same place on same chrom, I don't think it matters though because I do not utilize this data
        for ind in locus:
            if ind==['.', '.']:                                                     #Identifies individuals without genotype data and excludes them from the GENO file for this locus
                indcount+=1
                continue
            else:
                GENO.write('%s\t%s\t%d\n' % (loci_names[index], names[indcount], ind.count('0')))
                indcount+=1

    PAR = open('par.' + prefix, 'w')
    PAR.write('genotypename:    %s_Geno\nsnpname:         %s_SNP\nindivname:       %s_Indiv\nevecoutname:     %s_%dloci.evec\nevaloutname:     %s_%dloci.eval\nnumoutevec:      8\nnumoutlieriter:  3\nusenorm:         YES\nsnpweightoutname: %s_snpweights' % (prefix, prefix, prefix, prefix, len(loci_names), prefix, len(loci_names), prefix))
    PAR.close()
    GENO.close()
    SNP.close()
    INDIV.close()


#-------> The Actual Program Starts Here <-------

pops={'FR32_ATCACG_before': '1', 'FR33_TTAGGC_before': '1', 'FR34_ACTTGA_before': '1', 'FR35_GATCAG_before': '1', 'FR36_TAGCTT_before': '1', 'FR37_GGCTAC_before': '1', 'FR77_ATCACG_after': '2', 'FR78_TTAGGC_after': '2', 'FR79_ACTTGA_after': '2', 'FR80_GATCAG_after': '2', 'FR81_TAGCTT_after': '2', 'FR82_GGCTAC_after': '2'}

#Opens an infile specified by the user. Should be a .vcf file 
IN = open(sys.argv[1], 'r')
quality_cutoff=float(sys.argv[3])                    #for individual genotypes (not the SNP as a whole)

linecount=0
contigs=[]
pos=[]
ref=[]
alt=[]
all_gens=[]
locus=[]
headers=[]

#Reads through the vcf file line by line
for line in IN:
    line=line.rstrip()
    cols=line.split('\t')
    gens=[]                           #Will hold individual level genotype data, resets list form previous SNP
    gencount=0                        #Will be used below to keep track of the individuals when appending new SNPs
    
    if cols[0][0] == '#' and cols[0][1] != '#':       #To pull out headers

        headers=cols[0:2] + cols[3:5] + cols[9:]      #Creates list of headers of interest
        print headers[4:]                             #Prints the IDs of the individuals with genotypes in the vcf

    if cols[0][0] != '#':                             #Specifies only SNP lines
        linecount+=1                                  #Keeps track of the numbers of SNPs
        
        gens=cols[9:]                  #Creates a list of all the genotype info
        desc=cols[8]
        parts=desc.split(':')			#Finds which value in the genotype info is Genotype quality
        count=0							#and stores the position of the GQ.
        for items in parts:
        	if items == 'GQ':
        		GQ_pos = count
        	count+=1        
        
        if linecount==1:                      #will only execute these commands once, when reading the first line of the file
            for gen in gens:
                all_gens.append([])           #Creates seperate list entry for each individual, will hold genotypes at all SNPs
                locus.append([])              #Creates seperate list entry for each individual, will hold genoypes at a specific locus 
        
        for gen in gens:                    #Steps through the genotype info for each individual
            parts=gen.split(':')            #Breaks apart the different sections of the genotype info
            genotype=parts[0]               #Most probable genotype for the individual, with the two alleles seperated by '/'
            genotypes=genotype.split('/')
            
            if genotypes != ['.','.']:         #Had to be added because in new version of GATK, .vcf files are written so that if there is no genotype info for an individual, their genotype is './.' as opposed to './.:0:0:0:0', or whatever
                gen_quality=float(parts[GQ_pos])
                
                if gen_quality<quality_cutoff:
                    genotypes=['.','.']
                
            locus[gencount]=genotypes           #Replace genotype data from last locus with genotype form this locus
            gencount+=1                         #Keeping track of the individual whose data is being examined
    
        if ['.', '.'] not in locus and locus != []:
            ##########################################################################
            #   Use these 4 lines to only output info for SNPs with no missing data  #
            ##########################################################################
            contigs.append(cols[0])   #Creates a list of all the contig names               
            pos.append(cols[1])       #Creates a list of all the SNP base positions
            ref.append(cols[3])       #Creates a list of all the reference bases
            alt.append(cols[4])       #Creates a list of all the alternative bases
            locuscount=0
            for each in locus:
                all_gens[locuscount].append(each)                  #Appends genotype data to the all_gens list
                locuscount+=1                                      #keeps track of the number of loci with no missing data

print "Number of SNPs in reference contigs = %d" % (linecount)

list_headers, list_lines = genos4smart(headers, all_gens, pops)
smartpca_input(list_headers, list_lines, sys.argv[2])