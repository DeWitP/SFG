#!/usr/bin/env python
import sys

#Written by Jason Ladner

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

def make_bayescan_input(vcf_file, gen_qual_cutoff, indivs2pops_dict, to_make_genofile=0):
    fin = open(vcf_file, "r")
    fout_snpkey = open("snpkey.txt", "w")
    if int(to_make_genofile): fout_genofile = open("used_snp_genos.txt", "w")
    #Reads through the vcf file line by line
    snpcount=0
    full_snpcount=0
    indiv_names = []
    output_lines_bypop = {}
    locus = []
        
    for line in fin:
        cols=line.rstrip().split('\t')
        genos=[]                                                             #Will hold individual level genotype data, resets list form previous SNP
        gencount=0                                                          #Will be used below to keep track of the individuals when appending new SNPs
        
        if cols[0][0] == '#' and cols[0][1] != '#':                         #To pull out headers
    
            headers=cols[0:2] + cols[3:5] + cols[9:]                        #Creates list of headers of interest
            indiv_names = headers[4:]                                       #Prints the IDs of the individuals with genotypes in the vcf
            
            if int(to_make_genofile): fout_genofile.write("%s\t%s\n" % (recursive_join(headers[:4]), recursive_join(headers[4:], "\t\t")))
            
        if cols[0][0] != '#':                                               #Specifies only SNP lines
            snpcount+=1                                                     #Keeps track of the numbers of SNPs
            genos=cols[9:]                                                   #Creates a list of all the genotype info
            desc=cols[8]
            if snpcount==1:                                                 #will only execute these commands once, when reading the first line of the file
                for gen in genos:                                            # This step is needed so that the script can handle vcf files with variable numbers of individuals
                    locus.append([])                                        #Creates seperate list entry for each individual, will hold genoypes at a specific locus 
            locus_genos = get_genos(genos, locus, desc, gen_qual_cutoff)
       
            merged_genos_dict = merge_genos_by_pop(locus_genos, indiv_names, indivs2pops_dict)
            if check_for_min_data(merged_genos_dict):                       # Checks to make sure that each population has at least the minimum amount of data allowed
                full_snpcount+=1
                fout_snpkey.write("%d\t%s_%s\n" % (full_snpcount, cols[0], cols[1]))
                for key, value in merged_genos_dict.items():
                    output_lines_bypop[key] = output_lines_bypop.get(key, [])
                    output_lines_bypop[key].append("%d\t%d\t2\t%d\t%d" % (full_snpcount, value.count('0') + value.count('1'), value.count('0'), value.count('1')))
                    
                    if int(to_make_genofile): fout_genofile.write("%s\t%s\t%s\t%s\t%s\n" % (cols[0], cols[1], cols[3], cols[4],recursive_join(locus_genos)))
                    
                    #Need to add something here that will print the genotype information for each individual to an output file, if wanted. This is what Steve wants for excel
    fin.close()
    fout_snpkey.close()
    write_bayes_input(output_lines_bypop, "bayes_input.txt")

def get_genos(full_genos_data_list, locus, desc, quality_cutoff):
    parts=desc.split(':')			#Finds which value in the genotype info is Genotype quality
    count=0							#and stores the position of the GQ.
    for items in parts:
        if items == 'GQ':
            GQ_pos = count
        count+=1   
    gencount=0
    for gen in full_genos_data_list:                                #Steps through the genotype info for each individual
        parts=gen.split(':')                                        #Breaks apart the different sections of the genotype info
        genotype=parts[0]                                           #Most probable genotype for the individual, with the two alleles seperated by '/'
        genotypes=genotype.split('/')
        if genotypes != ['.','.']:                                  #Had to be added because in new version of GATK, .vcf files are written so that if there is no genotype info for an individual, their genotype is './.' as opposed to './.:0:0:0:0', or whatever
            # The exact positions of these variables within a .vcf can change
            # Check your genotype format and change script accordingly
            gen_quality=float(parts[GQ_pos])
            
            if gen_quality<quality_cutoff:
                genotypes=['.','.']

        locus[gencount]=genotypes                                   #Replace genotype data from last locus with genotype from this locus
        gencount+=1                                                 #Keeping track of the individual whose data is being examined

    return locus

def merge_genos_by_pop(locus_genos, indiv_names, indivs2pops_dict):
    genosbypop_dict={}
    for index, value in enumerate(locus_genos):
        if indivs2pops_dict.has_key(indiv_names[index]):
            pop_of_origin = indivs2pops_dict[indiv_names[index]]
            genosbypop_dict[pop_of_origin] = genosbypop_dict.get(pop_of_origin, []) + value
    return genosbypop_dict

def write_bayes_input(output_lines_bypop, outfile):
    fout = open(outfile, "w")
    fout.write("[loci]=%d\n\n[populations]=%d\n\n" % (len(output_lines_bypop[output_lines_bypop.keys()[0]]), len(output_lines_bypop)))
    pop_count=0
    for each_pop in output_lines_bypop:
        pop_count+=1
        fout.write("\n[pop]=%d\n" % (pop_count))
        print pop_count, each_pop
        for line in output_lines_bypop[each_pop]:
            fout.write("%s\n" % (line))

def check_for_min_data(merged_genos_dict, min_data=5):
    for pop_info in merged_genos_dict.values():
        if len(pop_info) - pop_info.count(".") < min_data * 2: return False
    return True

def write_alleles(outfile_pointer, locus_name, locus_genos):
    return True

#---------------------------------------------------------------------------------------------------------------
###----------> This function will join a list of values using a specified delimiter, the default is tab
###----------> Works even if some of the values in the list are tuples and integers
###------------------------------------------------------------------------------------------------------------

def recursive_join(list, delimiter="\t"):				# Function definition: def NAME(list of parameters):
	ready_to_join = []									#							STATEMENTS
	for index, value in enumerate(list):
		if type(value) == type(()) or type(value) == type([]):
			ready_to_join.append(recursive_join(value))
		elif type(value) == type(1) or type(value) == type(1.0): 
			ready_to_join.append(str(value)) 
		else: 
			ready_to_join.append(value)					# end of the function calls itself, hence recursive
			
	joined=delimiter.join(ready_to_join)
	return joined


###------------> Actual script starts here:
#For pipelines class
indivs2pops_dict={'FR32_ATCACG_before': 'Before', 'FR33_TTAGGC_before': 'Before', 'FR34_ACTTGA_before': 'Before', 'FR35_GATCAG_before': 'Before', 'FR36_TAGCTT_before': 'Before', 'FR37_GGCTAC_before': 'Before', 'FR77_ATCACG_after': 'After', 'FR78_TTAGGC_after': 'After', 'FR79_ACTTGA_after': 'After', 'FR80_GATCAG_after': 'After', 'FR81_TAGCTT_after': 'After', 'FR82_GGCTAC_after': 'After'}

make_bayescan_input(sys.argv[1], int(sys.argv[2]), indivs2pops_dict, 1)
