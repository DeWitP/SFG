#!/usr/bin/perl -w

# Parsing BLAST reports with BioPerl's Bio::SearchIO module
# WI Biocomputing course - Bioinformatics for Biologists - October 2003

# See help at http://www.bioperl.org/wiki/HOWTO:SearchIO for all data that can be extracted
use lib '/Users/TomOliver/Pipeline/BioPerl-1.6.0/blib/lib';
use Bio::SearchIO;

# Prompt the user for the file name if it's not an argument
# NOTE: BLAST file must be in text (not html) format
if (! $ARGV[0])
{
   print "What is the BLAST file to parse? ";

   # Get input and remove the newline character at the end
   chomp ($inFile = <STDIN>);
}
else
{
   $inFile = $ARGV[0];
}

$report = new Bio::SearchIO(
         -file=>"$inFile",
              -format => "blastxml"); 

print "QueryAcc\tQuery_Length\tHitDescription\tHitName\tHItLength\tHitBits\tHSP_rank\t\%ID\teValue\tQuery_Start\tQuery_End\tHit_start\tHit_end\tHSP_length\n";

# Go through BLAST reports one by one              
while($result = $report->next_result) 
{
   # Go through each each matching sequence
   while($hit = $result->next_hit) 
   { 
      # Go through each each HSP for this sequence
        while ($hsp = $hit->next_hsp)
         { 
            # Print some tab-delimited data about this HSP
            # Check out http://www.bioperl.org/wiki/HOWTO:SearchIO for all the info you can add to the output
            
            print $result->query_description, "\t";
            print $result->query_length, "\t"; # Query length
#			print $result->query_accession, "\t";
            print $hit->description, "\t";  #commented out because returned nothing with nemvedatabase
            print $hit->name, "\t";    #added to see if name of hits would come up in parsed output it worked!! but need to edit column headers appropriately
 			print $hit->length, "\t"; 	# Length of the Hit sequence 
 #          print $hit->significance, "\t"; #redundant, same as evalue
 			print $hit->bits, "\t";
 			print $hsp->rank, "\t";
            print $hsp->percent_identity, "\t";
            print $hsp->evalue, "\t";
            print $hsp->start('query'), "\t"; 	 # query start position from alignment
		 	print $hsp->end('query'), "\t"; # query end position from alignment
 			print $hsp->start('hit'), "\t"; # hit start position from alignment
 			print $hsp->end('hit'), "\t"; 	# hit end position from alignment start and end of query as array 
			print $hsp->hsp_length, "\n";
      } 
   } 
}


		