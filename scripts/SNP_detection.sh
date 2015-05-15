#!/bin/bash

# Calling only very high-quality SNPs with the Haplotype caller (recommended for all diploid samples, if analyzing population samples, use Unified Genotyper instead):
# This will likely find only the most obvious variant sites. We will use this dataset as "true" sites for the Variant Score Quality Recalibrator (VQSR) in the next step.

	
java -Xmx16g -jar ../scripts/GenomeAnalysisTK_2_5.jar \
     -T HaplotypeCaller \
     -R ./assembly/Trinity.fasta \
     -I merged_realigned.bam \
     -stand_call_conf 20.0 \
     -stand_emit_conf 20.0 \
     -o raw_snps_indels_Q20.vcf


#Selecting an appropriate quality score threshold (from the Broad Institute's Wiki site):
#A common question is the confidence score threshold to use for variant detection. We recommend:
#Deep (> 10x coverage per sample) data 
#    we recommend a minimum confidence score threshold of Q30 with an emission threshold of Q10. These Q10-Q30 calls will be emitted filtered out as LowQual. 
#Shallow (< 10x coverage per sample) data 
#    because variants have by necessity lower quality with shallower coverage, we recommend a min. confidence score of Q4 and an emission threshold of Q3. 


#Now, we can call variants with a threshold appropriate for our sequence coverage

java -Xmx16g -jar ../scripts/GenomeAnalysisTK_2_5.jar \
     -T HaplotypeCaller \
     -R ./assembly/Trinity.fasta \
     -I merged_realigned.bam \
     -stand_call_conf 3.0 \
     -stand_emit_conf 3.0 \
     -o raw_snps_indels_Q3.vcf

#And finally, we annotate the list of variants for the VariantRecalibration to be run properly

java -Xmx16g -jar ../scripts/GenomeAnalysisTK_2_5.jar \
  -T VariantAnnotator \
  -R ./assembly/Trinity.fasta \
  -I merged_realigned.bam \
  -o raw_snps_indels_Q3_annotated.vcf \
  --variant raw_snps_indels_Q3.vcf \
  -L raw_snps_indels_Q3.vcf