#!/bin/bash

# Variant recalibrator. 
# The first command train a model for SNPs, the second for InDels.
# Due to the small size of out dataset, we might get this error message:

### ERROR MESSAGE: NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe. 
### Please consider raising the number of variants used to train the negative model (via --percentBadVariants 0.05, for example) 
### or lowering the maximum number of Gaussians to use in the model (via --maxGaussians 4, for example).

#If you do get this message, try changing the settings for -percentBad, -minNumBad and --maxGaussians in the following two commands.


java -Xmx4g -jar ../scripts/GenomeAnalysisTK_2_5.jar \
   -T VariantRecalibrator \
   -R ./assembly/Trinity.fasta \
   -input raw_snps_indels_Q3_annotated.vcf \
   -percentBad 0.05 -minNumBad 50 --maxGaussians 2 \
   -resource:concordantSet,known=true,training=true,truth=true,prior=10.0 raw_snps_indels_Q20.vcf \
   -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP \
   -mode SNP \
   -recalFile VQSR_SNP.recal \
   -tranchesFile VQSR_SNP.tranches


java -Xmx4g -jar ../scripts/GenomeAnalysisTK_2_5.jar \
   -T VariantRecalibrator \
   -R ./assembly/Trinity.fasta \
   -input raw_snps_indels_Q3_annotated.vcf \
   --maxGaussians 1 -percentBad 0.05 -minNumBad 50 \
   -resource:concordantSet,known=true,training=true,truth=true,prior=10.0 raw_snps_indels_Q20.vcf \
   -an DP -an FS -an ReadPosRankSum -an MQRankSum \
   -mode INDEL \
   -recalFile VQSR_INDEL.recal \
   -tranchesFile VQSR_INDEL.tranches



#Applying the recalibration
   
java -Xmx3g -jar ../scripts/GenomeAnalysisTK_2_5.jar \
   -T ApplyRecalibration \
   -R ./assembly/Trinity.fasta \
   -input raw_snps_indels_Q3_annotated.vcf \
   -tranchesFile VQSR_SNP.tranches \
   -recalFile VQSR_SNP.recal \
   -o recalibratedSNPs.rawIndels.vcf \
   --ts_filter_level 99.9 \
   -mode SNP
   
java -Xmx3g -jar ../scripts/GenomeAnalysisTK_2_5.jar \
   -T ApplyRecalibration \
   -R ./assembly/Trinity.fasta \
   -input recalibratedSNPs.rawIndels.vcf \
   -tranchesFile VQSR_INDEL.tranches \
   -recalFile VQSR_INDEL.recal \
   -o analysis_ready.vcf \
   --ts_filter_level 99.9 \
   -mode INDEL


#Finally, save all the SNPS that have passed the VQSR filter into a new vcf file:

grep "PASS\|^#" analysis_ready.vcf > VQSR_PASS_SNPS.vcf
