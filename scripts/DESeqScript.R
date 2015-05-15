### This is an R script that implements the program DESeq for gene expression analysis.
### Much more information on the program and specific function (particularly for checking quality) 
### can be found here: http://bioconductor.org/packages/2.8/bioc/vignettes/DESeq/inst/doc/DESeq.pdf
### The elements of this script were written by Melissa Pespeni and Daniel Barshis.

#setwd("YOURWORKINGDIRECTORYPATH")  # The drag and drop from finder works in R, too.
library(DESeq)

#useful functions
#head() - prints out the top 6 lines
#dim() - prints the dimensions of a variable
#nrow() - returns the number of rows in a vector or matrix
# ?[functionName] - opens documentation describing the function

#read in your data to make counts table
countsTable <- read.delim('combinedcounts.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(countsTable)

#define the conditions of each individual (e.g. "Cold" and "Hot")  There should be the same number of conditions described as there are samples in your data file, and in the same order
conds <- c("Cold","Cold","Cold","Cold","Cold","Cold","Hot","Hot","Hot","Hot","Hot","Hot") 

#make count data sets
cds <- newCountDataSet(countsTable, conds)
cds<-estimateSizeFactors(cds)
sizeFactors(cds)			# factor such that values in a column can be brought to a common scale by dividing by the corresponding size factor
cds<-estimateVarianceFunctions(cds)

#negative binomial testing (e.g. compare "Cold" vs. "Hot" for differences in gene expression)
res <- nbinomTest(cds, "Cold", "Hot")
head(res)

#count the number of significantly differentially expressed genes
nrow(res[res$padj<0.05 & !is.na(res$padj),]) # at the .05 level
nrow(res[res$padj<0.01 & !is.na(res$padj),]) # at the .01 level
											 
#filter for contigs with average(baseMean) >5
res5<-res[res$baseMean>5, ]
dim(res)
dim(res5)  # number of genes that have >5 counts

#p-value readjustment Benjamini and Hochberg after >5 filtering
res5$padj <- p.adjust(res5$pval, method="BH")

nrow(res5[res5$padj<0.01 & !is.na(res5$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes

#### Now filter for average counts and variance

# Make a counts table that is scaled by the size factors
temp = t(sizeFactors(cds))
sizematrix<-matrix(data=temp, nrow=nrow(countsTable), ncol=ncol(temp), byrow=TRUE)
scaledcounts = countsTable/sizematrix
head(scaledcounts)

# Calculate the standard deviation across all individuals for each gene
sdcol<-apply(scaledcounts,1,sd)
head(sdcol)
sdcolT<-data.frame(sdcol)
head(sdcolT)

# Combine the results data frame with the calculated standard deviation
ressd<-cbind(res,sdcolT)
head(ressd)

ressdFilter<-ressd[ressd$baseMean>5 & ressd$sdcol<ressd$baseMean, ]
dim(ressd)
dim(ressdFilter)  # number of genes that have >5 counts and where the mean is greater than the std deviation


#p-value readjustment Benjamini and Hochberg after >5 filtering
ressdFilter$padj <- p.adjust(ressdFilter$pval, method="BH")
head(ressdFilter)

nrow(ressdFilter[ressdFilter$padj<0.01 & !is.na(ressdFilter$padj),])  # Num significantly differentially expressed genes in the filtered data frame


#### Output files that can be generated:

id<-data.frame(rownames(scaledcounts))  # making a data frame of the rownames from the scaledcounts data frame
colnames(id)<-c("id")  # renaming the column header to "id"
scaledcountsID<-cbind(id,scaledcounts)  # combining the columns (column binding)

CountsAndResTable<-merge(scaledcountsID,res, sort=FALSE)  # merge the two tables; uses the common 'id' column to link the tables
head(CountsAndResTable)
dim(CountsAndResTable)

write.table(CountsAndResTable, file="ScaledcountsAndResultsOutput.txt", sep="\t", row.names=F)


CountsAndRessdTable<-merge(scaledcountsID,ressdFilter, sort=FALSE)  # merge the two tables; uses the common 'id' column to link the tables
head(CountsAndRessdTable)
dim(CountsAndRessdTable)

write.table(CountsAndRessdTable, file="ScaledcountsAndFilteredResultsOutput.txt", sep="\t", row.names=F)


write.table(ressdFilter, file="FilteredResultsOutput.txt", sep="\t", row.names=F)

CountsAndResTable<-merge(scaledcountsID,ressdFilter)
head(CountsAndResTable)
dim(CountsAndResTable)

write.table(CountsAndResTable, file="CountsAndResultsOutput.txt", sep="\t", row.names=F)

Sig<-data.frame(CountsAndResTable[ressdFilter$padj<0.01 & !is.na(ressdFilter$padj),])
head(Sig)
dim(Sig)
write.table(Sig, file="Sig_CountsAndResultsOutput.txt", sep="\t", row.names=F)

library(gplots)
heatmap.2(data.matrix(scaledcounts))
