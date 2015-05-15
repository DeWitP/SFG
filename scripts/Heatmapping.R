library(gplots)

# type "?heatmap.2" to get the R documentation for heatmap.2 that describes the parameters and provides examples (also search the interwebs)

mydata<-read.delim(file.choose(), row.names=1) # select your .txt file of normalized fold difference or counts data

head(mydata) # check out your data

pairs.breaks <- seq(-3.0, 3.0, by=0.1)
length(pairs.breaks)  # take this value minus one as the input for n in the line below

mycol <- colorpanel(n=60, low="blue",mid="black",high="yellow") #n needs to be 1 less than length pairs.breaks

heatmap.2(data.matrix(mydata), Rowv=T, Colv=F, dendrogram = c("row"), scale = "none", keysize=1, breaks = pairs.breaks, col=mycol, trace = "none", symkey = F, density.info = "none", labCol=c(""), colsep=c(24), sepcolor=c("white"), sepwidth=c(.1,.1))