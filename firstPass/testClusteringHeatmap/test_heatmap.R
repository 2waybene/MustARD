##==========================
##  simulat normal data
##===========================
library(gplots) 
y <- matrix(rnorm(500), 100, 5, dimnames=list(paste("g", 1:100, sep=""), paste("t", 1:5, sep=""))) 
heatmap.2(y) # Shortcut to final result


## Row- and column-wise clustering
hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete")
## Tree cutting
mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)]
## Plot heatmap
mycol <- colorpanel(40, "darkblue", "yellow", "white") # or try redgreen(75)
mycol <- redgreen(75)
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=mycol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc) 


##==========================
##  ESCC data on 229 genes
##===========================

setwd("/Users/jyli/myGit/MustARD/firstPass/testClusteringHeatmap/")
setwd("/Users/li11/myGit/MustARD/firstPass/testClusteringHeatmap/")
dt <- read.table ("RNASeq_ESCC_on_229_genes.txt",header = TRUE)

dt[dt == 0] = 0.001
which(dt == 0 )


head(dt[,-1])
str(dt[,-1])

y = as.matrix(dt[,-1])
y = log10(as.matrix(dt[,-1]))

heatmap(y)
heatmap.2(log10(y))

mycol <- redgreen(75)
heatmap.2(log10(y),  col=mycol , density.info="none", trace="none")


