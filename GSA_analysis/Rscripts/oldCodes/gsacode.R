library(shiny)
library(shinyFiles)
library(samr)
#setwd("C:/Users/Huang/Desktop/microarray")
setwd("X:/project2016/GSA_by_R/")
getwd()


x1 <- read.csv("normalization Pax9 siRNA vs. Control.csv")
str(x1)
str(matrix(rnorm(1000*9),ncol=9))
x <- x1[3:10]
y <- c(rep(1, 4), rep(2, 4))
y
samfit <- SAM(x, y, resp.type = "Two class unpaired")
plot(samfit)
str(x)
z <- x1
z$X <- as.character(z$X)
str(z)
z$GeneSymbol <- as.character(z$GeneSymbol)
str(z)
data <- list(x = x, y = y, geneid = z$X, genenames = z$GeneSymbol, logged2 = TRUE)

library(GSA)
set.seed(1234567)
#setwd("C:/Users/Huang/Desktop/microarray")
setwd("X:/project2016/GSA_by_R/")
#c1 <- GSA.read.gmt("Knowledge based gene sets (20140627).gmt")
c1 <- GSA.read.gmt("Knowledge based gene sets (20161005).gmt")
str(c1)
str(c1$genesets)
str(c1$geneset.names)
c1$geneset.names

c1$geneset.names[2] == c1$geneset.names[17]

gsafit <- GSA(x, y, c1$genesets, z$GeneSymbol, resp.type = "Two class unpaired", nperms = 1000)
str(gsafit)
plot(gsafit$GSA.scores)
gsafit$fdr.lo
gsafit$fdr.hi
gsafit$gene.scores
GSA.plot(gsafit)
GSA.listsets(gsafit, c1$geneset.names, FDRcut=.5)
GSA.correlate(c1$genesets, z$GeneSymbol)
gsa2 <- GSA.func(x, y, c1$genesets, z$GeneSymbol, resp.type = "Two class unpaired")
str(gsa2)

#for loop from 1:length(gsals$negative[,1])
GSA.genescores(5, c1$genesets, gsafit, z$GeneSymbol)
GSA.plot(gsafit, fac = 1, FDRcut = 1)

gsals <- GSA.listsets(gsafit, c1$geneset.names, FDRcut=.5)
str(gsals)
str(gsals$negative)
gsals$negative[,1]
length(gsals$negative[,1])
as.numeric(gsals$negative[,1][1])
GSA.genescores(as.numeric(gsals$negative[,1][1]), c1$genesets, gsafit, z$GeneSymbol)

sink("GSA-TF-negative.txt")
for (i in 1:length(gsals$negative[,1])){
  print(gsals$negative[,2][i])
  print(GSA.genescores(as.numeric(gsals$negative[,1][i]),
                 c1$genesets, gsafit, z$GeneSymbol))
  cat("\n")
  
}
sink()

sink("GSA-TF-positive.txt")
for (i in 1:length(gsals$positive[,1])){
  print(gsals$positive[,2][i])
  print(GSA.genescores(as.numeric(gsals$positive[,1][i]),
                       c1$genesets, gsafit, z$GeneSymbol))
  cat("\n")
  
}
sink()

sink("GSA-TF-output.txt")
print(GSA.listsets(gsafit, c1$geneset.names, FDRcut=.5))
sink()

sink("negativenames.txt")

for (i in as.numeric(gsals$negative[,1])){
  print(c1$geneset.names[i])
  
}
sink()

sink("positivenames.txt")

for (i in as.numeric(gsals$positive[,1])){
  print(c1$geneset.names[i])
  
}
sink()

