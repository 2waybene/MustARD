##=================================================
##  Script Name: gsacode_v2.R
##  Author: Jianying Li
##  Date: 10/27/2016
##  History: initally coded by Yicheng Li
##  Comment: Modified by JYL to clean geneset
##           and some cleaning work
##=================================================

##  Load required packages and set up working directory

library(shiny)
library(shinyFiles)
library(samr)

##  Set up working directory, i.e. 
setwd("X:/project2016/GSA_by_R/")


##  Working on Microarry data 

x1 <- read.csv("normalization Pax9 siRNA vs. Control.csv")

##  The following codes need to be stardardized
# str(x1)   ## Not necessary -- JYL
# str(matrix(rnorm(1000*9),ncol=9))   ## Not necessary -- JYL
x <- x1[3:10]
y <- c(rep(1, 4), rep(2, 4))
#y

##  Suggested modification -- JYL
# x <- x1[-c(1,2)]


##  Fitting SAM model NOT seems necessary though-- JYL

samfit <- SAM(x, y, resp.type = "Two class unpaired")
plot(samfit)
#str(x)
z <- x1
z$X <- as.character(z$X)
#str(z)
z$GeneSymbol <- as.character(z$GeneSymbol)
#str(z)
data <- list(x = x, y = y, geneid = z$X, genenames = z$GeneSymbol, logged2 = TRUE)



##  GSA analysis starts here-- JYL  
library(GSA)
set.seed(1234567)   ## EXTREMELY important!!  -- JYL 
setwd("C:/Users/Huang/Desktop/microarray")

##==========================================================
##  Function for clean and capitalize the gene set symbols
##=========================================================
listToUpper <-  function (listIN)
{
  list2return <- c(list())
  for (i in 1:length(listIN)){
    list2return[[i]] <- toupper(listIN[[i]])
  }
  return (list2return)
}

clean.gmt.data <- function (db)
{
  temp.geneset <- c(list())
  for (i in 1:length(db$genesets))
  {
    temp.geneset[[i]] <- db$genesets[[i]]
    if (length(which(temp.geneset[[i]]=="")) > 0 )
    {
      temp.geneset[[i]]  <- temp.geneset[[i]][-which(temp.geneset[[i]]=="")]
    }
  }
  db$genesets <- temp.geneset 
  return (db)
}
##=========================================================
##  Function end    -- JYL
##=========================================================


c1 <- GSA.read.gmt("Knowledge based gene sets (20140627).gmt")
c1 <- GSA.read.gmt("Knowledge based gene sets (20161005).gmt")

#str(c1)
#str(c1$genesets)
#str(c1$geneset.names)
#c1$geneset.names


cleaned.geneset       <- clean.gmt.data (c1)
genesets.in.upperCase <-listToUpper(cleaned.geneset$genesets)
set.seed(1234567) ## EXTREMELY important!!  -- JYL 
gsafit <- GSA(x, y,genesets.in.upperCase, toupper(z$GeneSymbol), resp.type = "Two class unpaired", nperms = 1000)

##  gsafit <- GSA(x, y, c1$genesets, z$GeneSymbol, resp.type = "Two class unpaired", nperms = 1000)  ## YL

##  Examine GSA results  ## Non-essential -- JYL
str(gsafit)
plot(gsafit$GSA.scores)
gsafit$fdr.lo
gsafit$fdr.hi
gsafit$gene.scores
GSA.plot(gsafit)
GSA.listsets(gsafit, c1$geneset.names, FDRcut=.5)
GSA.correlate(c1$genesets, z$GeneSymbol)


##  It is suggested the whenever GSA model is run, 
##  set seed to ensure reproducibility
set.seed(1234567) ## EXTREMELY important!!  -- JYL 
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

