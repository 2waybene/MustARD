Knowledge_based_gene_sets_20140627.gmt

c1 <- GSA.read.gmt("Knowledge_based_gene_sets_20140627.gmt")
c1 <- GSA.read.gmt("Knowledge_based_gene_sets_20161005.gmt")
c1 <- GSA.read.gmt("just2test.gmt")


#just to test




##  Remove the last TWO "questionable genesets" from 20161005.gmt"

str(c1)
c1$genesets[[18]] <- NULL
c1$genesets[[17]] <- NULL
c1$geneset.names  <- c1$geneset.names[-c(17,18)]  
c1$geneset.descriptions <- c1$geneset.descriptions[-c(17,18)]



##  Remove the last "questionable geneset" from 20140627.gmt

str(c1)
c1$genesets[[17]] <- NULL
c1$geneset.names  <- c1$geneset.names[-17]  
c1$geneset.descriptions <- c1$geneset.descriptions[-17]

cleaned.geneset       <- c1
genesets.in.upperCase <-cleaned.geneset$genesets

cleaned.geneset       <- clean.gmt.data (c1)
genesets.in.upperCase <-listToUpper(cleaned.geneset$genesets)



set.seed(1234567) ## EXTREMELY important!!  -- JYL
gsafit <- GSA(x, y,genesets.in.upperCase, toupper(z$GeneSymbol), resp.type = "Two class unpaired", nperms = 1000)
##  Examine GSA results  ## Non-essential -- JYL
str(gsafit)
plot(gsafit$GSA.scores)
gsafit$fdr.lo
gsafit$fdr.hi
gsafit$gene.scores
GSA.plot(gsafit)
GSA.listsets(gsafit, c1$geneset.names, FDRcut=.5)
GSA.correlate(c1$genesets, z$GeneSymbol)
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
print(GSA.listsets(gsafit, c1$geneset.names, FDRcut=.5))

      