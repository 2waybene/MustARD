library(GSA)
setwd("x:/project2016/GSA_by_R/")

##==================================================================================
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

##==================================================================================
geneset<- GSA.read.gmt( "cleaningData/Knowledge_based_gene_sets_20161005.gmt")
geneset<- GSA.read.gmt( "cleaningData/Knowledge_based_gene_sets_20140627.gmt")

c1 <- GSA.read.gmt( "cleaningData/Knowledge_based_gene_sets_20140627.gmt")
c1 <- GSA.read.gmt("Knowledge based gene sets (20161005).gmt")

str(geneset)
str(geneset$genesets)
geneset$genesets[[17]] <- toupper (geneset$genesets[[17]])



listToUpper <-  function (listIN)
{
  list2return <- c(list())
  for (i in 1:length(listIN)){
    list2return[[i]] <- toupper(listIN[[i]])
  }
  return (list2return)
}


##  Initial gene sets
set.seed(1234567)
gsafit <- GSA(x, y, geneset$genesets, z$GeneSymbol, resp.type = "Two class unpaired", nperms = 1000)
print(GSA.listsets(gsafit, geneset$geneset.names , FDRcut=.5))


##  Gene sets are capitalized
set.seed(1234567)
genesets.in.upperCase <-listToUpper (geneset$genesets)
gsafit <- GSA(x, y,genesets.in.upperCase, toupper(z$GeneSymbol), resp.type = "Two class unpaired", nperms = 1000)
print(GSA.listsets(gsafit, geneset$geneset.names , FDRcut=.5))


cleaned.geneset <- clean.gmt.data (geneset)
cleaned.geneset$geneset.names[[17]] <- "V2Sox2 target genes in esophagus"
genesets.in.upperCase <-listToUpper (cleaned.geneset$genesets)
set.seed(1234567)
gsafit <- GSA(x, y,genesets.in.upperCase, toupper(z$GeneSymbol), resp.type = "Two class unpaired", nperms = 1000)

print(GSA.listsets(gsafit, cleaned.geneset$geneset.names , FDRcut=.5))


#geneset $genesets[18]
print ("Gene sets in database file")
geneset$geneset.names


genesets.new <- clean.gmt.data(geneset)

genesets.new$geneset.names[[17]] <- "V2Sox2 target genes in esophagus"
str(genesets.new$genesets[[17]])


##======================================
##  Try to compare these two files
##======================================

geneset.1 <- GSA.read.gmt( "cleaningData/genset_1_from_20161005.txt")
geneset.2 <- GSA.read.gmt( "cleaningData/genset_1_mod_3.txt")

str(geneset.1)
geneset.1$geneset.names[[17]]

str(geneset.2)

geneset.2$geneset.names[[17]]



gsafit <- GSA(x, y, geneset.1$genesets[c(5,12,13,14,17)], z$GeneSymbol, resp.type = "Two class unpaired", nperms = 1000)
print(GSA.listsets(gsafit, geneset.1$geneset.names[c(5,12,13,14,17)] , FDRcut=.5))




geneset.1$geneset.names[c(5,1)]

c1$genesets[[18]]
c1.cleaned$genesets[[18]]


c1 <- GSA.read.gmt("cleaningData/Knowledge_based_gene_sets_20161005.gmt")
c2 <- GSA.read.gmt("cleaningData/temp2.txt")

str(c1$genesets)
str(c2$genesets)

c1.cleaned <-  clean.gmt.data(c1)
c2.cleaned <- clean.gmt.data(c2)

str(c1.cleaned$genesets)
str(c2.cleaned$genesets)


c1.cleaned$genesets[[17]][-1] %in% toupper (c2.cleaned$genesets[[17]])

c1.mod <- c(list())
c2.mod <- c(list())
for (i in 1:18)
{
  cat(i); cat("\n")
  cat(c1.cleaned$genesets[[i]] %in% c2.cleaned$genesets[[i]])
  cat("\n")
  c1.mod[[i]] = c1.cleaned$genesets[[i]]
  c2.mod[[i]] = c2.cleaned$genesets[[i]]
}

setNames <- c1.cleaned$geneset.names


gsafit <- GSA(x, y, c1.cleaned$genesets, z$GeneSymbol, resp.type = "Two class unpaired", nperms = 1000)
print(GSA.listsets(gsafit, c1$genesets. , FDRcut=.5))

gsafit <- GSA(x, y, c1.mod, z$GeneSymbol, resp.type = "Two class unpaired", nperms = 1000)

gsafit <- GSA(x, y, c2.cleaned$genesets, z$GeneSymbol, resp.type = "Two class unpaired", nperms = 1000)

gsafit <- GSA(x, y, c2.mod, z$GeneSymbol, resp.type = "Two class unpaired", nperms = 1000)
print(GSA.listsets(gsafit, setNames , FDRcut=.5))


i = 17

c1.cleaned$genesets[[i]]
c2.cleaned$genesets[[i]]


c2.mod[[i]] <- c2.cleaned$genesets[[i]]
c2.mod[[i]] <- NULL
c1.mod[[i]] <- NULL


for (i in 1:18)
{
  cat(i); cat("\n")
  cat(c1.mod[[i]] %in% c2.mod[[i]])
  cat("\n")
}

length(c1.mod)

c2.mod[[18]] <- toupper (c2.cleaned$genesets[[17]])
loweredSet <- tolower(c2.mod[[18]])
c2.mod[[18]] <- loweredSet

c1.mod[[18]] <- c1.cleaned$genesets[[17]]


