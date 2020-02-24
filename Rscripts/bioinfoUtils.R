##====================================================================
##  File  : bioinfoUtils.R
##  Author: Jianying Li
##  Comment: a collection of R utility scripts for bioinformatic
##           related data processing and manipulation
##====================================================================

##  dependent libraries:

require("biomaRt")

##  Help doc: 

##  Function #1:  convert.factors.to.strings.in.dataframe
##              Converting factors into strings in a standard dataframe
##  input: a dataframe of factors
##  output: a dataframe of characters


##  Function #2: mm2hs_entrezSymbol 
##           convert mouse gene symbols to human gene symbols
##  input: a vector of mouse gene symbols 
##  output: a vecor of human gene symbols


convert.factors.to.strings.in.dataframe <- function(dataframe)
{
  class.data  <- sapply(dataframe, class)
  factor.vars <- class.data[class.data == "factor"]
  for (colname in names(factor.vars))
  {
    dataframe[,colname] <- as.character(dataframe[,colname])
  }
  return (dataframe)
}


# Basic function to convert mouse to human gene names
mm2hs_entrezSymbol <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}
