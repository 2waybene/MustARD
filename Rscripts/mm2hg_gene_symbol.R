source ("x:/project2020/MustARD/Rscripts/bioinfoUtils.R")



mmEsGeneSets <- read.csv ("x:/project2020/MustARD/geneSets/mm_ES_NRF2-ChiP_geneset.csv")
system.time(
  hsGenes_mmEsGeneSets <- mm2hs_entrezSymbol(convert.factors.to.strings.in.dataframe(mmEsGeneSets)$MouseGene)
)

dim(mmEsGeneSets)
head(mmEsGeneSets)
length(hsGenes_mmEsGeneSets)
head(hsGenes_mmEsGeneSets)

library(readxl)
newGenelist <- read_excel("x:/project2020/MustARD/doc/mmNRF2ChIPseqES_1923_genes.xlsx") # All gene signature (high/low) table


mmList <- newGenelist$`Gene Name`
system.time(
  hsGenes_mmEsGeneSets <- mm2hs_entrezSymbol(mmList)
)

length(mmList)
#[1] 1923
length(hsGenes_mmEsGeneSets)
#[1] 1543

newGenelist <- read.csv("x:/project2020/MustARD/doc/mmNRF2ChIPseqES_1923_genesOnly.csv") # All gene signature (high/low) table
head(newGenelist)
