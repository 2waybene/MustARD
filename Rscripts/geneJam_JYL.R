load("x:/project2020/Rscripts/mouseAffyProbe.rda")

##===========================================
##  JYL
##===========================================

my.affyProb2Symbols <- function (aListOfProbes, db = mouse4302.db)
{
  require(magrittr)
  require(mouse4302.db)
  probeAnnoted <- AnnotationDbi::select(
    x       = db,
    keys    = aListOfProbes,
    columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
    keytype = "PROBEID"
  )
  return (probeAnnoted )
}

##===========================================

str(cleanID)
#chr [1:45101] "1415670_at" "1415671_at" "1415672_at" "1415673_at" "1415674_a_at" "1415675_at" "1415676_a_at" "1415677_at" "1415678_at" "1415679_at" ...

system.time(
  getSymbol.me <-my.affyProb2Symbols(cleanID )
)

##  'select()' returned 1:many mapping between keys and columns
##  user  system elapsed 
##  0.43    0.06    0.50 

dim(getSymbol.me)
#[1] 48348     4
head(getSymbol.me)

system.time(
  getSymbol.james <-genejam::freshenGenes(cleanID,
                                          ann_lib=c("mouse4302.db", "org.Mm.eg.db"),
                                          include_source=TRUE,
                                          try_list=c("ENTREZID","REFSEQ2EG","SYMBOL2EG","ACCNUM2EG","ALIAS2EG"),
                                          final=c("SYMBOL"))
)

# user  system elapsed 
# 9.75    0.00    9.75 
dim(getSymbol.james)
#[1] 45101     5
head(getSymbol.james)