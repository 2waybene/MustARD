dat <- read.table ("/Users/li11/myGit/MustARD/temp/CCLE_RNAseq_genes_rpkm_20180929.txt", header = TRUE, sep="\t")
colnames(dat) %in% "LUNG"

library(stringr)
which(str_detect(colnames(dat), "OESOPHAGUS"))

ESCC.dat <- dat[,c(1,2, which(str_detect(colnames(dat), "OESOPHAGUS")))]

dim(ESCC.dat)
# total 27 ESCC cell lines

unique(ESCC.dat$Description)

geneSet <- read.csv("/Users/li11/Downloads/mmNRF2ChIPseqES_1923_genes_to_hs.csv")

geneSet$HGNC
