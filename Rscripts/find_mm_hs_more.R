source ("/Users/li11/myGit/MustARD/Rscripts/utils.R")

##============================
##  
mmEsGeneSets <- read.csv ("/Users/li11/myGit/MustARD/doc/mmNRF2ChIPseqES.csv")
dim(mmEsGeneSets)


mmEsGeneSets <- read.csv ("/Users/li11/myGit/MustARD/doc/mmNRF2ChIPseqES_1923_genesOnly.csv")
dim(mmEsGeneSets)
head(mmEsGeneSets)
system.time(
  hsGenes_mmEsGeneSets <- convertMouseGeneList(mmEsGeneSets$Gene_Name)
)


converted <- as.data.frame(genesV2)
head(converted)
dim(converted)
colnames(converted ) <- c("MGI", "HGNC")

updated.hs.genes <- merge(mmEsGeneSets, converted, by.x ="Gene_Name", by.y="MGI", all = TRUE)
dim(mmEsGeneSets)
dim(updated.hs.genes)
      
length(which(mmEsGeneSets$Gene_Name %in% updated.hs.genes$Gene_Name))
length(which(updated.hs.genes$Gene_Name %in% mmEsGeneSets$Gene_Name))


length(updated.hs.genes$Gene_Name)
length(unique(updated.hs.genes$Gene_Name))
length(unique(updated.hs.genes$HGNC))
length(updated.hs.genes$HGNC)

df = as.character(updated.hs.genes$Gene_Name)
df = as.data.frame(df)
duplicated.rows <- which(duplicated(df) | duplicated(df[nrow(df):1, ])[nrow(df):1])
updated.hs.genes[duplicated.rows ,]

write.csv(updated.hs.genes, file = "/Users/li11/myGit/MustARD/doc/mmNRF2ChIPseqES_1923_genes_to_hs.csv", row.names = FALSE)
