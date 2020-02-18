##======================================================================================================
##  mm2hg_gene_symbol.R
##  credit: https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
##=====================================================================================================
source ("x:/project2020/MustARD/Rscripts/utils.R")
source("x:/R-project/customPackages/plotTools.R")
source("x:/R-project/customPackages/dataManipTools.R")

##============================
##  
mmEsGeneSets <- read.csv ("x:/project2020/MustARD/doc/mmNRF2ChIPseqES.csv")
dim(mmEsGeneSets)
#[1] 1940   26
mmEsGeneSets.wt          <- mmEsGeneSets$Gene.Name[(mmEsGeneSets$X1_WT_NRF2..1.Present == 1)]
mmEsGeneSets.Nrf2dMinus  <- mmEsGeneSets$Gene.Name[(mmEsGeneSets$X2_Nrf2_dbminus_NRF2..1.Present == 1)]
mmEsGeneSets.Keep1dMinus <- mmEsGeneSets$Gene.Name[(mmEsGeneSets$X3_Keap1_dbminus_NRF2..1.Present == 1)]

length(mmEsGeneSets.wt)
length(mmEsGeneSets.Nrf2dMinus)
length(mmEsGeneSets.Keep1dMinus)

mmEsGeneSets.Nrf2dMinus <- mmEsGeneSets$Gene.Name

##======================================================================================================
##  Testing procedures here
##==============================

#mmEsGeneSets <- read.csv ("x:/project2020/MustARD/geneSets/mm_ES_NRF2-ChiP_geneset.csv")
#system.time(
#  hsGenes_mmEsGeneSets <- convertMouseGeneList(convert.factors.to.strings.in.dataframe(mmEsGeneSets)$MouseGene)
#)
#write.csv (hsGenes_mmEsGeneSets, file ="x:/project2020/MustARD/geneSets/hsGenes_mm_ES_NRF2-ChiP_geneset.csv", row.names = FALSE)

hsGenes_mmEsGeneSets <- read.csv("x:/project2020/MustARD/geneSets/hsGenes_mm_ES_NRF2-ChiP_geneset.csv")$x

length(hsGenes_mmEsGeneSets)
which(is.na(hsGenes_mmEsGeneSets))
#integer(0)
dim(mmEsGeneSets)

##===============================
##  compare to CRISPR paper
##===============================
load("X:/project2020/MustARD/learningFromBigWigs/Behan_nature_CRISPR-Cas9/Behan_CRISPR_ESCC_hit_genes.rda")
length(ESCC_hit_gene_sets$FitAll19)
length(ESCC_hit_gene_sets$Fit18More)
length(ESCC_hit_gene_sets$Fit15More)
length(ESCC_hit_gene_sets$Fit10More)
length(ESCC_hit_gene_sets$FitAny)

length(hsGenes_mmEsGeneSets)
#[1] 1559
length(intersect(as.character(ESCC_hit_gene_sets$FitAll19), hsGenes_mmEsGeneSets))
#[1] 27
length(intersect(as.character(ESCC_hit_gene_sets$Fit18More), hsGenes_mmEsGeneSets))
#[1] 43
length(intersect(as.character(ESCC_hit_gene_sets$Fit15More), hsGenes_mmEsGeneSets))
#[1] 66
length(intersect(as.character(ESCC_hit_gene_sets$Fit10More), hsGenes_mmEsGeneSets))
#[1] 111
length(intersect(as.character(ESCC_hit_gene_sets$FitAny), hsGenes_mmEsGeneSets))
#[1] 283


humanNrf2High <- read.csv ("X:/project2020/MustARD/geneSets/hg_NRF2-high_ESCC_geneset.csv")
dim(humanNrf2High)

length(intersect(as.character(humanNrf2High$GeneName), hsGenes_mmEsGeneSets))
#[1] 58
length(intersect(as.character(humanNrf2High$GeneName), as.character(ESCC_hit_gene_sets$FitAny)))
#[1] 42


VennDiagram <- draw.three.list.venndigram(as.character(humanNrf2High$GeneName),hsGenes_mmEsGeneSets, as.character(ESCC_hit_gene_sets$FitAny),
                                          "HSNrf2H",  "mm2hsChIpSeq",  "FittnessESCC")

VennDiagram <- draw.three.list.venndigram(as.character(humanNrf2High$GeneName),hsGenes_mmEsGeneSets, as.character(ESCC_hit_gene_sets$FitAll19),
                                          "HSNrf2H",  "mm2hsChIpSeq",  "FittnessESCC")

grid.draw(VennDiagram$figure)
       


intersect(as.character(humanNrf2High$GeneName), hsGenes_mmEsGeneSets)
write.csv (intersect(as.character(humanNrf2High$GeneName), hsGenes_mmEsGeneSets), file ="x:/project2020/MustARD/geneSets/Set1.csv", row.names = FALSE)
write.csv (intersect(as.character(humanNrf2High$GeneName), as.character(ESCC_hit_gene_sets$FitAny)), file ="x:/project2020/MustARD/geneSets/Set2.csv", row.names = FALSE)
write.csv (intersect(as.character(ESCC_hit_gene_sets$FitAny), hsGenes_mmEsGeneSets), file ="x:/project2020/MustARD/geneSets/Set3.csv", row.names = FALSE)
