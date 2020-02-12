#=====================================================================
##  ESCC_hit_genes.R
##  Paper: https://www.nature.com/articles/s41586-019-1103-9#Sec44  
##====================================================================

##  Download supplemental table 1 & 2

setwd("X:/project2020/MustARD/learningFromBigWigs/Behan_nature_CRISPR-Cas9/")


##===============================================================================
##  27 Esophagus cell lines, of which 20 are Esophageal Squamous Cell Carcinoma
##===============================================================================


ESCC.cell.lines <- read.csv("esophagus_CellLines.csv")
fit.genes <- read.csv("fitness_genes_all.csv")
dim(fit.genes)

ESCC <- ESCC.cell.lines$CMP_id[ ESCC.cell.lines$CancerType %in% "Esophageal Squamous Cell Carcinoma"]
length(which(colnames (fit.genes) %in% ESCC))
#19, missing "SIDM00249"
ESCC[-which(ESCC %in% colnames (fit.genes))]
#[1] SIDM00249
##  ONLY 19 out of 20 have hit genes

length(fit.genes$CMP.id)
#[1] 7470

ESCC.fit.genes <- fit.genes[ , c(1, (which(colnames (fit.genes) %in% ESCC)))]
dim(ESCC.fit.genes)
#[1] 7470   20

length(which(rowSums (ESCC.fit.genes[,-1]) == 19))
#[1] 374
length(which(rowSums (ESCC.fit.genes[,-1]) == 18))
#[1] 187
length(which(rowSums (ESCC.fit.genes[,-1]) >= 18))
#[1] 561

length(which(rowSums (ESCC.fit.genes[,-1]) >= 10))
#[1] 1353
length(which(rowSums (ESCC.fit.genes[,-1]) >= 15))
#[1] 945


ESCC.genes.all <- ESCC.fit.genes[,1][which(rowSums (ESCC.fit.genes[,-1]) == 19)]
ESCC.genes.18.more <- ESCC.fit.genes[,1][which(rowSums (ESCC.fit.genes[,-1]) >= 18)]
ESCC.genes.15.more <- ESCC.fit.genes[,1][which(rowSums (ESCC.fit.genes[,-1]) >= 15)]
ESCC.genes.10.more <- ESCC.fit.genes[,1][which(rowSums (ESCC.fit.genes[,-1]) >= 10)]
ESCC.genes.1.more <- ESCC.fit.genes[,1][which(rowSums (ESCC.fit.genes[,-1]) >= 1)]

ESCC_hit_gene_sets <- list ("FitAll19" = ESCC.genes.all , "Fit18More" = ESCC.genes.18.more  , 
                            "Fit15More" = ESCC.genes.15.more, "Fit10More" = ESCC.genes.10.more, 
                            "FitAny" = ESCC.genes.1.more)
save (ESCC_hit_gene_sets , file = "Behan_CRISPR_ESCC_hit_genes.rda")
load("Behan_CRISPR_ESCC_hit_genes.rda")
ESCC_hit_gene_sets$FitAll19
ESCC_hit_gene_sets$Fit18More
ESCC_hit_gene_sets$Fit15More
ESCC_hit_gene_sets$Fit10More
ESCC_hit_gene_sets$FitAny
##===============================================================================
##  Compare to human NRF2 high ESCC geneset
##===============================================================================


humanNrf2High <- read.csv ("X:/project2020/MustARD/geneSets/hg_NRF2-high_ESCC_geneset.csv")
dim(humanNrf2High)

length(intersect(as.character(humanNrf2High$GeneName), as.character(ESCC_hit_gene_sets$Fit10More)))
#[1] 6


length(intersect(as.character(humanNrf2High$GeneName), as.character(ESCC.genes.1.more)))


length(intersect(as.character(humanNrf2High$GeneName),as.character(ESCC.fit.genes[,1])))
intersect(as.character(humanNrf2High$GeneName),as.character(ESCC.fit.genes[,1]))
length(as.character(ESCC.fit.genes[,1]))

ESCC.fit.genes[which(ESCC.fit.genes[,1] == "KEAP1"),]


