setwd("X:/project2020/MustARD/learningFromBigWigs/Behan_nature_CRISPR-Cas9/")

ESCC.cell.lines <- read.csv("esophagus_CellLines.csv")
fit.genes <- read.csv("fitness_genes_all.csv")
dim(fit.genes)
ESCC <- ESCC.cell.lines$CMP_id[ ESCC.cell.lines$CancerType %in% "Esophageal Squamous Cell Carcinoma"]


length(which(colnames (fit.genes) %in% ESCC))
#19, missing "SIDM00249"

ESCC[-which(ESCC %in% colnames (fit.genes))]
#[1] SIDM00249


ESCC.fit.genes <- fit.genes[, (which(colnames (fit.genes) %in% ESCC))]

dim(ESCC.fit.genes)
apply (ESCC.fit.genes[,-1], SUM)

ESCC.genes <- ESCC.fit.genes[,1][which(rowSums (ESCC.fit.genes[,-1]) == 18)]


