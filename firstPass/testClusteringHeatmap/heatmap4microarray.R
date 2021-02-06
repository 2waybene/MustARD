##=====================================================================================================================
##  heatmap4microarray.R
##  Credit: https://bioramble.wordpress.com/2015/08/03/heatmaps-part-3-how-to-create-a-microarray-heatmap-with-r/
##=====================================================================================================================


library(ALL)
data(ALL)

class(ALL)
# how much data are we dealing with?
dim(ALL)


# get samples with either no cytogenetic abnormalities (NEG)
# or the BCR-ABL translocation (BCR/ABL)
neg_bcrabl <- ALL$mol.biol %in% c("NEG", "BCR/ABL")
# get indices cancers originating from B-cells
bcell <- grepl("^B", ALL$BT)
# subset the ALL data set
all <- ALL[, bcell & neg_bcrabl]
# adjust the factor levels to reflect the subset
all$mol.biol <- droplevels(all$mol.biol)
all$mol.biol <- relevel(all$mol.biol, ref = "NEG")
# how much data are we left with?
dim(all)
# determine the standard deviation for all genes across the samples
# note that this is essentially an optimized version of
# apply(exprs(all), 1, sd)
library(genefilter)
all_sd <- rowSds(exprs(all))
# get the 200 most variable genes
top200 <- names(sort(all_sd, decreasing = TRUE))[1:200]
all_var <- all[top200, ]

dist_cor <- function(x) {
  as.dist(1 - cor(t(x), method = "pearson"))
}

clus_wd2 <- function(x) {
  hclust(x, method = "ward.D2")
}

library(RColorBrewer)
redblackgreen <- colorRampPalette(c("green", "black", "red"))(n = 100)


class_labels <- ifelse(all_var$mol.biol == "NEG", "grey80", "grey20")


library(gplots)

heatmap.2( exprs(all_var), 
          
           ##===my test data=====
           #log10(y),
           ##====================
          # clustering
          distfun = dist_cor, 
          hclust = clus_wd2,
          # scaling (genes are in rows)
          scale = "row",
          # color
          col = redblackgreen, 
          # labels
          labRow = "", 
    #      ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none")



library(genefilter)
# the shortest interval containing half of the data
# reasonable estimate of the "peak" of the distribution
sh <- shorth(all_sd)
# we take only genes that have a standard deviation
# greater than "sh"
all_sh <- all[all_sd >= sh, ]
# how many genes do we have left?
dim(all_sh)

tt <- rowttests(all_sh, all_sh$mol.biol) 

# use the Benjamini-Hochberg method to adjust 
tt$p.adj <- p.adjust(tt$p.value, method = "BH")
# subset the pre-filtered "all_sh" for genes
# with an adjusted p-value smaller or equal to 0.05
all_sig <- all_sh[tt$p.adj <= 0.05, ]
# how many genes are we left with?
dim(all_sig)

heatmap.2(
  #exprs(all_sig), 
          ##===my test data=====
          log10(y),
          ##====================
          
          # clustering
        #  distfun = dist_cor, 
        #  hclust = clus_wd2,
          # scaling (genes are in rows)
          scale = "row",
          # color
         # col = redblackgreen, 
          col = yellowblackblue, 
          # labels
          labRow = "", 
   #       ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none")

yellowblackblue <- colorRampPalette(c("dodgerblue", "black", "gold"))(n = 100)




