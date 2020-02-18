library(ROntoTools)
kpg <- keggPathwayGraphs("hsa")
# to update the pathway cache for human run:
# kpg <- keggPathwayGraphs("hsa", updateCache = TRUE)
# this is time consuming and depends on the available bandwith.
head(names(kpg))
kpg[["path:hsa04110"]]
head(nodes(kpg[["path:hsa04110"]]))
head(edges(kpg[["path:hsa04110"]]))


kpn <- keggPathwayNames("hsa")
# to update the pathway cache for human run:
# kpn <- keggPathwayNames("hsa", updateCache = TRUE)
# this is time consuming and depends on the available bandwidth.
head(kpn)


# The table contains the expression fold change and signficance of each
# probe set in peripheral blood mononuclear cells (PBMC) from 12 MS patients
# and 15 controls.
load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
head(top)
# select differentially expressed genes at 1% and save their fold change in a
# vector fc and their p-values in a vector pv
fc <- top$logFC[top$adj.P.Val <= .01]
names(fc) <- top$entrez[top$adj.P.Val <= .01]
pv <- top$P.Value[top$adj.P.Val <= .01]
names(pv) <- top$entrez[top$adj.P.Val <= .01]
# alternativly use all the genes for the analysis
# NOT RUN:
# fc <- top$logFC
# names(fc) <- top$entrez
# pv <- top$P.Value
# names(pv) <- top$entrez
# get the reference
ref <- top$entrez
# load the set of pathways
kpg <- keggPathwayGraphs("hsa")
# set the beta information (see the citated documents for meaning of beta)
kpg <- setEdgeWeights(kpg)
# inlcude the significance information in the analysis (see Voichita:2012
# for more information)
# set the alpha information based on the pv with one of the predefined methods
kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)
# perform the pathway analysis
# in order to obtain accurate results the number of boostraps, nboot, should
# be increase to a number like 2000
peRes <- pe(fc, graphs = kpg, ref = ref, nboot = 100, verbose = TRUE)
# obtain summary of results
head(Summary(peRes))

