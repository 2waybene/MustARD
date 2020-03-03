# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Mon Mar 2 20:43:51 EST 2020

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE36133", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL15308", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXX3XX3XXXXXXXXXXXXXXXX5XXXXXX",
               "XXXXXXXXXXXXXXXXXXXX3XXXXXXXX405XXXXXXXXXX5XXXXXXX",
               "XXXXXXX1XXX0XXXXXXXX0211111XXXXX1XXXXXXXXXXXX11111",
               "1XXXX0XX4XXXXXXXXXXXXX4XXXXXXXXXXXXXXXXXXXXXXXX4X3",
               "X5XXXXXX4XXXXXX0315X00XX004XXXXXXXXXXXXXXXXXXXXXXX",
               "4XXXXXXX0XXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXX2XXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX4XXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXX422X4XXXXX4XXXXXXXX324",
               "0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXX0XXXXXXXXXX11121311004031010302100505010005031",
               "3412131130013X1100133100141311041000303XX05012X143",
               "11555X21X5X213101XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXX3XXXXXXXXXXXXXXXXXXXXXXXXXX004",
               "334XXXXXXXXXXXXXXXXXX1XXXXX1XXXXXX1XXXXXXXXXX0XXXX",
               "XXXXX4XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXX4XXXXXXXXXXXXXXXXX1XXX4XXXXXXXXXX4XXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXX")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G5-G0, G1-G0, G2-G1, G3-G2, G4-G3, G5-G4, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

data.lung <- list ("dt"=gset, "dn"= design, "f"= fit, "f2"= fit2)
save (data.lung, file = "/Users/li11/myGit/MustARD/ESCC/workingWithMicroarray_GSE36133_CCLE/lung_cancer.rda")

tT <- topTable(fit2, adjust="fdr", sort.by="B")


tT <- topTable(fit2, adjust="fdr", sort.by="B", number=25000)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","ORF"))
head(tT)
write.table(tT, file=stdout(), row.names=F, sep="\t")


################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE36133", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL15308", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXX3XX3XXXXXXXXXXXXXXXX5XXXXXX",
               "XXXXXXXXXXXXXXXXXXXX3XXXXXXXX405XXXXXXXXXX5XXXXXXX",
               "XXXXXXX1XXX0XXXXXXXX0211111XXXXX1XXXXXXXXXXXX11111",
               "1XXXX0XX4XXXXXXXXXXXXX4XXXXXXXXXXXXXXXXXXXXXXXX4X3",
               "X5XXXXXX4XXXXXX0315X00XX004XXXXXXXXXXXXXXXXXXXXXXX",
               "4XXXXXXX0XXXXXXXXXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXX2XXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX4XXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXX422X4XXXXX4XXXXXXXX324",
               "0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXX0XXXXXXXXXX11121311004031010302100505010005031",
               "3412131130013X1100133100141311041000303XX05012X143",
               "11555X21X5X213101XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXX3XXXXXXXXXXXXXXXXXXXXXXXXXX004",
               "334XXXXXXXXXXXXXXXXXX1XXXXX1XXXXXX1XXXXXXXXXX0XXXX",
               "XXXXX4XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXX4XXXXXXXXXXXXXXXXX1XXX4XXXXXXXXXX4XXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXX")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  # set group names

# eliminate samples marked as "X"
sel <- which(sml != "GX")
sml <- sml[sel]
gset <- gset[ ,sel]

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("ADC","ADSC","LCC","NSCLC","SCC","Others")
table(fl)

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf","#f2cb98","#dcdaa5","#dff4e4","#f4dff4", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE36133", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")



dim(ex)
