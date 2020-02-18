setwd("X:/project2020/MustARD/GSA_analysis/htmlReport/")
library(rmarkdown)
library(knitr)


render("x:/project2016/woychik-ma-rnaseq/scripts/comparing_DEGs_v1.rmd", 
       output_file = "x:/project2016/woychik-ma-rnaseq/comparison/comparing_DEGs_v1.html")

