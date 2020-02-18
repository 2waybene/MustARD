setwd("x:/project2016/GSA_by_R/")
library(rmarkdown)




# This produce GSA analysis using the original 16 gene sets from David
render("dissecting_gmt_file.rmd",  params = list(dbFile = "cleaningData/Knowledge_based_gene_sets_20140627.gmt"),
       output_file = "Investigate_20140627.html")


# This produce GSA analysis using the original 16 gene sets from David but with slighly modification

# First,
# Adding a new line at the end of gene set file
render("dissecting_gmt_file.rmd",  params = list(dbFile = "cleaningData/genset_1_mod.txt"),
       output_file = "Investigate_data1_mod.html")

##  Next, 
# By adding "na" category
# Modifying the same Sox2xx to V2Sox2

# The results from two lists (dataset 1) seemed to be the same!
render("dissecting_gmt_file.rmd",  params = list(dbFile = "cleaningData/genset_1_mod_2.txt"),
       output_file = "Investigate_data1_mod_2.html")


##  Okay, at this moment, we found out that missing a new line character does NOT seem to cause any trouble!!


##============================================================
# Focus on the second database file
##===========================================================



render("dissecting_gmt_file.rmd",  params = list(dbFile = "cleaningData/Knowledge_based_gene_sets_20161005.gmt"),
       output_file = "Investigate_20161005.html")
## The results are all messed up!!



##============================================================
# two-things can be done here
##===========================================================

##  create a gene set 1 from 20161005.gmt by removing the last gene list genset_1_from_20161005.txt
##  As predicted, the problem are allover!!

render("dissecting_gmt_file.rmd",  params = list(dbFile = "cleaningData/genset_1_from_20161005.txt"),
       output_file = "Investigate_magic_1.html")


##  create a gene set 2 from 20140627.gmt by ADDING the last gene list: genset_1_mod_3.txt
##  As predicted, NO problem here!
render("dissecting_gmt_file.rmd",  params = list(dbFile = "cleaningData/genset_1_mod_3.txt"),
       output_file = "Investigate_magic_2.html")

##====================================================================================================

## Manually appended the gene set on top of 20140627.gmt, Missing the last one
render("dissecting_gmt_file.rmd",  params = list(dbFile = "cleaningData/temp.txt"),
       output_file = "Investigate_data2_added_last_geneset.html")


## Manually modified temp.txt by creating a new line, NO problem
render("dissecting_gmt_file.rmd", params = list(dbFile = "cleaningData/temp2.txt"),
       output_file = "Investigate_v2_w_newline.html")



## Manually modified temp.txt by creating a new line, 
render("dissecting_gmt_file.rmd", params = list(dbFile = "cleaningData/temp3.txt"),
       output_file = "Investigate_v2_w_temp3.html")




render("dissecting_gmt_file.rmd",  params = list(dbFile = "cleaningData/temp4.txt") , 
       output_file = "GSA_4_Chen_v2_cleaned_appended4.html") 



