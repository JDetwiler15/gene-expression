# Exploring the data from the Allen Brain Institute
# Written by Arianna Krinos, last edits 2 April 2019

pacman::p_load(LMGene,ggplot2,dplyr,limma,edgeR)
setwd("/Users/akrinos/Documents/Spring 2019/Computing the Brain/bgenes/data/human_MTG_gene_expression_matrices_2018-06-14")


##### Set-up / reading in data - use as needed #####
#firsttime = Sys.time()
#exonmatrix = read.csv("human_MTG_2018-06-14_exon-matrix.csv")
#intronmatrix = read.csv("human_MTG_2018-06-14_intron-matrix.csv")
#secondtime = Sys.time()
#totaltime = secondtime - firsttime
#print(totaltime)

#save(exonmatrix, file = "exonmatrix.RData")
#save(intronmatrix, file = "intronmatrix.RData")

#load("exonmatrix.RData")
#load("intronmatrix.RData")

#samplecolumns = read.csv("human_MTG_2018-06-14_samples-columns.csv")
#alldata = bind_rows(exonmatrix, intronmatrix)

# First, let's do it with only exon data, as in Bakken et al. 2018's first cut

# we need to calculation log2(CPM + 1) for each of the genes, and then use the 
# limma package on these transformed data

#' 1) use edgeR to compute the expression data; 
#' actually, step 2) (logexpressdata = glog(alldata, lambda = 0)) can 
#' also be done in this step using log = TRUE

logexpressdatacpm = cpm(exonmatrix, lib.size = NULL, log = TRUE, prior.count = 1) # we're taking log of CPM+1

##### Use edgeR and limma to compute log data and differential expression #####
# Using log2(CPM + 1)

# We use the filtering procedure defined in the limma tutorial
# http://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf

dge <- DGEList(counts=exonmatrix) # creates a handy data structure containing important info about these expression data
keep <- filterByExpr(dge, design = NULL) # filters out the zero expression rows 
dge <- dge[keep,,keep.lib.sizes=FALSE] # only takes the ones you want w reasonable expression
dge <- calcNormFactors(dge) # calculates the TNN normalization 
logCPM <- cpm(dge, log=TRUE, prior.count=1) # calculates the log_2(CPM+1) metric 

fit <- lmFit(logCPM, design = NULL)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))
treat(fit,lfc=log2(2))