# Exploring the data from the Allen Brain Institute
# Written by Arianna Krinos, last edits 9 April 2019

##### Needed functions and setup ##### 
pacman::p_load(LMGene,ggplot2,dplyr,limma,edgeR,tsne,dplyr,irlba,Rtsne,devtools,igraph,ggplot2)
devtools::install_github("JinmiaoChenLab/Rphenograph")
setwd("/Users/akrinos/Documents/Spring 2019/Computing the Brain/bgenes/data/human_MTG_gene_expression_matrices_2018-06-14")

CalcBrokenStick <- function(x) {
  m <- 0
  out <- matrix(NA, ncol = 3, nrow = length(x))
  colnames(out) <- c("pc", "var_pct", "bstick_thresh")
  for (i in 1:length(x)) {
    for (k in i:length(x)) {
      m <- m + ((1 / length(x)) * (1 / k))
    }
    out[i, ] <- c(i, (x[i] / sum(x)) * 100, m * 100)
    m <- 0
  }
  pc.var.order <- order(out[, "var_pct"], decreasing = TRUE)
  out <- out[pc.var.order, ]
  out[, "bstick_thresh"] <- sort(out[, "bstick_thresh"], decreasing = TRUE)
  return(out)
}


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
dge_fortsne = dge # we want to save the matrix here for use by tSNE
dge <- calcNormFactors(dge) # calculates the TNN normalization 
logCPM <- cpm(dge, log=TRUE, prior.count=1) # calculates the log_2(CPM+1) metric 

fit <- lmFit(logCPM, design = NULL)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))
treat(fit,lfc=log2(2))

##### Compute PCA using irlba (a quickie solution) #####
pc_dge_2 <- irlba(t(dge$counts)[-1,], nv=20, nu=0, center=colMeans(t(dge$counts)[-1,])) #, right_only=TRUE)
newmatrix = t(dge$counts)[-1,] %*% pc_dge_2$v

##### Data processing steps from AllenInstitute github #####
# Keep cells with at least M genes detected
min.genes <- 1000
keep.samp <- which(apply(expr, 2, function(x) sum(x > 0)) > min.genes)
# Keep genes with some variability
keep.genes <- which(apply(expr[, keep.samp], 1, sd) > 0)
# Subset and log-transform counts per million
expr.cnt <- expr[keep.genes, keep.samp]
expr.cpm <- sweep(expr.cnt, 2, colSums(expr.cnt), "/") * 1e6
expr.log <- log2(expr.cpm + 1)


##### Pull in PhenoGraph results from Python source (rename file) #####
phenograph = read.csv("../../code/communities_output.csv")

##### Run PhenoGraph & tSNE in R using package from https://github.com/AllenInstitute/RNAseq_cluster #####
pal1 <- c("#d77a7f", "#8eda48", "#7340cd", "#d6c847", "#ce4cc5", "#64db8e", 
          "#432876", "#509140", "#7171cd", "#d1863a", "#79acd9", "#d24530", 
          "#6dc7b7", "#d23e70", "#c6d394", "#8d3870", "#827f38", "#cd90cb", 
          "#3a4e32", "#c9c5c6", "#3e263b", "#ae8875", "#556983", "#753627")
palette(pal1)

nn.num <- 15  # Number of nearest cells to compare in building graph
phenographres_pca = Rphenograph::Rphenograph(newmatrix, k = nn.num)

tsne1 <- Rtsne(newmatrix, perplexity = 20)$Y

## Reduce data by some amount 

meta <- samplecolumns
expr <- exonmatrix
genes <- read.csv("/Users/akrinos/Documents/Spring 2019/Computing the Brain/bgenes/data/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_genes-rows.csv") # gene data

# modify the number of genes you save 
min.genes <- 1000 # decide how many genes need to be detected
keep.samp <- which(apply(expr, 2, function(x) sum(x > 0)) > min.genes) # select only the cell samples that have this level of expression 

# Keep genes with some variability
keep.genes <- which(apply(expr[, keep.samp], 1, sd) > 0) # from the exon matrix, select the columns from the cells you wish to play with

# Subset and log-transform counts per million
expr.cnt <- expr[keep.genes, keep.samp]
expr.cpm <- sweep(expr.cnt, 2, colSums(expr.cnt), "/") * 1e6
expr.log <- log2(expr.cpm + 1)
meta.subset <- droplevels(meta[keep.samp, ])
rownames(meta.subset) <- colnames(expr.cnt)
genes.subset <- droplevels(genes[keep.genes, ])

# this can be equal to all of meta, or droplevels (meta[])
meta.subset <- droplevels(meta[keep.samp, ]) #meta - so far I've only kept these, will do more processing in future

## Now lets look at the genes instead of the samples

expr.mean <- apply(expr.cpm, 1, mean) # get the mean CPMs of the gene expression matrix
expr.mean.bin <- cut(expr.mean, breaks = 10^seq(-5, 5, 0.1), include.lowest = TRUE) # bin up the expression data
expr.sd <- apply(expr.cpm, 1, sd) # get the standard deviation of the CPM data
expr.cv <- expr.sd / expr.mean # get the coefficient of variation of the data (sd/mean)
expr.cv.z.list <- tapply(expr.cv, expr.mean.bin, scale) # creates a pairing of the coefficients of variation for the binned expression data
expr.cv.z <- do.call("rbind", expr.cv.z.list) # creates a combined non-list (matrix) of this binning
expr.cv.z <- expr.cv.z[match(names(expr.mean), rownames(expr.cv.z)), 1]
expr.cv.z[expr.cv.z > 8] <- 8 # scale down to 8 if expression is greater than 8

# Select z-score threshold for variable genes
z.thresh <- 2
top.var.genes <- which(expr.cv.z > z.thresh & expr.mean > 0.1)
par(mfrow = c(1, 2))
plot(expr.mean, expr.cv, cex = 0.5, log = "xy", 
     xlab = "Average expression", 
     ylab = "Coefficient of variation (CV)", 
     main = paste(length(top.var.genes), "variable genes"),pch=20)
points(expr.mean[top.var.genes], expr.cv[top.var.genes], col = "orange", pch = 20)
plot(expr.mean, expr.cv.z, cex = 0.5, log = "x", 
     xlab = "Average expression", 
     ylab = "CV (z-score)", pch = 20,
     main = paste(length(top.var.genes), "variable genes"))
points(expr.mean[top.var.genes], expr.cv.z[top.var.genes], col = "orange", pch = 20)

## Full PCA procedure taken from https://github.com/AllenInstitute/RNAseq_cluster/blob/master/R/cluster_cells.Rmd
# Center and z-score genes across samples
expr.scaled <- scale(expr.cpm) # this takes forever
# Calculate PCA using variable genes
pca1 <- prcomp(expr.scaled[top.var.genes, ]) # perform PCA on the scaled expression data, selecting for the top variable genes

## This is the up to 20 principal components thing 
# Find PCs with more variance explained than broken stick distribution
brstick1 <- CalcBrokenStick(pca1$sdev^2)
max.pcs <- 20 # maximum number of principal components to use
sig.pcs <- NULL
for (i in 1:max.pcs) {
  var.pct <- brstick1[i, "var_pct"]
  var.thresh <- brstick1[i, "bstick_thresh"]
  if (var.pct > var.thresh) {
    sig.pcs <- c(sig.pcs, brstick1[i, "pc"])
  } else {
    break
  }
}

# Select PCs to retain for further analysis
keep.pcs <- sig.pcs
pca.var <- summary(pca1)$importance[1, ]^2
par(mfrow = c(1,1))
plot(pca.var[1:max.pcs], type = "b", xlab = "PC", ylab = "Explained variance")
abline(v = max(sig.pcs), col = "light blue")
text(x = max(sig.pcs), y = 0.9 * max(pca.var), labels = "# significant PCs", pos = 2)
abline(v = max(keep.pcs), col = "blue")
text(x = max(keep.pcs), y = 0.8 * max(pca.var), labels = "Keep PCs", pos = 2)
expr.pcs <- pca1$rotation[, keep.pcs]

# Graph-based clustering (Jaccard/Louvain)
nn.num <- 15  # Number of nearest cells to compare in building graph
graph.cl <- Rphenograph::Rphenograph(expr.pcs, k = nn.num)
plot.lab <- paste(length(unique(membership(graph.cl[[2]]))), "clusters; ",
                  round(modularity(graph.cl[[2]]), 2), "modularity")
plot(graph.cl[[1]], vertex.size = 6, vertex.label = NA, edge.color = NA, #edge.color = "grey80", 
     vertex.color = pal1[membership(graph.cl[[2]])], main = plot.lab)

tsne2 <- Rtsne(expr.pcs, perplexity = 20)$Y


select.cl <- "graph"
meta.subset$cluster <- as.character(membership(graph.cl[[2]]))
fills = factor(meta.subset$cluster, levels = c(as.character(1:length(unique(meta.subset$cluster)))))

ggplot(data.frame(tsne2), aes(x = X1, y = X2, fill = fills)) + 
  xlab("tSNE1") +  ylab("tSNE2") + 
  ggtitle(paste(length(unique(meta.subset$cluster)), select.cl, "clusters")) + 
  geom_point(colour = "black", pch = 21) + 
  scale_fill_discrete(name = "Cluster") + theme_classic(base_size = 16)

##### This next section is same as above but mostly my own code ##### 

par(mfrow = c(1,1))
plot.lab <- paste(length(unique(membership(phenographres_pca[[2]]))), "clusters; ",
                  round(modularity(phenographres_pca[[2]]), 2), "modularity")
plot(phenographres_pca[[1]], vertex.size = 6, vertex.label = NA, edge.color = "grey80", 
     vertex.color = pal1[membership(phenographres_pca[[2]])], main = plot.lab)

meta.subset$cluster <- as.character(membership(phenographres_pca[[2]]))
select.cl <- "graph"
fills = factor(meta.subset$cluster, levels = c(as.character(1:length(unique(meta.subset$cluster)))))
ggplot(data.frame(tsne1), aes(x = X1, y = X2, fill = fills)) + 
         xlab("tSNE1") +  ylab("tSNE2") + 
         ggtitle(paste(length(unique(meta.subset$cluster)), select.cl, "clusters")) + 
         geom_point(colour = "black", pch = 21) + 
         scale_fill_discrete(name = "Cluster") + theme_classic(base_size = 16)

plot(tsne1, col = phenographres_pca[[2]], xlab = "tSNE1", ylab = "tSNE2")

##### Run tsne on entire dataset #####

#tsneresults = tsne(dge_fortsne)
#save(tsneresults, file = "tsneresults_8Apr19.RData")

##### Run tsne on PCA reduced dataset #####
tsneresults = tsne(dge_fortsne)
save(tsneresults, file = paste0("tsneresults_pca_", Sys.date(), ".RData"))

###### Visualize results of tsne #####

load("tsneresults_8Apr19.RData")

overall = cbind(tsneresults, phenograph[,2])
colnames(overall) = c("tsne1", "tsne2", "phenocluster")
overall = data.frame(overall)

ggplot(overall) + geom_point(aes(x = tsne1, y = tsne2, fill = factor(phenocluster)), pch = 21, size = 3) + 
  theme_classic(base_size = 18) + scale_fill_discrete(name = "PhenoGraph community")

dge_mess = t(data.frame(rbind(dge$counts[,2:ncol(dge$counts)], samplecolumns$cluster)))

colnames(dge_mess) = c(rownames(dge$counts)[1:nrow(dge$counts)], "cluster")
dge_mess = data.frame(dge_mess) %>% group_by(cluster) %>% summarize_if(is.numeric, funs(mean))#,median,mode,max,min))

##### Look into ways to perform hierarchical clustering #####
library(cluster)
DM = as.matrix(dist(dge$counts))
HC = hclust(as.dist(DM), method="single")
