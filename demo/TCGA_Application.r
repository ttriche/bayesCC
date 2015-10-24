# Load the data.  Each file can be downloaded from the TCGA Data Portal via
# https://tcga-data.nci.nih.gov/docs/publications/brca_2012/ 
# (accessed 12/14/2012) 
#
# For formatting reasons, we convert the GE data (BRCA.exp.348.med)
# to a .csv file be for loading. 
# Open the file in Microsoft Excel and save it as .csv.  
# All other files can be read as they are givein in the Data Portal.
#
require(impute) # if running standalone, i.e. for testing
require(gtools) # if running standalone 

# raw downloads from TCGA portal are in inst/extdata/
setwd(system.file("extdata", "", package="bayesCC", mustWork=TRUE))

# cleaned up the following a little bit from the original... --tjt 
#
# mRNA expression
GE <- read.csv("BRCA.exp.348.med.csv.gz", row.names = 1)
names(GE) <- substr(names(GE), 1, 16)

# Impute missing mRNA expression values 
Exp.mat <- impute.knn(as.matrix(GE))$data 

# microRNA expression
miRNA <- read.csv("BRCA.348.precursor.txt.gz", row.names = 1)
names(miRNA) <- substr(names(miRNA), 1, 16)
miRNA.mat <- as.matrix(miRNA[,names(GE)])

# Remove miRNAs with more than 50% missing values
miRNA.mat <- miRNA.mat[rowSums(miRNA.mat == 0) < ncol(miRNA.mat) * 0.5, ]

# protein expression via RPPA 
RPPA <- read.table("rppaData-403Samp-171Ab-Trimmed.txt.gz", head=T, row.names=1)
names(RPPA) <- substr(names(RPPA), 1, 16)
RPPA.mat <- as.matrix(RPPA[,names(GE)])

# DNA methylation
Meth <- read.table("BRCA.Methylation.574probes.802.txt.gz", row.names = 1)
names(Meth) <- substr(names(Meth), 1, 16)
Meth.mat <- as.matrix(Meth[,names(GE)])

allTheSame <- function(...) {
  eq1 <- function(x) base::identical(x, list(...)[[1]])
  return(all(unlist(lapply(list(...), eq1)) == TRUE) )
}
stopifnot(allTheSame(colnames(miRNA.mat), colnames(RPPA.mat), 
                     colnames(Meth.mat), colnames(Exp.mat)))

################### Data processing#############################
X1 <- Exp.mat[apply(Exp.mat, 1, sd) > 1.5, ] # Select most variable genes
X2 <- sqrt(Meth.mat)  # take square root of methylation data ## WTF?! --tjt
X3 <- log(1 + miRNA.mat) # take log of miRNA data (i.e. transcript counts) 
X4 <- scale(RPPA.mat, center = TRUE, scale = TRUE) # center & scale RPPA 

# list of data matrices for input to bayesCC 
X <- list(X1, X2, X3, X4) # X1 = GE, X2 = Meth, X3 = miRNA, X4 = RPPA
BRCAData <- X 
if (FALSE) {
  save(BRCAData, compress="xz", file="../../data/BRCAData.rda")
}

# Perform BCC
#
# Fit the model for K=2, ..., 5 to compare the mean adjusted adherence
# This can take a while (30 minutes to 1 hour for each K)
# To reduce computation time, lower the number of MCMC draws (or parallelize)
#
require(parallel)
Ks <- 2:5 # to reduce computation, don't try an absurd number of clusters!
names(Ks) <- paste0("K", Ks)
Results <- mclapply(Ks, function(k) {
                      bayesCC(BRCAData, K=k, IndivAlpha=T, maxiter=10000)
                    })

# compute adherences
meanAlphaStar <- c()
upperAlphaStar <- c()
lowerAlphaStar <- c()

for (k in Ks) {
  i <- names(Ks)[which(Ks == k)]
  AlphaRaw <- (Results[[i]]$AlphaVec[,2000:10000] - 1 / k) / (1 - 1 / k)
  AlphaStar <- apply(AlphaRaw, 2, mean)
  meanAlphaStar[i] <- mean(AlphaStar)
  upperAlphaStar[i] <- quantile(AlphaStar, 0.975)
  lowerAlphaStar[i] <- quantile(AlphaStar, 0.025)
}

# Plot AlphaStar (mean adjusted adherence) values for each K
plotCI(Ks, meanAlphaStar, ui = upperAlphaStar, li = lowerAlphaStar, 
       xlab = "Number of clusters (K)", ylab = "Mean adjusted adherence") 
# Maximized when K=3

# Fit model for K=3 (can take 30 min - 1 hour)
K3 <- Results[[2]]

# Get alpha values
K3$Alpha

# Make principal component plots of clusters
par(mfrow = c(2, 2))

# mRNA expression
PCs <- prcomp(t(X1))
plot(PCs$x[, 1], PCs$x[, 2], cex = 0, xlab = "PC 1", ylab = "PC 2", 
     main = expression(paste("GE  (", alpha, " = 0.91)")))
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[1]][, 1] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[1]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[1]][, 2] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[1]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[1]][, 3] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[1]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[1]][, 1] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[1]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[1]][, 2] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[1]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[1]][, 3] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[1]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[1]][, 1] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[1]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[1]][, 2] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[1]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[1]][, 3] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[1]][, 3] == 1, 2], pch = 8)

# DNA methylation
m <- 2
PCs <- prcomp(t(X2))
plot(PCs$x[, 1], PCs$x[, 2], cex = 0, xlab = "PC 1", ylab = "PC 2", 
     main = expression(paste("ME  (", alpha, " = 0.69)")))
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)

# miRNA expression
PCs <- prcomp(t(X3))
plot(PCs$x[, 1], PCs$x[, 2], cex = 0, xlab = "PC 1", ylab = "PC 2", 
     main = expression(paste("miRNA  (", alpha, " = 0.56)")))
m <- 3
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "red", 
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "red", 
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "blue", 
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)

m <- 4
PCs <- prcomp(t(X4))
plot(PCs$x[, 1], PCs$x[, 2], cex = 0, xlab = "PC 1", ylab = "PC 2", 
     main = expression(paste("RPPA  (", alpha, " = 0.70)")))
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "black", 
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)

# Find cluster matching matrices The TCGA clusterings were obtained from
# Supplemental Table 1 in the 2012 Nature article 'Comprehensive molecular
# portraits of human breast tumours' The table can be downloaded from (URL):
# http://www.nature.com/nature/journal/v490/n7418/full/nature11412.html
# *accessed 12/14/2012 
# The table should be saved as a 'csv' file before loading
#
Table <- read.csv("Table1Nature.csv", header = TRUE)
Subtype <- Table$Integrated.Clusters..with.PAM50.  ##Comprehensive subtypes

### Match IDs
TableIDs <- Table$Complete.TCGA.ID
TableIDs <- gsub("-", ".", as.character(TableIDs))
namesExp2 <- substr(namesExp, 1, 12)
MatchSub <- match(namesExp2, TableIDs, nomatch = 0)
sum(TableIDs[MatchSub] == namesExp2)
Subtype <- Subtype[MatchSub]

# get and match source-specific subtypes
GESubtype <- Table$PAM50.mRNA
GESubtype <- GESubtype[MatchSub]

MethSubtype <- Table$methylation.Clusters
MethSubtype <- MethSubtype[MatchSub]

miRNASubtype <- Table$miRNA.Clusters
miRNASubtype <- miRNASubtype[MatchSub]

RPPASubtype <- Table$RPPA.Clusters
RPPASubtype <- RPPASubtype[MatchSub]
# Group 'LumA' and 'LumA/B'
RPPASubtype[as.numeric(RPPASubtype) == 3] <- levels(RPPASubtype)[4]

Clusters <- rep(1, 348)
Clusters[K3$Cbest[, 2] == 1] <- 2
Clusters[K3$Cbest[, 3] == 1] <- 3
# ConfusionMat: All
confmat <- matrix(nrow = 4, ncol = 3)
for (i in 1:4) {
  for (j in 1:3) {
    confmat[i, j] <- sum(Subtype == i & Clusters == j)
  }
}

# ConfusionMat: Exp
ClustersGE <- rep(1, 348)
ClustersGE[K3$Lbest[[1]][, 2] == 1] <- 2
ClustersGE[K3$Lbest[[1]][, 3] == 1] <- 3
confmat <- matrix(nrow = 5, ncol = 3)
for (i in 1:5) {
  for (j in 1:3) {
    confmat[i, j] <- sum(as.numeric(GESubtype) == i & ClustersGE == j)
  }
}

# ConfusionMat: Meth
ClustersMeth <- rep(1, 348)
ClustersMeth[K3$Lbest[[2]][, 2] == 1] <- 2
ClustersMeth[K3$Lbest[[2]][, 3] == 1] <- 3
confmat <- matrix(nrow = 5, ncol = 3)
for (i in 1:5) {
  for (j in 1:3) {
    confmat[i, j] <- sum(MethSubtype == i & ClustersMeth == j)
  }
}

# ConfusionMat: miRNA
ClustersmiRNA <- rep(1, 348)
ClustersmiRNA[K3$Lbest[[3]][, 2] == 1] <- 2
ClustersmiRNA[K3$Lbest[[3]][, 3] == 1] <- 3
confmat <- matrix(nrow = 7, ncol = 3)
for (i in 1:7) {
  for (j in 1:3) {
    confmat[i, j] <- sum(miRNASubtype == i & ClustersmiRNA == j)
  }
}

# ConfusionMat: RPPA
ClustersRPPA <- rep(1, 348)
ClustersRPPA[K3$Lbest[[4]][, 2] == 1] <- 2
ClustersRPPA[K3$Lbest[[4]][, 3] == 1] <- 3
confmat <- matrix(nrow = 6, ncol = 3)
for (i in c(1:6)) {
  ### ignore 'Lum A' and 'X' factors
  for (j in 1:3) {
    confmat[i, j] <- sum(as.numeric(RPPASubtype) == i & ClustersRPPA == j)
  }
}


# Heatmaps: Expression
p <- dim(X1)[1]
pre.gene <- t(X1)

subtype <- ClustersGE
n <- length(subtype)

############################ 2.centering datasets ###

gene <- pre.gene - matrix(rep(apply(pre.gene, 2, mean, na.rm = T), each = n), nrow = n)

################################################# Functions for hierarchical clustering ###

row.cluster <- function(data) {
  d.data <- dist(data)
  h.data <- hclust(d.data)
  ordered.data <- data[h.data$order, ]
  return(ordered.data)
}

col.cluster <- function(data) {
  d.data <- dist(t(data))
  h.data <- hclust(d.data)
  ordered.data <- data[, h.data$order]
  return(ordered.data)
}

###################################################### 3.Gene data ordered by subtypes and clustering ###

types <- levels(as.factor(subtype))
K <- length(types)
ordered.x <- c()
sub.n <- c()
for (i in 1:K) {
  typedx <- gene[subtype == types[i], ]
  sub.n <- c(sub.n, dim(typedx)[1])
  ordered.typedx <- row.cluster(typedx)
  ordered.x <- rbind(ordered.x, ordered.typedx)
}


clustered.gene <- col.cluster(ordered.x)
final.x <- clustered.gene


## thresholding
final.x[final.x > 4] <- 4
final.x[final.x < (-4)] <- -4

########################################### 4.Heatmap of gene ###

layout(matrix(c(1:2), 2, 1), heights = c(1, 6), respect = FALSE)

### color key ###
grid <- seq(-4, 4, by = 8 / 100)
par(mar = c(1.5, 0.5, 2, 20))
image(x = grid, y = 1:5, matrix(rep(grid, 5), length(grid), 5), col = greenred(100), 
  axes = FALSE, main = "Color Key", xlab = "", ylab = "", cex.main = 1)
axis(1, seq(-4, 4, by = 1), seq(-4, 4, by = 1))

### title of heatmap ###
mtext("GE with subtypes", side = 4, line = 4, las = 1, cex = 1.7, font = 2)

### Heatmap ###

par(mar = c(5, 0.5, 1, 3))
image(x = 1:n, y = 1:p, final.x[, p:1], axes = FALSE, col = greenred(100), 
      xlab = "", ylab = "")
axis(1, cumsum(sub.n) - sub.n / 2, types, tick = F, cex.axis = 1.3,
     col.axis = "blue", font.axis = 3)
axis(1, c(0, cumsum(sub.n)) + 0.5, rep(" ", K + 1))
abline(v = cumsum(sub.n) + 0.5, col = "white")
mtext("Samples", side = 1, line = 3, cex = 1.5)
mtext("Genes", side = 4, line = 1, cex = 1.5)

# Heatmaps: Meth
hist(X2)
p <- dim(X2)[1]
pre.gene <- t(X2)
subtype <- ClustersMeth
n <- length(subtype)

# 2.centering datasets ###
gene <- pre.gene - matrix(rep(apply(pre.gene, 2, mean, na.rm = T), each = n), nrow = n)

# 3.Gene data ordered by subtypes and clustering ###
types <- levels(as.factor(subtype))
K <- length(types)
ordered.x <- c()
sub.n <- c()
for (i in 1:K) {
  typedx <- gene[subtype == types[i], ]
  sub.n <- c(sub.n, dim(typedx)[1])
  ordered.typedx <- row.cluster(typedx)
  ordered.x <- rbind(ordered.x, ordered.typedx)
}
clustered.gene <- col.cluster(ordered.x)
final.x <- clustered.gene

## thresholding
final.x[final.x > 0.4] <- 0.4
final.x[final.x < (-0.4)] <- -0.4

# 4.Heatmap of gene ###
layout(matrix(c(1:2), 2, 1), heights = c(1, 6), respect = FALSE)

# color key ###
grid <- seq(-0.4, 0.4, by = 0.08 / 100)
par(mar = c(1.5, 0.5, 2, 20))
image(x = grid, y = 1:5, matrix(rep(grid, 5), length(grid), 5), col = greenred(100), 
  axes = FALSE, main = "Color Key", xlab = "", ylab = "", cex.main = 1)
axis(1, seq(-0.4, 0.4, by = 0.1), seq(-0.4, 0.4, by = 0.1))

### title of heatmap ###
mtext("ME with subtypes", side = 4, line = 4, las = 1, cex = 1.7, font = 2)

# Heatmap
par(mar = c(5, 0.5, 1, 3))
image(x = 1:n, y = 1:p, final.x[, p:1], axes = FALSE, col = greenred(100), 
      xlab = "", ylab = "")
axis(1, cumsum(sub.n) - sub.n / 2, types, tick = F, 
     cex.axis = 1.3, col.axis = "blue", font.axis = 3)
axis(1, c(0, cumsum(sub.n)) + 0.5, rep(" ", K + 1))
abline(v = cumsum(sub.n) + 0.5, col = "white")
mtext("Samples", side = 1, line = 3, cex = 1.5)
mtext("Probes", side = 4, line = 1, cex = 1.5)

# Heatmaps: miRNA
p <- dim(X3)[1]
pre.gene <- t(X3)
subtype <- ClustersmiRNA
n <- length(subtype)

# 2.centering datasets ###
gene <- pre.gene - matrix(rep(apply(pre.gene, 2, mean, na.rm = T), each = n), nrow = n)

# 3.Gene data ordered by subtypes and clustering ###
types <- levels(as.factor(subtype))
K <- length(types)
ordered.x <- c()
sub.n <- c()
for (i in 1:K) {
  typedx <- gene[subtype == types[i], ]
  sub.n <- c(sub.n, dim(typedx)[1])
  ordered.typedx <- row.cluster(typedx)
  ordered.x <- rbind(ordered.x, ordered.typedx)
}

clustered.gene <- col.cluster(ordered.x)
final.x <- clustered.gene

## thresholding
final.x[final.x > 2] <- 2
final.x[final.x < (-2)] <- -2

# 4.Heatmap by gene ###
layout(matrix(c(1:2), 2, 1), heights = c(1, 6), respect = FALSE)

# color key ###
grid <- seq(-2, 2, by = 5 * 0.08 / 100)
par(mar = c(1.5, 0.5, 2, 20))
image(x = grid, y = 1:5, matrix(rep(grid, 5), length(grid), 5), col = greenred(100), 
  axes = FALSE, main = "Color Key", xlab = "", ylab = "", cex.main = 1)
axis(1, seq(-2, 2, by = 5 * 0.1), seq(-2, 2, by = 5 * 0.1))

# title of heatmap ###
mtext("miRNA with subtypes", side = 4, line = 4, las = 1, cex = 1.7, font = 2)

# Heatmap ###
par(mar = c(5, 0.5, 1, 3))
image(x = 1:n, y = 1:p, final.x[, p:1], axes = FALSE, col = greenred(100), 
      xlab = "", ylab = "")
axis(1, cumsum(sub.n) - sub.n / 2, types, tick = F, 
     cex.axis = 1.3, col.axis = "blue", font.axis = 3)
axis(1, c(0, cumsum(sub.n)) + 0.5, rep(" ", K + 1))
abline(v = cumsum(sub.n) + 0.5, col = "white")
mtext("Samples", side = 1, line = 3, cex = 1.5)
mtext("miRNAs", side = 4, line = 1, cex = 1.5)

# Heatmaps: RPPA
p <- dim(X4)[1]
pre.gene <- t(X4)
subtype <- ClustersRPPA
n <- length(subtype)

# 2.centering datasets ###
gene <- pre.gene - 
        matrix(rep(apply(pre.gene, 2, mean, na.rm = T), each = n), nrow = n)

# 3.Gene data ordered by subtypes and clustering ###
types <- levels(as.factor(subtype))
K <- length(types)
ordered.x <- c()
sub.n <- c()
for (i in 1:K) {
  typedx <- gene[subtype == types[i], ]
  sub.n <- c(sub.n, dim(typedx)[1])
  ordered.typedx <- row.cluster(typedx)
  ordered.x <- rbind(ordered.x, ordered.typedx)
}
clustered.gene <- col.cluster(ordered.x)
final.x <- clustered.gene

## thresholding
final.x[final.x > 2] <- 2
final.x[final.x < (-2)] <- -2

# 4.Heatmap of gene ###
layout(matrix(c(1:2), 2, 1), heights = c(1, 6), respect = FALSE)

### color key ###
grid <- seq(-2, 2, by = 5 * 0.08 / 100)
par(mar = c(1.5, 0.5, 2, 20))
image(x = grid, y = 1:5, matrix(rep(grid, 5), length(grid), 5),
      col = greenred(100), axes = FALSE, main = "Color Key", 
      xlab = "", ylab = "", cex.main = 1)
axis(1, seq(-2, 2, by = 5 * 0.1), seq(-2, 2, by = 5 * 0.1))

### title of heatmap ###
mtext("RPPA with subtypes", side = 4, line = 4, las = 1, cex = 1.7, font = 2)

# Heatmap ###
par(mar = c(5, 0.5, 1, 3))
image(x = 1:n, y = 1:p, final.x[, p:1], axes = FALSE, col = greenred(100), 
      xlab = "", ylab = "")
axis(1, cumsum(sub.n) - sub.n / 2, types, tick=F, 
     cex.axis=1.3, col.axis="blue", font.axis=3)
axis(1, c(0, cumsum(sub.n)) + 0.5, rep(" ", K + 1))
abline(v = cumsum(sub.n) + 0.5, col = "white")
mtext("Samples", side = 1, line = 3, cex = 1.5)
mtext("RPPA", side = 4, line = 1, cex = 1.5)
