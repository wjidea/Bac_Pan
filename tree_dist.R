# load dependencies
library(ape)
library(ggplot2)

# read in data
setwd("/Users/wjidea/Google\ Drive/Graduate_Study/4_Collaborations/2014_Acidovorax_Quan/Results/trees_mod/")
files <- dir("./", pattern = ".*tre")
trees_list <- lapply(files, read.tree)

Dist2Genome <- vector()
for (i in 1:length(trees_list)){
  Dist2Genome[i] <- dist.topo(trees_list[[i]], PhyloTree, method="PH85")
}
hist(Dist2Genome)
rug(Dist2Genome)
summary(Dist2Genome)
which.min(Dist2Genome)
which(Dist2Genome == 9)

par(mfrow=c(2,2))
head(trees_list)

for (i in 1:100){
  plot(trees_list[[i]],main = files[i],cex.main=0.8, cex=0.8)
}
iaa = 9

plot(trees_list[[iaa]],main = files[iaa],cex.main=0.8, cex=0.8)

# Sandbox
PhyloTree <- read.tree("/Users/wjidea/Google\ Drive/Graduate_Study/4_Collaborations/2014_Acidovorax_Quan/Results/core_gene.tre")
