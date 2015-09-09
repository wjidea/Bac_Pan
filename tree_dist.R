# load dependencies
library(ape)
library(ggplot2)
library(ggtree)
library(gridExtra)
#vignette("ggtree", package = "ggtree")


# read in data
setwd("/Users/wjidea/GoogleDrive/Graduate_Study/4_Collaborations/2014_Acidovorax_Quan/Results/trees_mod/")
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
PhyloTree <- read.tree("/Users/wjidea/GoogleDrive/Graduate_Study/4_Collaborations/2014_Acidovorax_Quan/Results/core_gene.tre")
plot(PhyloTree)
#ggplot(PhyloTree, aes(x, y)) + geom_tree() + theme_tree() + xlab("") + ylab("")
ggtree(PhyloTree) %>% add_legend()
p <- ggtree(PhyloTree) + geom_tiplab()
p
annotation_clade2(p, "QH1", "Aa99-2","selected clade1") %>%
  annotation_clade2(p, "QHB1", "Cat98-1", "selected clade2")
data(chiroptera)
gzoom(chiroptera, grep("Plecotus", chiroptera$tip.label))
gzoom(PhyloTree, c(1:9))
gzoom(PhyloTree, c(10:16))
