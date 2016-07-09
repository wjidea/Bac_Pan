# load dependencies
library(ape)
library(ggplot2)
library(ggtree)
library(gridExtra)
library("colorspace")
#vignette("ggtree", package = "ggtree")

# read in data
setwd("/Users/wjidea/Desktop/alnTre 2/trees/")
files <- dir("./", pattern = ".*tre")
trees_list <- lapply(files, function(x) ggtree::read.tree(x, tree.names = paste(x)))

cls <- list(maize=c("AA19", "AA78", "AA81", "Aa99"),
            Turf1=c("INV", "KL3", "URI9", "QHB1", "MDB1", "NCT3"),
            Turf2=c("COLB", "INDB", "SF12", "QH1", "MOR","SH7", "MD5"),
            Others=c("AC01", "Cat9", "Sa2"))
lenL <- length(trees_list) - 2
j = 0
treeL <- list()
for (i in 1:lenL){
  j = j + 1
  trees_list[[i]]$tip.label <- gsub("_1_[0-9]{2,5}$", "",trees_list[[i]]$tip.label, perl = TRUE)
  tree <- groupOTU(trees_list[[i]], cls)
  plotTree <- ggtree(tree) + geom_tiplab(aes(color = group)) +
              scale_color_manual(values=c("green", "red", "blue", "black")) +
    ggtitle(paste(files[i]))
  treeL[[j]] <- plotTree
  if (j %% 9 == 0 ){
    finalOut <- arrangeGrob(treeL[[1]], treeL[[2]],treeL[[3]],treeL[[4]],
                            treeL[[5]],treeL[[6]],treeL[[7]],treeL[[8]],
                            treeL[[9]], ncol=3)
    ggsave(filename = paste("trees_group",i,".png",sep=""), device = "png",
           plot = finalOut, width = 15, height = 15, 
           path = "/Volumes/SP128g/TEMP/treeFiles")
    j <- 0
    treeL <- list()
  }
}

# draw spectified tree
treeName <- "group_4637.tre"
specTree <- read.tree(treeName)
specTree <- groupOTU(specTree, cls)
ggtree(specTree) + geom_tiplab(aes(color = group)) +
  scale_color_manual(values=c("green", "red", "blue", "black")) +
  ggtitle(paste(specTree))



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
