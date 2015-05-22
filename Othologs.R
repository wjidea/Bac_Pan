#!/usr/bin/env Rscript
# ABSTRACT: Create R plots
# PODNAME: create_plots.R

# Take the output files from the pan genome pipeline and create nice plots.
setwd("/Users/wjidea/Google\ Drive/Graduate_Study/Collaborations/Acidovorax/Scripts/data")

mydata = read.table("number_of_new_genes.Rtab")
boxplot(mydata, data=mydata, main="Number of new genes",xlab="Number of genomes",
        ylab="Number of genes",varwidth=TRUE, ylim=c(0,max(mydata)), outline=FALSE)

mydata = read.table("number_of_conserved_genes.Rtab")
boxplot(mydata, data=mydata, main="Number of conserved genes",
        xlab="Number of genomes",ylab="Number of genes",varwidth=TRUE, 
        ylim=c(0,max(mydata)), outline=FALSE)

mydata = read.table("number_of_genes_in_pan_genome.Rtab")
boxplot(mydata, data=mydata, main="Number of genes in the pan-genome",
        xlab="Number of genomes", ylab="Number of genes",varwidth=TRUE, 
        ylim=c(0,max(mydata)), outline=FALSE)

mydata = read.table("number_of_unique_genes.Rtab")
boxplot(mydata, data=mydata, main="Number of unique genes",
        xlab="Number of genomes", ylab="Number of genes",varwidth=TRUE, 
        ylim=c(0,max(mydata)), outline=FALSE,)

mydata = read.table("blast_identity_frequency.Rtab")
plot(mydata,main="Number of blastp hits with different percentage identity",  
     xlab="Blast percentage identity", ylab="No. blast results")