############### Install packages and load in data ############################

# install and load necessary packages
doInstall <- FALSE  # Change to FALSE if you don't want packages installed.
toInstall <- c("ggplot2","magrittr", "dplyr", "vegan", "RColorBrewer")
if(doInstall){install.packages(toInstall, repos = "http://cran.r-project.org")}
lapply(toInstall, library, character.only = TRUE)


# Load data 
data_genes <- read.csv("../2015_05_core_pan_genome_1.csv",header = TRUE, 
                       na.strings = c("","NA"),stringsAsFactors = FALSE)

# data cleanning
data_sub <- select(data_genes, c(Gene, No..isolates, AA38,AA78.5,Aa99.2,COLB1,
                                 Cat98.1, INDB2, INV, KL3, MDB1, MOR, NCT3, QH1, 
                                 QHB1, SF12, SH7, Sa2))
data_sub <- rename(data_sub, Isolate_number = No..isolates, AA78 = AA78.5, 
                   Aa99 = Aa99.2, Cat98 = Cat98.1 )
head(data_sub)

# data transform
# str(data_sub)
# data_sub[1:2,1:10]
 
for (j in 3:ncol(data_sub)){
  for (i in 1:nrow(data_sub)) {
    if (!(is.na(data_sub[i,j]))) {
      data_sub[i,j] <- data_sub[i,1]
    }
  }
}
# Archive the transformed data
work_data <- data_sub

# function to calculate the intersect among groups
x1 <- intersect(work_data[,5],intersect(work_data[,3],work_data[,4]))
length(x1)

core_genes <- function(s1 = 3,s2 = 4,s3 = 5){
  core_set <- intersect(work_data[,s1],intersect(work_data[,s2],work_data[,s3]))
  return(length(core_set))
}

# resampling calculate core genome from 3 randomly sampled strains
# 1. maize pathogens (3 strains) AA78, Aa99, AA38
# 2. all turf pathogens (12 strains)
# 3. group1 turf (KL3, INV, QHB1 MDB1, NCT3)
# 4. group2 turf (SH7, QH1, MOR, COLB1, INDB2, SF12)
# 5. all turf and maize pathogens (15 strains)
# 6. all strains 15 strains + Cat98-1

#colnames(work_data)[3:18]
turf_strains_index <- c(6,8:18)
colnames(work_data)[turf_strains_index]
maize_strains_index <- c(3:5)
colnames(work_data)[maize_strains_index]
turf_g1_index <- c(9,10,11,13,15)
colnames(work_data)[turf_g1_index]
turf_g2_index <- c(17,14,12,6,8,16)
colnames(work_data)[turf_g2_index]
turf_maize_index <- c(3:6,8:18)
colnames(work_data)[turf_maize_index]
all_index <- c(3:18)
colnames(work_data)[all_index]
g1_maize_index <- c(9,10,11,13,15,3:5)
colnames(work_data)[g1_maize_index]
g2_maize_index <- c(17,14,12,6,8,16,3:5)
colnames(work_data)[g2_maize_index]
# task 1 maize pathogens (3 strains) AA78, Aa99, AA38
T1 <- core_genes()
T1

# task 2 - all turf pathogens (12 strains)
store_coreset_2 <- c()
i <- 0
for (i in 1:1000){
  sub_strains <- turf_strains_index[sample(1:length(turf_strains_index),3,replace = FALSE)]
  store_coreset_2[i] <- core_genes(sub_strains[1],sub_strains[2],sub_strains[3])
}
summary(store_coreset_2)
boxplot(store_coreset_2, main = "12 Turf strain")

# task 3 group1 turf (KL3, INV, QHB1 MDB1, NCT3)
store_coreset_3 <- c()
i <- 0
for (i in 1:1000){
  sub_strains <- turf_g1_index[sample(1:length(turf_g1_index),3,replace = FALSE)]
  store_coreset_3[i] <- core_genes(sub_strains[1],sub_strains[2],sub_strains[3])
}
summary(store_coreset_3)
boxplot(store_coreset_3, main = "group 1 Turf strain")

# Task 4 group2 turf (SH7, QH1, MOR, COLB1, INDB2, SF12)
store_coreset_4 <- c()
i <- 0
for (i in 1:1000){
  sub_strains <- turf_g2_index[sample(1:length(turf_g2_index),3,replace = FALSE)]
  store_coreset_4[i] <- core_genes(sub_strains[1],sub_strains[2],sub_strains[3])
}
summary(store_coreset_4)
boxplot(store_coreset_4, main = "group 2 Turf strain")

# task 5 all turf and maize pathogens (15 strains)
store_coreset_5 <- c()
i <- 0
for (i in 1:1000){
  sub_strains <- turf_maize_index[sample(1:length(turf_maize_index),3,replace = FALSE)]
  store_coreset_5[i] <- core_genes(sub_strains[1],sub_strains[2],sub_strains[3])
}
summary(store_coreset_5)
boxplot(store_coreset_5, main = "turf and maize strains")

# task 6 all strains 15 strains + Cat98-1
store_coreset_6 <- c()
i <- 0
for (i in 1:1000){
  sub_strains <- all_index[sample(1:length(all_index),3,replace = FALSE)]
  store_coreset_6[i] <- core_genes(sub_strains[1],sub_strains[2],sub_strains[3])
}
summary(store_coreset_6)
boxplot(store_coreset_6, main = "all strains")

# task 7 group1 and maize strains
store_coreset_7 <- c()
i <- 0
for (i in 1:1000){
  sub_strains <- g1_maize_index[sample(1:length(g1_maize_index),3,replace = FALSE)]
  store_coreset_7[i] <- core_genes(sub_strains[1],sub_strains[2],sub_strains[3])
}
summary(store_coreset_7)
boxplot(store_coreset_7, main = "group1 and maize strains")
# task 8 group2 and maize strains
store_coreset_8 <- c()
i <- 0
for (i in 1:1000){
  sub_strains <- g2_maize_index[sample(1:length(g2_maize_index),3,replace = FALSE)]
  store_coreset_8[i] <- core_genes(sub_strains[1],sub_strains[2],sub_strains[3])
}
summary(store_coreset_8)
boxplot(store_coreset_8, main = "group1 and maize strains")

# store all data into a dataframe
all_strains <- cbind(store_coreset_1 = c(2330),store_coreset_3,store_coreset_4,
                     store_coreset_2, store_coreset_5,store_coreset_7,
                     store_coreset_8,store_coreset_6)
all_strains <- rename(as.data.frame(all_strains),maize_strains = store_coreset_1, 
                      turf_strains = store_coreset_2, 
                      group1_turf = store_coreset_3,group2_turf = store_coreset_4,
                      maize_turf = store_coreset_5,all_strains = store_coreset_6,
                      g1_maize = store_coreset_7, g2_maize = store_coreset_8)
boxplot(all_strains, main = "Core genome size resampling", ylab = "Number of genes",
        xlab = "group of strains")
apply(all_strains,2,summary)
