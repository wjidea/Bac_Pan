############### Install packages and load in data ############################

# install and load necessary packages
doInstall <- FALSE  # Change to FALSE if you don't want packages installed.
toInstall <- c("ggplot2","magrittr", "dplyr", "vegan", "RColorBrewer")
if(doInstall){install.packages(toInstall, repos = "http://cran.r-project.org")}
lapply(toInstall, library, character.only = TRUE)


# Load data 
data_genes <- read.csv("../2015_05_core_pan_genome_1.csv",header = TRUE, 
                       na.strings = c("","NA"),stringsAsFactors = FALSE)
# head(data_genes)

# data cleanning
data_sub <- select(data_genes, c(Gene, No..isolates, AA38,AA78.5,Aa99.2,COLB1,Cat98.1,
                   INDB2, INV, KL3, MDB1, MOR, NCT3, QH1, QHB1, SF12, SH7, Sa2))
data_sub <- rename(data_sub, Isolate_number = No..isolates, AA78 = AA78.5, Aa99 = Aa99.2, Cat98 = Cat98.1 )
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

