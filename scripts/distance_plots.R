## Script to check euclidean distance between 
setwd("Doctorado/Pedilanthus_project/scripts/")
## Loading euclidean distance files
euc_dist2 <- read.csv("../Data/euclidean_distance2.csv")
euc_dist3 <- read.csv("../Data/euclidean_distance3.csv")
euc_dist4 <- read.csv("../Data/euclidean_distance4.csv")
euc_dist5 <- read.csv("../Data/euclidean_distance5.csv")

##
matrix(euc_dist$Euc_dist, dimnames = )
list(euc_dist$File1)
matrix_euc <- xtabs(Euc_dist ~ File1+ File2,euc_dist)

xtabs(value~var1+var2, df)
matrix_euc

library(tidyr)
library(dplyr)
library(reshape2)
euc_dist2_acast <-acast(euc_dist2, File1~File2, value.var="Euc_dist")


euc_dist %>% pivot_wider(names_from = File2, values_from = Euc_dist)

hclust(matrix_euc,method = 'complete')
as.factor(euc_dist$File1)
