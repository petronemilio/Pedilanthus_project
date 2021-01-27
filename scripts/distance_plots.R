##### Distance measures #####
# In this script, plots and tables are made 
# to analize different measures of word diversity 
# and complexity

##### Load useful libraries #####
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
##### Load files
## Loading euclidean distance files
euc_dist2 <- read.csv("Data/euclidean_distance2.csv")
euc_dist3 <- read.csv("Data/euclidean_distance3.csv")
euc_dist4 <- read.csv("Data/euclidean_distance4.csv")
euc_dist5 <- read.csv("Data/euclidean_distance5.csv")

##
matrix(euc_dist$Euc_dist, dimnames = )
list(euc_dist$File1)
matrix_euc <- xtabs(Euc_dist ~ File1+ File2,euc_dist)

xtabs(value~var1+var2, df)
matrix_euc

euc_dist2_acast <-acast(euc_dist2, File1~File2, value.var="Euc_dist")


euc_dist %>% pivot_wider(names_from = File2, values_from = Euc_dist)

hclust(matrix_euc,method = 'complete')
as.factor(euc_dist$File1)


##### Load files of lempelziv measurments
lempel.by.cell <- read.csv("Data/lemplzivbyfile.csv", row.names=1)
#See a exploratory boxplot
boxplot(lempel.by.cell$Value ~ lempel.by.cell$Name)
###Add species identifiers to files
#Make data frame too match species
files <- levels(as.factor(lempel.by.cell$Name))
species<- c("E. bracteata","E. lomelli","E. colligata","E. coalcomanensis",
            "E. calcarata","E. calcarata","E. finkii","E. conzattii","E. cyri",
           "E. peritropoides", "E. cymbifera","E. tehuacana","E. diazlunana",
           "E. diazlunana","E. diazlunana","E. tithymaloides","E. tithymaloides",
           "E. personata","E. personata")
habit <- c("xeric","xeric","mesic","mesic",
           "mesic","mesic","mesic","mesic","xeric",
           "mesic","xeric","xeric","xeric",
           "xeric","xeric","xeric","xeric","xeric","xeric")

species.id <- as.data.frame(cbind(files,species,habit))

match.id <- match(lempel.by.cell$Name,species.id$files)

lempel.by.cell$species <- as.character(species.id$species[match.id])
###boxplot with species
boxplot(lempel.by.cell$Value ~ lempel.by.cell$species)
#Make palete
boxplot(lempel.by.cell$Value ~ lempel.by.cell$Name, boxfill= lempel.by.cell$species)
#
ggplot(lempel.by.cell, aes(x = Name, y = Value, color = species)) + geom_boxplot()
ggplot(lempel.by.cell, aes(x = species, y = Value, color = species)) + geom_boxplot()

#Now plot distribution of compression
bract<-subset(lempel.by.cell, lempel.by.cell$species == "E. bracteata")
#Make loop to plot each file
plot(lempel.by.cell$Value[lempel.by.cell$species=="E. tithymaloides"], type ="l")

par(mfrow = c(2,2))
for(i in files) {
  temp <- subset(lempel.by.cell, lempel.by.cell$Name == i)
  plot(temp$Value, type ="l", xlab="Compression distance", )
}
dev.off()

lempel.by.cell <- transform(lempel.by.cell, Sequence=ave(seq_along(Name), Name, FUN=seq_along))

ggplot(lempel.by.cell, aes(x = Sequence,y = Value,color = species)) +
  geom_line()
#Add identifier of habit

match.id <- match(lempel.by.cell$Name,species.id$files)

lempel.by.cell$habit <- as.character(species.id$habit[match.id])

ggplot(lempel.by.cell, aes(x = Sequence,y = Value,color = habit)) +
  geom_line()
