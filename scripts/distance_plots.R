##### Distance measures #####
# In this script, plots and tables are made 
# to analize different measures of word diversity 
# and complexity

##### Load useful libraries #####
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(corrplot)
library(gplots)
##### Load files
## Loading euclidean distance files
euc_dist2 <- read.csv("Data/euclidean_distance2.csv")
euc_dist3 <- read.csv("Data/euclidean_distance3.csv")
euc_dist4 <- read.csv("Data/euclidean_distance4.csv")
euc_dist5 <- read.csv("Data/euclidean_distance5.csv")

##Euc_dist objects are not properly ordered to make matrix
euc_dist2[c("File1", "File2")] <- lapply(euc_dist2[c("File1", "File2")], factor, 
                           levels=unique(unlist(euc_dist2[c("File1", "File2")]))) 
#
euc_dist3[c("File1", "File2")] <- lapply(euc_dist3[c("File1", "File2")], factor, 
                                         levels=unique(unlist(euc_dist3[c("File1", "File2")]))) 
#
euc_dist4[c("File1", "File2")] <- lapply(euc_dist4[c("File1", "File2")], factor, 
                                         levels=unique(unlist(euc_dist4[c("File1", "File2")]))) 
#
euc_dist5[c("File1", "File2")] <- lapply(euc_dist5[c("File1", "File2")], factor, 
                                         levels=unique(unlist(euc_dist5[c("File1", "File2")]))) 

#ggplot(data = euc_dist2, aes(x=File1, y=File2, fill=Euc_dist, 
                               #label= value))

matrix_euc2 <- xtabs(Euc_dist ~ File1+ File2,euc_dist2)
corrplot(matrix_euc2, is.corr = FALSE)

##
matrix_euc3 <- xtabs(Euc_dist ~ File1+ File2,euc_dist3)
corrplot(matrix_euc3, is.corr = FALSE,type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05, 
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=0.5, outline=T)
hclust(matrix_euc3, method = "complete")
##
matrix_euc4 <- xtabs(Euc_dist ~ File1+ File2,euc_dist4)
corrplot(matrix_euc4, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05, 
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
##
matrix_euc5 <- xtabs(Euc_dist ~ File1+ File2,euc_dist5)
corrplot(matrix_euc5, is.corr = FALSE)
corrplot(matrix_euc5, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05, 
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
#####Load word counts to calc euc distance in R #####
wc.2 <- as.matrix(read.csv("Data/wordcounts2.csv",row.names=1))
wc.3 <- as.matrix(read.csv("Data/wordcounts3.csv", row.names = 1))
wc.3[is.na(wc.3)] <- 0
wc.4 <- as.matrix(read.csv("Data/wordcounts4.csv", row.names = 1))
wc.4[is.na(wc.4)] <- 0
wc.5 <- as.matrix(read.csv("Data/wordcounts5.csv", row.names = 1))
wc.5[is.na(wc.5)] <- 0

euc.dist.r.2 <- dist(wc.2)
euc.dist.r.3 <- dist(wc.3)
euc.dist.r.4 <- dist(wc.4)
euc.dist.r.5 <- dist(wc.5)

#Euclidean distances calculated with r are the same as the ones in python
#Try to solve problem of different total of cells
boxplot(rowSums(wc.2))
hist(rowSums(wc.2))
rowSums(wc.3)
#Try to count cell number for each and calc the counts per 100,000
rowSums(wc.2)/10000
dosporcadadiezmil <- rowSums(wc.2/2)/10000
wc.2.pordiezmil <- round((dosporcadadiezmil* wc.2),0)
#para 3
tresporcadadiezmil <- rowSums(wc.3/3)/10000
wc.3.pordiezmil <- round((tresporcadadiezmil* wc.3),0)
#para 4
cuatroporcadadiezmil <- rowSums(wc.4/4)/10000
wc.4.pordiezmil <- round((cuatroporcadadiezmil* wc.4),0)
#para 5
cincoporcadadiezmil <- rowSums(wc.5/5)/10000
wc.5.pordiezmil <- round((cincoporcadadiezmil* wc.5),0)
#####
euc.dist.r.2.std <- dist(wc.2.pordiezmil)
euc.dist.r.3.std <- dist(wc.3.pordiezmil)
euc.dist.r.4.std <- dist(wc.4.pordiezmil)
euc.dist.r.5.std <- dist(wc.5.pordiezmil)
####
plot(hclust(euc.dist.r.5.std))

pdf("Figures/eucdistwordlength4.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.4))
plot(hclust(euc.dist.r.4.std))
dev.off()
png("Figures/eucdistwordlength4.png", height = 480, width = 480) # Para guardar en PNG
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.4))
plot(hclust(euc.dist.r.4.std))
dev.off()
#
pdf("Figures/eucdistwordlength2.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.2))
plot(hclust(euc.dist.r.2.std))
dev.off()
####
heatmap.2(as.matrix(euc.dist.r.2),key=TRUE,scale = "row", 
            margins = c(10,12), cexRow=0.5)
dev.off()
heatmap(as.matrix(euc.dist.r.3))
heatmap(as.matrix(euc.dist.r.3.std))

heatmap(as.matrix(euc.dist.r.4))
heatmap(as.matrix(euc.dist.r.5))
heatmap(as.matrix(euc.dist.r.5.std))

ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile()
##### Load files of lempelziv measurments
lempel.by.cell <- read.csv("Data/lemplzivbyfile.csv", row.names=1)
#See a exploratory boxplot
boxplot(lempel.by.cell$Value ~ lempel.by.cell$Name)
###Add species identifiers to files
#Make data frame too match species
files <- levels(as.factor(lempel.by.cell$Name))
species<- c("E. bracteata","E. lomelli","E. colligata","E. coalcomanensis",
            "E. calcarata","E. calcarata","E. finkii","E. conzattii","E. cyri",
            "E. peritropoides", "E. cymbifera","E. tehuacana","L-systemMesic",
            "L-systemXeric","E. diazlunana","E. diazlunana","E. diazlunana",
            "E. tithymaloides","E. tithymaloides","E. personata","E. personata",
            "L-system")
habit <- c("xeric","xeric","mesic","mesic",
           "mesic","mesic","mesic","mesic","xeric",
           "mesic","xeric","xeric","L-systemCSM","L-systemCSX","xeric",
           "xeric","xeric","xeric","xeric","xeric","xeric","L-system")
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
#
par(mfrow = c(2,1))
plot(lempel.by.cell$Value[lempel.by.cell$species=="L-system"], type ="l", xlab="Compression distance")
plot(lempel.by.cell$Value[lempel.by.cell$Name=="EPM5_edited_cells.txt"], type ="l", xlab="Compression distance")
dev.off()
#
lempel.by.cell <- transform(lempel.by.cell, Sequence=ave(seq_along(Name), Name, FUN=seq_along))

ggplot(lempel.by.cell, aes(x = Sequence,y = Value,color = species)) +
  geom_line()
#Add identifier of habit

match.id <- match(lempel.by.cell$Name,species.id$files)

lempel.by.cell$habit <- as.character(species.id$habit[match.id])

ggplot(lempel.by.cell, aes(x = Sequence,y = Value,color = habit)) +
  geom_line()

######## Shannon entropy #######
shannonentropy.by.cell <- read.csv("Data/shannonentropy.csv", row.names=1)

#Add species and habit
match.id <- match(shannonentropy.by.cell$Name,species.id$files)
shannonentropy.by.cell$species <- as.character(species.id$species[match.id])
shannonentropy.by.cell$habit <- as.character(species.id$habit[match.id])

#
boxplot(shannonentropy.by.cell$Value ~ shannonentropy.by.cell$Name)
boxplot(shannonentropy.by.cell$Value ~ shannonentropy.by.cell$species)
boxplot(shannonentropy.by.cell$Value ~ shannonentropy.by.cell$habit)

plot(shannonentropy.by.cell$Value[lempel.by.cell$Name=="974_edited_cells.txt"])
plot(shannonentropy.by.cell$Value[lempel.by.cell$Name=="EPM6_S2-1_edited_cells.txt"])
