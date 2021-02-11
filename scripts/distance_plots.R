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
library(wesanderson)
##### Load files
## Loading euclidean distance files
euc_dist2 <- read.csv("Data/euclidean_distance2.csv")
euc_dist3 <- read.csv("Data/euclidean_distance3.csv")
euc_dist4 <- read.csv("Data/euclidean_distance4.csv")
euc_dist5 <- read.csv("Data/euclidean_distance5.csv")
euc_dist10 <- read.csv("Data/euclidean_distance10.csv")

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
#
euc_dist10[c("File1", "File2")] <- lapply(euc_dist10[c("File1", "File2")], factor, 
                                         levels=unique(unlist(euc_dist10[c("File1", "File2")]))) 

#ggplot(data = euc_dist2, aes(x=File1, y=File2, fill=Euc_dist, 
                               #label= value))

matrix_euc2 <- xtabs(Euc_dist ~ File1+ File2,euc_dist2)
corrplot(matrix_euc2, is.corr = FALSE)

##
matrix_euc3 <- xtabs(Euc_dist ~ File1+ File2,euc_dist3)
corrplot(matrix_euc3, is.corr = FALSE,type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05, 
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=0.5, outline=T)

##
matrix_euc4 <- xtabs(Euc_dist ~ File1+ File2,euc_dist4)
corrplot(matrix_euc4, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05, 
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
##
matrix_euc5 <- xtabs(Euc_dist ~ File1+ File2,euc_dist5)
corrplot(matrix_euc5, is.corr = FALSE)
corrplot(matrix_euc5, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05, 
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
###
matrix_euc10 <- xtabs(Euc_dist ~ File1+ File2,euc_dist10)
corrplot(matrix_euc10, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05, 
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)

#####Load word counts to calc euc distance in R #####
wc.2 <- as.matrix(read.csv("Data/wordcounts2.csv",row.names=1))
wc.3 <- as.matrix(read.csv("Data/wordcounts3.csv", row.names = 1))
wc.3[is.na(wc.3)] <- 0
wc.4 <- as.matrix(read.csv("Data/wordcounts4.csv", row.names = 1))
wc.4[is.na(wc.4)] <- 0
wc.5 <- as.matrix(read.csv("Data/wordcounts5.csv", row.names = 1))
wc.5[is.na(wc.5)] <- 0
wc.10 <- as.matrix(read.csv("Data/wordcounts10.csv", row.names = 1))
wc.10[is.na(wc.10)] <- 0
#
euc.dist.r.2 <- dist(wc.2)
euc.dist.r.3 <- dist(wc.3)
euc.dist.r.4 <- dist(wc.4)
euc.dist.r.5 <- dist(wc.5)
euc.dist.r.10 <- dist(wc.10)
#Get proportion of counts
wc.2.freq <- wc.2/rowSums(wc.2)
euc.dist.r2.freq <- dist(wc.2.freq)
plot(hclust(euc.dist.r2.freq))
#
wc.3.freq <- wc.3/rowSums(wc.3)
euc.dist.r3.freq <- dist(wc.3.freq)
plot(hclust(euc.dist.r3.freq))
#
wc.4.freq <- wc.4/rowSums(wc.4)
euc.dist.r4.freq <- dist(wc.4.freq)
plot(hclust(euc.dist.r4.freq))
#
wc.5.freq <- wc.5/rowSums(wc.5)
euc.dist.r5.freq <- dist(wc.5.freq)
plot(hclust(euc.dist.r5.freq))
#
wc.10.freq <- wc.10/rowSums(wc.10)
euc.dist.r10.freq <- dist(wc.10.freq)
plot(hclust(euc.dist.r10.freq))
#Euclidean distances calculated with r are the same as the ones in python
#Try to solve problem of different total of cells
boxplot(rowSums(wc.2))
hist(rowSums(wc.2))
rowSums(wc.3)
#Try to count cell number for each and calc the counts per 20,000
rowSums(wc.2)/20000
dosporcadaveintemil <- rowSums(wc.2)/20000
wc.2.porveintemil <- round((wc.2/dosporcadaveintemil),0)
wc.2.porveintemil.freq <- wc.2.porveintemil/rowSums(wc.2.porveintemil)
euc.dist.r.2.freq.std <- dist(wc.2.porveintemil.freq)
plot(hclust(euc.dist.r.2.freq.std))
#para 3
tresporcadaveintemil <- rowSums(wc.3)/20000
wc.3.porveintemil <- round((wc.3/tresporcadaveintemil),0)
wc.3.porveintemil.freq <- wc.3.porveintemil/rowSums(wc.3.porveintemil)
euc.dist.r.3.freq.std <- dist(wc.3.porveintemil.freq)
plot(hclust(euc.dist.r.3.freq.std))

#para 4
cuatroporcadaveintemil <- rowSums(wc.4)/20000
wc.4.porveintemil <- round((wc.4/cuatroporcadaveintemil),0)
wc.4.porveintemil.freq <- wc.4.porveintemil/rowSums(wc.4.porveintemil)
euc.dist.r.4.std.pr <- dist(wc.4.porveintemil)
euc.dist.r.4.freq.std <- dist(wc.4.porveintemil.freq)
plot(hclust(euc.dist.r.4.freq.std))
plot(hclust(euc.dist.r.4.std.pr))

#para 5
cincoporcadaveintemil <- rowSums(wc.5)/20000
wc.5.porveintemil <- round((wc.5/cincoporcadaveintemil),0)
wc.5.porveintemil.freq <- wc.5.porveintemil/rowSums(wc.5.porveintemil)
euc.dist.r.5.freq.std <- dist(wc.5.porveintemil.freq)
plot(hclust(euc.dist.r.5.freq.std))

#para 10
diezporcadaveintemil <- rowSums(wc.10)/20000
wc.10.porveintemil <- round((wc.10/diezporcadaveintemil),0)
wc.10.porveintemil.freq <- wc.10.porveintemil/rowSums(wc.10.porveintemil)
euc.dist.r.10.freq.std <- dist(wc.10.porveintemil.freq)
plot(hclust(euc.dist.r.10.freq.std))

#####
euc.dist.r.2.std <- dist(wc.2.porveintemil)
euc.dist.r.3.std <- dist(wc.3.porveintemil)
euc.dist.r.4.std <- dist(wc.4.porveintemil)
euc.dist.r.5.std <- dist(wc.5.porveintemil)
euc.dist.r.10.std <- dist(wc.10.porveintemil)
####

plot(hclust(euc.dist.r.5.std))

pdf("Figures/eucdistwordlength2.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r2.freq))
plot(hclust(euc.dist.r.2.freq.std))
dev.off()
pdf("Figures/eucdistwordlength2abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.2))
plot(hclust(euc.dist.r.2.std))
dev.off()
pdf("Figures/eucdistwordlength3.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r3.freq))
plot(hclust(euc.dist.r.3.freq.std))
dev.off()
pdf("Figures/eucdistwordlength3abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.3))
plot(hclust(euc.dist.r.3.std))
dev.off()
pdf("Figures/eucdistwordlength4.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r4.freq))
plot(hclust(euc.dist.r.4.freq.std))
dev.off()
pdf("Figures/eucdistwordlength4abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.4))
plot(hclust(euc.dist.r.4.std))
dev.off()
pdf("Figures/eucdistwordlength5.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r5.freq))
plot(hclust(euc.dist.r.5.freq.std))
dev.off()
pdf("Figures/eucdistwordlength5abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.5))
plot(hclust(euc.dist.r.5.std))
dev.off()
pdf("Figures/eucdistwordlength10.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r10.freq))
plot(hclust(euc.dist.r.10.freq.std))
dev.off()
pdf("Figures/eucdistwordlength10abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.10))
plot(hclust(euc.dist.r.10.std))
dev.off()

####
heatmap.2(as.matrix(euc.dist.r.2),key=TRUE,scale = "row", 
            margins = c(10,12), cexRow=0.5)
dev.off()
pdf("Figures/heatmaplength2.pdf") # Para guardar en PDF
heatmap(as.matrix(euc.dist.r2.freq))
heatmap(as.matrix(euc.dist.r.2.freq.std))
heatmap(as.matrix(euc.dist.r.2))
heatmap(as.matrix(euc.dist.r.2.std))
dev.off()
pdf("Figures/heatmaplength3.pdf") # Para guardar en PDF
heatmap(as.matrix(euc.dist.r3.freq))
heatmap(as.matrix(euc.dist.r.3.freq.std))
heatmap(as.matrix(euc.dist.r.3))
heatmap(as.matrix(euc.dist.r.3.std))
dev.off()

pdf("Figures/heatmaplength4.pdf") # Para guardar en PDF
heatmap(as.matrix(euc.dist.r4.freq))
heatmap(as.matrix(euc.dist.r.4.freq.std))
heatmap(as.matrix(euc.dist.r.4))
heatmap(as.matrix(euc.dist.r.4.std))
dev.off()
pdf("Figures/heatmaplength5.pdf") # Para guardar en PDF
heatmap(as.matrix(euc.dist.r5.freq))
heatmap(as.matrix(euc.dist.r.5.freq.std))
heatmap(as.matrix(euc.dist.r.5))
heatmap(as.matrix(euc.dist.r.5.std))
dev.off()
pdf("Figures/heatmaplength10.pdf") # Para guardar en PDF
heatmap(as.matrix(euc.dist.r10.freq))
heatmap(as.matrix(euc.dist.r.10.freq.std))
heatmap(as.matrix(euc.dist.r.10))
heatmap(as.matrix(euc.dist.r.10.std))
dev.off()

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
match.id <- match(lempel.by.cell$Name,species.id$files)
lempel.by.cell$habit <- as.character(species.id$habit[match.id])
###boxplot with species
boxplot(lempel.by.cell$Value ~ lempel.by.cell$species)
#Make palete
colores<- c("#71982d","#7464d7","#5aba46","#d369d3","#44be7a","#a2409a","#b8b436","#696fba",
            "#d68937","#5f9ed7","#c94c34","#3fc1bf","#da477c","#67b88c","#a14761","#91b869",
            "#c987c4","#4a772f","#dc8279","#327e58","#c2a864","#846a2b")

pdf("Figures/lempelzivFile.pdf") # Para guardar en PDF
boxplot(lempel.by.cell$Value ~ lempel.by.cell$Name, col=colores,
        las=1.5, xlab="Cell-file Lempel-Ziv algorithm",ylab = "Score", cex=0.4,
        cex.names=0.5, cex.axis = 0.3, cex.lab=1.2)
dev.off()
pdf("Figures/lempelzivsp.pdf")
boxplot(lempel.by.cell$Value ~ lempel.by.cell$species, col=colores,
        las=1.5, xlab = "Cell-file Lempel-Ziv algorithm",ylab = "Score", cex=0.4,
        cex.names=0.5, cex.axis=0.6, cex.lab=1.2)
dev.off()
pdf("Figures/lempelzivhabit.pdf")
boxplot(lempel.by.cell$Value ~ lempel.by.cell$habit, col=colores,
        xlab = "Cell-file Lempel-Ziv algorithm")
dev.off()

lempel.cellhabit.lm <-lm(lempel.by.cell$Value~lempel.by.cell$habit)
summary(lempel.cellhabit.lm)
lempel.cellsp.lm <-lm(lempel.by.cell$Value~lempel.by.cell$species)
summary(lempel.cellsp.lm)
#
ggplot(lempel.by.cell, aes(x = Name, y = Value, color = species)) + geom_boxplot()
ggplot(lempel.by.cell, aes(x = species, y = Value, color = species)) + geom_boxplot()

#Now plot distribution of compression
bract<-subset(lempel.by.cell, lempel.by.cell$species == "E. bracteata")
#Make loop to plot each file
plot(lempel.by.cell$Value[lempel.by.cell$species=="E. tithymaloides"], type ="l")

#####Make palete
pal1 <- wes_palette("BottleRocket2")
pal2 <- wes_palette("Rushmore1")
pal3 <- wes_palette("Darjeeling1")
pal4 <- wes_palette("FantasticFox1")
pal5 <- wes_palette("IsleofDogs1")
my_palette <- c(pal1,pal2,pal3,pal4,pal5)
#
pdf("Figures/lempelzivbyfile.pdf")
par(mfrow = c(2,1))
z=1
for(i in files) {
  temp <- subset(lempel.by.cell, lempel.by.cell$Name == i)
  plot(temp$Value, type ="l", xlab = "Cell file position", ylab="Lempel-ziv value", 
       main= i, col= my_palette[z])
  z = z +1
}
dev.off()
##
plot(lempel.by.cell$Value[lempel.by.cell$species=="L-system"], type ="b", xlab="Compression distance",
     col=my_palette[1], lwd=2)
z=1
for(i in files) {
  plot(lempel.by.cell$Value[lempel.by.cell$Name==i], type ="b", xlab="Compression distance",
       col=my_palette[1], lwd=2,  main= i)
  lines(lempel.by.cell$Value[lempel.by.cell$Name=="probLsystem_edited_cells.txt"], type ="b", xlab="Compression distance",
        col = my_palette[z], lwd = 2)
  z=z+1
}
#
pdf("Figures/lempelziv_bycell.pdf")
plot(lempel.by.cell$Value[lempel.by.cell$Name=="974_edited_cells.txt"], type ="b", xlab="Compression distance",
     col=my_palette[1], lwd=2,  main= "Lempel-ziv")
lines(lempel.by.cell$Value[lempel.by.cell$Name=="probLsystem_edited_cells.txt"], type ="b", xlab="Compression distance",
      col = my_palette[2], lwd = 2)
lines(lempel.by.cell$Value[lempel.by.cell$Name=="EPM10_edited_cells.txt"], type ="b", xlab="Compression distance",
      col = my_palette[3], lwd = 2)
legend("bottom", legend = c("E. peritropoides","L-system","E. diazlunana"),
       col =levels(as.factor(my_palette[1:3])) , pch = 16, inset =-0.15,xpd=TRUE,horiz=TRUE)
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
ggplot(lempel.by.cell, aes(x = Sequence,y = Value,color = habit)) +
  geom_line()

######## Shannon entropy #######
shannonentropy.by.cell <- read.csv("Data/shannonentropy.csv", row.names=1)

#Add species and habit
match.id <- match(shannonentropy.by.cell$Name,species.id$files)
shannonentropy.by.cell$species <- as.character(species.id$species[match.id])
shannonentropy.by.cell$habit <- as.character(species.id$habit[match.id])

#
pdf("Figures/shannonbyfile.pdf") # Para guardar en PDF
boxplot(shannonentropy.by.cell$Value ~ shannonentropy.by.cell$Name, col=colores,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        cex.names=0.5, cex.axis=0.3, cex.lab=1.2)
dev.off()
pdf("Figures/shanonbysp.pdf")
boxplot(shannonentropy.by.cell$Value ~ shannonentropy.by.cell$species, col=colores,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        names=c( expression(italic("E. bracteata")),expression(italic("E. calcarata")) ,
                 expression(italic("E. coalcomanensis")),expression(italic("E. colligata")), 
                 expression(italic("E. conzattii")),expression(italic("E. cymbifera")),
                 expression(italic("E. cyri")), expression(italic("E. diazlunana")),
                 expression(italic("E. finkii")),expression(italic("E. lomelli")),
                 expression(italic("E. peritropoides")),  expression(italic("E. personata")),
                 expression(italic("E. tehuacana")),expression(italic("E. tithymaloides")),
                 expression("L-system"),expression("L-systemM"), expression("LsystemX")),
                  cex.names=0.5, cex.axis=0.8, cex.lab=1.2)
dev.off()
pdf("Figures/shanonbyhabit.pdf") # Para guardar en PDF
boxplot(shannonentropy.by.cell$Value ~ shannonentropy.by.cell$habit,col=colores,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4)
dev.off()
shannondiff.lm <- lm(shannonentropy.by.cell$Value~shannonentropy.by.cell$species)
summary(shannondiff.lm)

pdf("Figures/entropybyfile.pdf")
par(mfrow = c(2,1))
z=1
for(i in files) {
  temp <- subset(shannonentropy.by.cell, shannonentropy.by.cell$Name == i)
  plot(temp$Value, type ="l", xlab = "Cell file position", ylab="Shannon-Entropy value", 
       main= i, col= my_palette[z], lwd = 2)
  z = z +1
}
dev.off()
######
#Small subset
shannon_subset<-subset(shannonentropy.by.cell, Name == "probLsystem_edited_cells.txt" | Name=="974_edited_cells.txt" |
                               Name == "EPM10_edited_cells.txt")
shannon_subset <- transform(shannon_subset, Sequence=ave(seq_along(species), 
                                                         species, FUN=seq_along))

#Shannon
pdf("Figures/shannon_bycell.pdf")
ggplot(shannon_subset, aes(x = Sequence,y = Value,color = species)) +
  geom_line()
dev.off()

#######Checkout files with ray cells########
shannonentropy.by.cell.with.ray <- read.csv("Data/shannonentropywithR.csv", row.names=1)

#Add species and habit
match.id <- match(shannonentropy.by.cell.with.ray$Name,species.id$files)
shannonentropy.by.cell.with.ray$species <- as.character(species.id$species[match.id])
shannonentropy.by.cell.with.ray$habit <- as.character(species.id$habit[match.id])

boxplot(shannonentropy.by.cell.with.ray$Value ~ shannonentropy.by.cell.with.ray$Name, 
        col=colores,las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        cex.names=0.5, cex.axis=0.3, cex.lab=1.2)

plot(shannonentropy.by.cell.with.ray$Value[shannonentropy.by.cell.with.ray$Name=="974_edited_cells.txt"], 
     type ="l", xlab="Compression distance")
lines(shannonentropy.by.cell$Value[shannonentropy.by.cell$Name=="974_edited_cells.txt"],
     type ="l", col =my_palette[2])
