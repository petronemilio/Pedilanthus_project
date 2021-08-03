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
library(pheatmap)
##### Load files
## Loading euclidean distance files
euc_dist2 <- read.csv("Data/euclidean_distance2.csv")
euc_dist3 <- read.csv("Data/euclidean_distance3.csv")
euc_dist4 <- read.csv("Data/euclidean_distance4.csv")
euc_dist5 <- read.csv("Data/euclidean_distance5.csv")
euc_dist6 <- read.csv("Data/euclidean_distance6.csv")
euc_dist7 <- read.csv("Data/euclidean_distance7.csv")
euc_dist8 <- read.csv("Data/euclidean_distanceR8.csv")
euc_dist9 <- read.csv("Data/euclidean_distance9.csv")
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
euc_dist6[c("File1", "File2")] <- lapply(euc_dist6[c("File1", "File2")], factor, 
                                         levels=unique(unlist(euc_dist6[c("File1", "File2")]))) 
#
euc_dist7[c("File1", "File2")] <- lapply(euc_dist7[c("File1", "File2")], factor, 
                                         levels=unique(unlist(euc_dist7[c("File1", "File2")]))) 
#
euc_dist8[c("File1", "File2")] <- lapply(euc_dist8[c("File1", "File2")], factor, 
                                         levels=unique(unlist(euc_dist8[c("File1", "File2")]))) 
#
euc_dist9[c("File1", "File2")] <- lapply(euc_dist9[c("File1", "File2")], factor, 
                                         levels=unique(unlist(euc_dist9[c("File1", "File2")]))) 
#
euc_dist10[c("File1", "File2")] <- lapply(euc_dist10[c("File1", "File2")], factor, 
                                         levels=unique(unlist(euc_dist10[c("File1", "File2")]))) 
#####Transform to matrix format#####
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
corrplot(matrix_euc5, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05, 
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
###
matrix_euc6 <- xtabs(Euc_dist ~ File1+ File2,euc_dist6)
corrplot(matrix_euc6, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05, 
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
#
matrix_euc7 <- xtabs(Euc_dist ~ File1+ File2,euc_dist7)
corrplot(matrix_euc7, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05, 
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
#
matrix_euc8 <- xtabs(Euc_dist ~ File1+ File2,euc_dist8)
corrplot(matrix_euc8, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05, 
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
#
matrix_euc9 <- xtabs(Euc_dist ~ File1+ File2,euc_dist9)
corrplot(matrix_euc9, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05, 
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
#
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
wc.6 <- as.matrix(read.csv("Data/wordcounts6.csv", row.names = 1))
wc.6[is.na(wc.6)] <- 0
wc.7 <- as.matrix(read.csv("Data/wordcounts7.csv", row.names = 1))
wc.7[is.na(wc.7)] <- 0
wc.8 <- as.matrix(read.csv("Data/wordcounts8.csv", row.names = 1))
wc.8[is.na(wc.8)] <- 0
wc.9 <- as.matrix(read.csv("Data/wordcounts9.csv", row.names = 1))
wc.9[is.na(wc.9)] <- 0
wc.10 <- as.matrix(read.csv("Data/wordcounts10.csv", row.names = 1))
wc.10[is.na(wc.10)] <- 0
#
euc.dist.r.2 <- dist(wc.2)
euc.dist.r.3 <- dist(wc.3)
euc.dist.r.4 <- dist(wc.4)
euc.dist.r.5 <- dist(wc.5)
euc.dist.r.6 <- dist(wc.6)
euc.dist.r.7 <- dist(wc.7)
euc.dist.r.8 <- dist(wc.8)
euc.dist.r.9 <- dist(wc.9)
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
wc.6.freq <- wc.6/rowSums(wc.6)
euc.dist.r6.freq <- dist(wc.6.freq)
plot(hclust(euc.dist.r6.freq))
#
wc.7.freq <- wc.7/rowSums(wc.7)
euc.dist.r7.freq <- dist(wc.7.freq)
plot(hclust(euc.dist.r7.freq))
#
wc.8.freq <- wc.8/rowSums(wc.8)
euc.dist.r8.freq <- dist(wc.8.freq)
plot(hclust(euc.dist.r8.freq))
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
rowSums(wc.2)
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
#para 6
seisporcadaveintemil <- rowSums(wc.6)/20000
wc.6.porveintemil <- round((wc.6/seisporcadaveintemil),0)
wc.6.porveintemil.freq <- wc.6.porveintemil/rowSums(wc.6.porveintemil)
euc.dist.r.6.freq.std <- dist(wc.6.porveintemil.freq)
plot(hclust(euc.dist.r.6.freq.std))
#Para 7
sieteporcadaveintemil <- rowSums(wc.7)/20000
wc.7.porveintemil <- round((wc.7/sieteporcadaveintemil),0)
wc.7.porveintemil.freq <- wc.7.porveintemil/rowSums(wc.7.porveintemil)
euc.dist.r.7.freq.std <- dist(wc.7.porveintemil.freq)
plot(hclust(euc.dist.r.7.freq.std))
#para 8
ochoporcadaveintemil <- rowSums(wc.8)/20000
wc.8.porveintemil <- round((wc.8/ochoporcadaveintemil),0)
#Para 9
nueveporcadaveintemil <- rowSums(wc.9)/20000
wc.9.porveintemil <- round((wc.9/nueveporcadaveintemil),0)
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
euc.dist.r.6.std <- dist(wc.6.porveintemil)
euc.dist.r.7.std <- dist(wc.7.porveintemil)
euc.dist.r.8.std <- dist(wc.8.porveintemil)
euc.dist.r.9.std <- dist(wc.9.porveintemil)
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
pdf("Figures/eucdistwordlength6.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r6.freq))
plot(hclust(euc.dist.r.6.freq.std))
dev.off()
pdf("Figures/eucdistwordlength6abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.6))
plot(hclust(euc.dist.r.6.std))
dev.off()
pdf("Figures/eucdistwordlength7.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r7.freq))
plot(hclust(euc.dist.r.7.freq.std))
dev.off()
pdf("Figures/eucdistwordlength7abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.7))
plot(hclust(euc.dist.r.7.std))
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

pheatmap(as.matrix(euc.dist.r.2.std), display_numbers = F)

pdf("Figures/heatmaplength2.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r2.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.2.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.2), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.2.std), display_numbers = F)
dev.off()
pdf("Figures/heatmaplength3.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r3.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.3.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.3), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.3.std), display_numbers = F)
dev.off()

pdf("Figures/heatmaplength4.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r4.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.4.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.4), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.4.std), display_numbers = F)
dev.off()
pdf("Figures/heatmaplength5.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r5.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.5.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.5), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.5.std), display_numbers = F)
dev.off()
pdf("Figures/heatmaplength6.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r6.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.6.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.6), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.6.std), display_numbers = F)
dev.off()
pdf("Figures/heatmaplength7.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r7.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.7.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.7), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.7.std), display_numbers = F)
dev.off()
pdf("Figures/heatmaplength10.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r10.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.10.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.10), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.10.std), display_numbers = F)
dev.off()
##Falta hacer los análisis de words incluyendo el sistema radial

##### Load files of lempelziv measurments #####
lempel.by.cell <- read.csv("Data/lemplzivbyfile.csv", row.names=1)
#See a exploratory boxplot
boxplot(lempel.by.cell$Value ~ lempel.by.cell$Name)
###Add species identifiers to files
#Make data frame too match species
files <- levels(as.factor(lempel.by.cell$Name))
species<- c("E. bracteata","E. lomelli","E. colligata","E. coalcomanensis",
            "E. calcarata","E. calcarata","E. finkii","E. tithymaloides","E. conzattii",
            "E. cyri","E. peritropoides", "E. cymbifera","E. tehuacana","L-systemMesic",
            "L-systemXeric","E. diazlunana","E. diazlunana","E. diazlunana","E. lomelli",
            "E. cyri", "E. cymbifera", "E. tithymaloides","E. tithymaloides","E. personata",
            "E. personata","probL-system","probetaL-system","RayL-system")
habit <- c("xeric","xeric","mesic","mesic",
           "mesic","mesic","mesic","xeric","mesic",
           "xeric","mesic","xeric","xeric","L-systemCSM",
           "L-systemCSX","xeric","xeric","xeric","xeric",
           "xeric","xeric","xeric","xeric","xeric",
           "xeric","probL-system","probetaL-system","RayL-system")

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
            "#c987c4","#4a772f","#dc8279","#327e58","#c2a864","#846a2b","#71982d")

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
plot(lempel.cellhabit.lm)
library(emmeans)
emmeans(lempel.cellhabit.lm, ~ habit)
lempel.cellsp.lm <-lm(lempel.by.cell$Value~lempel.by.cell$species)
summary(lempel.cellsp.lm)
##Remove rlsystem, and context sensitive lsystems
lempel.by.cell_subset <-subset(lempel.by.cell, lempel.by.cell$species != "RayL-system")
lempel.by.cell_subset <-subset(lempel.by.cell_subset, lempel.by.cell_subset$species != "L-systemMesic")
lempel.by.cell_subset <-subset(lempel.by.cell_subset, lempel.by.cell_subset$species != "L-systemXeric")
lempel.by.cell_subset$species <- factor(lempel.by.cell_subset$species)
#
lempel.cellsp.lm <-lm(lempel.by.cell_subset$Value ~ factor(lempel.by.cell_subset$species))
emm.lempel <- emmeans(lempel.cellsp.lm, specs = pairwise ~ species,adjust="tukey")
emm.lempel$contrasts
multcomp::cld(emm.lempel$emmeans, alpha = 0.10, Letters=LETTERS)
#
pdf("Figures/lempelzivsp_reduced.pdf")
boxplot(lempel.by.cell_subset$Value ~ lempel.by.cell_subset$species, 
        las=1.5, xlab = "Cell-file Lempel-Ziv algorithm",ylab = "Score", cex=0.4,
        cex.names=0.5, cex.axis=0.6, cex.lab=1.2)
dev.off()
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
pal6 <- wes_palette("Cavalcanti1")
my_palette <- c(pal1,pal2,pal3,pal4,pal5,pal6)
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
plot(lempel.by.cell$Value[lempel.by.cell$species=="probL-system"], type ="b", xlab="Compression distance",
     col=my_palette[1], lwd=2)
plot(lempel.by.cell$Value[lempel.by.cell$species=="probetaL-system"], type ="b", xlab="Compression distance",
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
plot(lempel.by.cell$Value[lempel.by.cell$species=="probL-system"], type ="l", xlab="Compression distance")
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
                 expression("L-systemM"),expression("L-systemX"), expression("probBetaLsystem"),
                 expression("probL-system"),expression("RayL-system")),
                  cex.names=0.5, cex.axis=0.8, cex.lab=1.2)
dev.off()
pdf("Figures/shanonbyhabit.pdf") # Para guardar en PDF
boxplot(shannonentropy.by.cell$Value ~ shannonentropy.by.cell$habit,col=colores,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4)
dev.off()


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
#######
##Remove rlsystem, and context sensitive lsystems
shannonentropy.by.cell_subset <-subset(shannonentropy.by.cell, shannonentropy.by.cell$species != "RayL-system")
shannonentropy.by.cell_subset <-subset(shannonentropy.by.cell_subset, 
                                       shannonentropy.by.cell_subset$species != "L-systemMesic")
shannonentropy.by.cell_subset <-subset(shannonentropy.by.cell_subset, 
                                       shannonentropy.by.cell_subset$species != "L-systemXeric")
shannonentropy.by.cell_subset$species <- factor(shannonentropy.by.cell_subset$species)
#
shannon.cellsp.lm <-lm(shannonentropy.by.cell_subset$Value ~ 
                      factor(shannonentropy.by.cell_subset$species))
emm.shannon <- emmeans(shannon.cellsp.lm, specs = pairwise ~ species,adjust="tukey")
emm.shannon$contrasts
multcomp::cld(emm.shannon$emmeans, alpha = 0.10, Letters=LETTERS)
#
pdf("Figures/shanonbysp_subset.pdf")
boxplot(shannonentropy.by.cell_subset$Value ~ shannonentropy.by.cell_subset$species,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        names=c( expression(italic("E. bracteata")),expression(italic("E. calcarata")) ,
                 expression(italic("E. coalcomanensis")),expression(italic("E. colligata")), 
                 expression(italic("E. conzattii")),expression(italic("E. cymbifera")),
                 expression(italic("E. cyri")), expression(italic("E. diazlunana")),
                 expression(italic("E. finkii")),expression(italic("E. lomelli")),
                 expression(italic("E. peritropoides")),  expression(italic("E. personata")),
                 expression(italic("E. tehuacana")),expression(italic("E. tithymaloides")),
                expression("probBetaLsystem"),expression("probL-system")),
        cex.names=0.5, cex.axis=0.8, cex.lab=1.2)
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
######Check out Shanon by window #####
shannonentropy.by.window <- read.csv("Data/shannonentropy_window.csv", row.names=1)
library(dplyr)
shannonentropy.by.window <- shannonentropy.by.window %>% 
                            group_by(Sample,Lineage) %>% mutate(id = row_number())


window947 <- subset(shannonentropy.by.window, Sample == "974_edited_cells.txt")
window883 <- subset(shannonentropy.by.window, Sample == "883_edited_cells.txt")
window883 <- subset(shannonentropy.by.window, Sample == "883_edited_cells.txt")
prbeta <-subset(shannonentropy.by.window, Sample == "probLsystembeta_edited_cells.txt")

pdf("Figures/shannon_window947.pdf")
plot(window947$id[window947$Lineage==2],window947$Value[window947$Lineage==2],type = "l",lwd = 3)
points(window947$id[window947$Lineage==3],window947$Value[window947$Lineage==3], col="green",type = "l",lwd = 3)
points(window947$id[window947$Lineage==4],window947$Value[window947$Lineage==4], col="red",type = "l",lwd = 3)
points(window947$id[window947$Lineage==5],window947$Value[window947$Lineage==5], col=pal1[1],type = "l",lwd = 3)
points(window947$id[window947$Lineage==6],window947$Value[window947$Lineage==6], col=pal1[2],type = "l",lwd = 3)
#points(window947$id[window947$Lineage==7],window947$Value[window947$Lineage==7], col=pal1[3],pch=19)
#points(window947$id[window947$Lineage==8],window947$Value[window947$Lineage==8], col=pal2[1],pch=19)
#points(window947$id[window947$Lineage==9],window947$Value[window947$Lineage==9], col=pal2[2],pch=19)
#points(window947$id[window947$Lineage==10],window947$Value[window947$Lineage==10], col=pal2[3],pch=19)
dev.off()

pdf("Figures/shannon_window883.pdf")
plot(window883$id[window883$Lineage==3],window883$Value[window883$Lineage==3],type = "l",lwd = 3)
points(window883$id[window883$Lineage==2],window883$Value[window883$Lineage==2], col="green",type = "l",lwd = 3)
points(window883$id[window883$Lineage==6],window883$Value[window883$Lineage==6], col="red",type = "l",lwd = 3)
points(window883$id[window883$Lineage==4],window883$Value[window883$Lineage==4], col=pal1[1],type = "l",lwd = 3)
points(window883$id[window883$Lineage==7],window883$Value[window883$Lineage==7], col=pal1[2],type = "l",lwd = 3)
dev.off()

pdf("Figures/probetaS.pdf")
plot(prbeta$id[prbeta$Lineage==6],prbeta$Value[prbeta$Lineage==6],type = "l",lwd = 3)
points(prbeta$id[prbeta$Lineage==2],prbeta$Value[prbeta$Lineage==2], col="green",type = "l",lwd = 3)
points(prbeta$id[prbeta$Lineage==1],prbeta$Value[prbeta$Lineage==1], col="red",type = "l",lwd = 3)
points(prbeta$id[prbeta$Lineage==4],prbeta$Value[prbeta$Lineage==4], col=pal1[1],type = "l",lwd = 3)
points(prbeta$id[prbeta$Lineage==7],prbeta$Value[prbeta$Lineage==7], col=pal1[2],type = "l",lwd = 3)
dev.off()

summary(cell_lengths)

ggplot(window947, aes(x= id, y = Value, color=Lineage)) +
  geom_point()

ggplot(lempel.by.cell, aes(x = Name, y = Value, color = species)) + geom_boxplot()
