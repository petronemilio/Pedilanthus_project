##### Distance measures #####
# In this script, plots and tables are made
# to analize different measures of word diversity
# and complexity
pal1<-c("#a5c533", "#9658d8", "#59d261", "#9c3db2", "#75bb38", "#5f77f3", "#bebe35", "#3659cb",
                 "#40af43", "#ca37a3", "#39c477", "#ee68d2", "#4f901e", "#cf6fe4", "#8fa837", "#8163d6",
                 "#e6ac35", "#614bae", "#7dbe68", "#a33791", "#4dc38d", "#d32f80", "#49852d", "#ac84e6",
                 "#bca53a", "#507fe3", "#e47f2e", "#4a9de3", "#e55632", "#3ec7da", "#b82e25", "#59ceb8",
                 "#de3953", "#308949", "#e861ac", "#3d733d", "#ea93e1", "#606b16", "#82479d", "#90be7e",
                 "#6b4fa0", "#cb8e35", "#345fa9", "#b3511f", "#58b8e3", "#c33361", "#54a27a", "#eb658b",
                 "#287e63", "#c972bd", "#7c9451", "#5259a9", "#928432", "#8897ec", "#7e5f16", "#8775be",
                 "#bcb66e", "#a15195", "#1a6447", "#b04a79", "#33a29e", "#c24f4e", "#2e73a9", "#eb8665",
                 "#5e99c8", "#af6938", "#788acb", "#8e4b27", "#c0a8e8", "#5b672f", "#775696", "#d69e6c",
                 "#576196", "#8f6f3d", "#b17db7", "#9e4d51", "#e18db3", "#92465f", "#dc8484", "#8b4c76")
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
## Loading species id file
species.id <- read.csv("meta/speciesID.csv")
insert <- "_NotConverge"
species.id$filesconverge <- sub("(?<=\\w)(?=\\.)", insert, species.id$files, perl=TRUE)

##Maybe include radial files
##### Load files of lempelziv measurments #####
lempel.by.cell <- read.csv("Data/lemplzivbyfile.csv", row.names=1)
#See a exploratory boxplot
boxplot(lempel.by.cell$Value ~ lempel.by.cell$Name)
###Add species identifiers to files
match.id <- match(lempel.by.cell$Name,species.id$files)
lempel.by.cell$species <- as.character(species.id$species[match.id])
match.id <- match(lempel.by.cell$Name,species.id$files)
lempel.by.cell$habit <- as.character(species.id$habit[match.id])
###boxplot with species
boxplot(lempel.by.cell$Value ~ lempel.by.cell$species)
#Remove Lsystems:
lempel.by.cell <- subset(lempel.by.cell, lempel.by.cell$Name != "contextmesicLsystem_edited_cells_NotConverge.txt" &
                           lempel.by.cell$Name != "contextxericLsystem_edited_cells_NotConverge.txt" &
                           lempel.by.cell$Name != "Ray_Lsystem_edited_cells_NotConverge.txt")
hist(lempel.by.cell$Value)
hist(log10(lempel.by.cell$Value))
#
lempel.aov <- aov(log10(lempel.by.cell$Value)~ lempel.by.cell$habit)#Factor is habit
lempel.species.aov <- aov(log10(lempel.by.cell$Value)~ lempel.by.cell$species)#factor is sp.
summary(lempel.aov) #habit
plot(lempel.aov)
TukeyHSD(lempel.aov)
summary(lempel.species.aov)
plot(lempel.species.aov)
lempel.species.tukey <- TukeyHSD(lempel.species.aov)
plot(lempel.species.tukey, las = 1)
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
boxplot(lempel.by.cell$Value ~ lempel.by.cell$habit, 
        col=c("#CB2314","#FAD510","#273046"),
        xlab = "Cell-file Lempel-Ziv algorithm")
dev.off()
pdf("Figures/lempelzivhabit_log10.pdf")
boxplot((log10(lempel.by.cell$Value)) ~ lempel.by.cell$habit, 
        col=c("#CB2314","#FAD510","#273046"),
        xlab = "Cell-file Lempel-Ziv algorithm")
dev.off()
#
#lempel.cellhabit.lm <-lm(lempel.by.cell$Value~lempel.by.cell$habit)
#summary(lempel.cellhabit.lm)
#plot(lempel.cellhabit.lm)
library(emmeans)
library(multcompView)
emm1 <- emmeans(lempel.species.aov, specs = pairwise ~ species,adjust="tukey")
emm1$contrasts
multcomp::cld(emm1$emmeans, alpha = 0.10, Letters=LETTERS)
#
#lempel.cellsp.lm <-lm(lempel.by.cell$Value~lempel.by.cell$species)
#summary(lempel.cellsp.lm)
##
emm.habit.lempel <- emmeans(lempel.aov, specs = pairwise ~ habit, adjust="tukey")
emm.habit.lempel$contrasts
multcomp::cld(emm.habit.lempel$emmeans, alpha = 0.10, Letters=LETTERS)
#
pdf("Figures/lempelzivsp_reduced.pdf")
boxplot(lempel.by.cell$Value ~ lempel.by.cell$species,
        col=c("#273046", "#FAD510","#FAD510","#FAD510","#FAD510","#273046","#273046","#273046",
                       "#FAD510","#273046","#FAD510","#273046","#273046","#273046","#CB2314","#CB2314"),
                       las=1.5, xlab = "Cell-file Lempel-Ziv algorithm",ylab = "Score", cex=0.4,
        cex.names=0.5, cex.axis=0.6, cex.lab=1.2)
dev.off()
pdf("Figures/lempelzivsp_reduced_log.pdf")
boxplot(log10(lempel.by.cell$Value) ~ lempel.by.cell$species,
        col=c("#273046", "#FAD510","#FAD510","#FAD510","#FAD510","#273046","#273046","#273046",
                       "#FAD510","#273046","#FAD510","#273046","#273046","#273046","#CB2314","#CB2314"),
                       las=1.5, xlab = "Cell-file Lempel-Ziv algorithm",ylab = "Score", cex=0.4,
        cex.names=0.5, cex.axis=0.6, cex.lab=1.2)
dev.off()

ggplot(lempel.by.cell, aes(x = Name, y = Value, color = species)) + geom_boxplot()
ggplot(lempel.by.cell, aes(x = species, y = Value, color = species)) + geom_boxplot()
#
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
#
z=1
for(i in files) {
   plot(lempel.by.cell$Value[lempel.by.cell$Name==i], type ="b", xlab="Compression distance",
        col=my_palette[1], lwd=2,  main= i)
   lines(lempel.by.cell$Value[lempel.by.cell$Name=="probLsystem_edited_cells.txt"], type ="b", xlab="Compression distance",
         col = my_palette[z], lwd = 2)
   z=z+1
}
##
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

################################################
################# Shannon entropy #######################################################
########################################################
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
#######
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
#########Remove rlsystem, and context sensitive lsystems
shannonentropy.by.cell_subset <-subset(shannonentropy.by.cell, shannonentropy.by.cell$Name != "contextmesicLsystem_edited_cells_NotConverge.txt" &
                                         shannonentropy.by.cell$Name != "contextxericLsystem_edited_cells_NotConverge.txt" &
                                         shannonentropy.by.cell$Name != "Ray_Lsystem_edited_cells_NotConverge.txt" )
#
hist(shannonentropy.by.cell_subset$Value)
#Check mean and var
plot(tapply(shannonentropy.by.cell_subset$Value, shannonentropy.by.cell_subset$species, var),
       tapply(shannonentropy.by.cell_subset$Value,shannonentropy.by.cell_subset$species, mean))
#check by habit
plot(tapply(shannonentropy.by.cell_subset$Value, shannonentropy.by.cell_subset$habit, var),
     tapply(shannonentropy.by.cell_subset$Value,shannonentropy.by.cell_subset$habit, mean))
#
shannon_morethanzerolessonefive <- subset(shannonentropy.by.cell_subset, shannonentropy.by.cell_subset$Value > 0.1 & 
                                            shannonentropy.by.cell_subset$Value < 1.5)
hist(shannon_morethanzerolessonefive$Value)
shannon.cellsp.aov <- aov(shannonentropy.by.cell_subset$Value ~
                            factor(shannonentropy.by.cell_subset$species))
plot(shannon.cellsp.aov)
shannon.cellsp.filtered.aov <- aov(shannon_morethanzerolessonefive$Value ~
                                     shannon_morethanzerolessonefive$species)
plot(shannon.cellsp.filtered.aov)
shannon.cellsp.filtered.aov
shannon.cellsp.kruskal <- kruskal.test(shannonentropy.by.cell_subset$Value ~
                                                  shannonentropy.by.cell_subset$species)
kruskal.test(shannonentropy.by.cell_subset$Value ~ shannonentropy.by.cell_subset$species)
shannon.cellhabit.kruskal <- kruskal.test(shannonentropy.by.cell_subset$Value ~
                                                     shannonentropy.by.cell_subset$habit)
kruskal.test(shannonentropy.by.cell_subset$Value ~shannonentropy.by.cell_subset$habit)

pairwise.wilcox.test.habit <- pairwise.wilcox.test(shannonentropy.by.cell_subset$Value,shannonentropy.by.cell_subset$habit,
                     p.adjust.method = "BH")
pairwise.wilcox.test.species <- pairwise.wilcox.test(shannonentropy.by.cell_subset$Value,shannonentropy.by.cell_subset$species,
                     p.adjust.method = "BH")
#
library(rcompanion)
#Compare habit
shannon.wilcox.habit.pvalue <- pairwise.wilcox.test.habit$p.value
shannon.wilcox.full.habit.pvalue <- fullPTable(shannon.wilcox.habit.pvalue)
multcompLetters(shannon.wilcox.full.habit.pvalue)
#Species
shannon.wilcox.pvalue <- pairwise.wilcox.test.species$p.value
shannon.wilcox.full.pvalue <- fullPTable(shannon.wilcox.pvalue)
multcompLetters(shannon.wilcox.full.pvalue)
#
aggregate(lempel.by.cell$Value, list(lempel.by.cell$habit), FUN=mean) 
aggregate(shannonentropy.by.cell_subset$Value, list(shannonentropy.by.cell_subset$habit), FUN=median) 

####
summary(shannon.cellsp.kruskal)
summary(shannon.cellhabit.kruskal)
#Check if dunn test should follow.
library(npmc) #previosly run devtools::install_github("cran/npmc")
#dat <- data.frame(var = shannonentropy.by.cell_subset$Value, 
 #                 class = shannonentropy.by.cell_subset$species)
#summary(npmc(dat), type = "Steel")
#Make paiwise wilcox

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

#                
pdf("Figures/shanonbysp_subset_filtered.pdf")
boxplot(shannonentropy.by.cell_subset$Value ~ shannonentropy.by.cell_subset$species,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        col=c("#273046", "#FAD510","#FAD510","#FAD510","#FAD510","#273046","#273046","#273046",
                       "#FAD510","#273046","#FAD510","#273046","#273046","#273046","#CB2314","#CB2314"),
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
#by category
pdf("Figures/shanonbyhabit.pdf") # Para guardar en PDF
boxplot(shannonentropy.by.cell_subset$Value ~ shannonentropy.by.cell_subset$habit,
        col=c("#CB2314","#FAD510","#273046"),
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4)
dev.off()    

