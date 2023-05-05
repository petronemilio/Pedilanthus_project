#Para poner tablas en formato wide agregar bibliotecas 
library(tidyr)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(ggstance)
library(vioplot)
library(plot3D)
library(scatterplot3d)
library(data.table)
library(fitdistrplus)
library(stringr)
library(scales)
##sET WORKING DIRECTORY
getwd()
###

#####Save palettes for graphs
pal1 <- wes_palette("BottleRocket2")
pal2 <- wes_palette("Rushmore1")
pal3 <- wes_palette("Darjeeling1")
pal4 <- wes_palette("FantasticFox1")
my_palette <- c(pal1,pal2,pal3,pal4)
############################################
#Graficar la distribución de las longitudes de las series. Including rays
cell_lengths <- read.csv("Data/cell_lengths_notConverge.csv")
files <- levels(as.factor(cell_lengths$Sample))
species<- c("E. peritropoides","E. peritropoides","E. bracteata","E. lomelli","E. colligata","E. coalcomanensis",
            "E. coalcomanensis","E. calcarata","E. calcarata","E. finkii","E. calcarata","E. tithymaloides",
            "E. conzattii","E. cyri","E. peritropoides", "E. cymbifera","E. tehuacana","L-systemMesic",
            "L-systemXeric","E. diazlunana","E. diazlunana","E. diazlunana","E. lomelli",
            "E. cyri", "E. cymbifera", "E. tithymaloides","E. tithymaloides","E. personata",
            "E. personata","probL-system","probetaL-system","RayL-system")
habit <- c("mesic","mesic","xeric","xeric","mesic","mesic",
           "mesic","mesic","mesic","mesic","mesic","xeric",
           "mesic","xeric","mesic","xeric","xeric","L-systemCSM",
           "L-systemCSX","xeric","xeric","xeric","xeric",
           "xeric","xeric","xeric","xeric","xeric",
           "xeric","probL-system","probetaL-system","RayL-system")
species.id <- as.data.frame(cbind(files,species,habit))
insert <- "_NotConverge"
#species.id$filesconverge <- sub("(?<=\\w)(?=\\.)", insert, species.id$files, perl=TRUE)
#write.csv(species.id,"meta/speciesID.csv")
#####
match.id <- match(cell_lengths$Sample,species.id$files)
cell_lengths$species <- as.character(species.id$species[match.id])
#
cell_lengths$habit <- as.character(species.id$habit[match.id])
agg <- aggregate(Number.of.cells ~ species, cell_lengths, function(x){
  qq <- quantile(x, probs = c(1, 3)/4)
  iqr <- diff(qq)
  lo <- qq[1] - 1.5*iqr
  hi <- qq[2] + 1.5*iqr
  c(Mean = mean(x), IQR = unname(iqr), lower = lo, high = hi)
}) 
agg
#####
p  <- ggplot(cell_lengths, aes(Number.of.cells, colour=species, fill=species))
p  <- p + geom_density(alpha=0.2)
p
cell_lengths_subset<-subset(cell_lengths, habit != "RayL-system")
cell_lengths_subset<-subset(cell_lengths_subset, species != "L-systemMesic")
cell_lengths_subset<-subset(cell_lengths_subset, species != "L-systemXeric")
q  <- ggplot(cell_lengths_subset, aes(Number.of.cells, colour=species, fill=species))
q  <- q + geom_density(alpha=0.2)
pdf("Figures/lengths_density.pdf")
q
dev.off()
pal<-c(wes_palette("Cavalcanti1"),wes_palette("GrandBudapest1"),wes_palette("GrandBudapest2"),
       wes_palette("FantasticFox1"))
par(mar=c(7,5,1,1))
new_pal<-c("#cb57aa","#6db744","#8760cf","#c2af45","#7d7fc5","#dd8d4a","#45b0cf","#cf483c",
           "#5cc08c","#c56179","#3b824e","#976530","#76853a")

boxplot(cell_lengths$Number.of.cells~species, cell_lengths,horizontal=TRUE,
        col=new_pal)
boxplot(cell_lengths$Number.of.cells~cell_lengths$species,
        col=new_pal)
#Remove cell lengths minor to 3 cells
cell_lengths <- subset(cell_lengths, cell_lengths$Number.of.cells > 3)
cell_lengths$species <- with(cell_lengths, reorder(species,Number.of.cells,mean))

boxplot(Number.of.cells~species, cell_lengths,horizontal = TRUE, las=1.5,
        col=c(pal), xlab="Longitud de las filas",ylab = NULL, cex=0.4,
        names=c(expression(italic("E. lomelli")),expression(italic("E. cymbifera")),
                expression(italic("E. tithymaloides")), expression(italic("E. diazlunana")),
                expression(italic("E. personata")), expression(italic("E. finkii")),
                expression(italic("E. bracteata")),expression(italic("E. coalcomanensis")),
                expression(italic("E. peritropoides")),expression(italic("E. colligata")),
                expression(italic("E. cyri")),expression(italic("E. tehuacana")),
                expression(italic("E. calcarata")),expression(italic("E. conzattii"))),
        cex.names=1.1, cex.axis=0.8, cex.lab=1.2)
##
boxplot(cell_lengths$Number.of.cells~cell_lengths$habit,
        col=new_pal)
##
summary(cell_lengths$Number.of.cells)
cellfilemean<-aggregate(cell_lengths$Number.of.cells, list(cell_lengths$species), FUN=mean)
cellfilemean[with(cellfilemean, order(x)),]

numbercells<-aggregate(cell_lengths$Number.of.cells, list(cell_lengths$Sample), FUN=mean)
##
hist(cell_lengths$Number.of.cells[cell_lengths$habit=="xeric"])
hist(cell_lengths$Number.of.cells[cell_lengths$habit=="mesic"])
hist(cell_lengths$Number.of.cells[cell_lengths$habit=="L-systemCSM"])
hist(cell_lengths$Number.of.cells[cell_lengths$habit=="probetaL-system"])

q  <- ggplot(cell_lengths, aes(Number.of.cells, colour=habit, fill=habit))
q  <- q + geom_density(alpha=0.2)
q

descdist(cell_lengths$Number.of.cells[cell_lengths$habit=="xeric"], discrete = FALSE)
descdist(cell_lengths$Number.of.cells[cell_lengths$habit=="mesic"], discrete = FALSE)

ggplot(cell_lengths, aes(x =habit,y=Number.of.cells, fill = habit)) +
  geom_violin()+ coord_flip()+
  scale_y_continuous(breaks = c(0,50,100,150,200,250,300,350,400,500,600))+
  scale_fill_manual(values=c(pal))+
  theme(legend.position = "none") 
###
descdist(cell_lengths$Number.of.cells)
# Draw the plot
calca <- subset(cell_lengths, species=="E. calcarata", select=c(species,Number.of.cells))
tithy <- subset(cell_lengths, species=="E. tithymaloides", select=c(species, Number.of.cells))
cyri <- subset(cell_lengths, species=="E. cyri", select=c(species, Number.of.cells))
brac <- subset(cell_lengths, species=="E. bracteata", select=c(species, Number.of.cells))
finkii<- subset(cell_lengths,species=="E. finkii", select=c(species, Number.of.cells))
colli<- subset(cell_lengths,species=="E. colligata", select=c(species, Number.of.cells))
diazlunana<- subset(cell_lengths,species=="E. diazlunana", select=c(species, Number.of.cells))
coalco<- subset(cell_lengths,species=="E. coalcomanensis", select=c(species, Number.of.cells))
lomelli<- subset(cell_lengths,species=="E. lomelli", select=c(species, Number.of.cells))
conza<- subset(cell_lengths,species=="E. conzattii", select=c(species, Number.of.cells))
peri<- subset(cell_lengths,species=="E. peritropoides", select=c(species, Number.of.cells))
tehua<- subset(cell_lengths,species=="E. tehuacana", select=c(species, Number.of.cells))
cymbi<- subset(cell_lengths,species=="E. cymbifera", select=c(species, Number.of.cells))
perso<- subset(cell_lengths,species=="E. personata", select=c(species, Number.of.cells))
#
descdist(calca$Number.of.cells)
descdist(tithy$Number.of.cells)
descdist(cyri$Number.of.cells)
descdist(brac$Number.of.cells)
descdist(finkii$Number.of.cells)
descdist(colli$Number.of.cells)
descdist(diazlunana$Number.of.cells)
descdist(coalco$Number.of.cells)
descdist(lomelli$Number.of.cells)
descdist(conza$Number.of.cells)
descdist(peri$Number.of.cells)
descdist(colli$Number.of.cells)
descdist(peri$Number.of.cells)
descdist(tehua$Number.of.cells)
descdist(cymbi$Number.of.cells)
descdist(perso$Number.of.cells)

pdf("Figures/cell_lengths.pdf")
vioplot(brac$Number.of.cells,lomelli$Number.of.cells,tithy$Number.of.cells,
        cyri$Number.of.cells,cymbi$Number.of.cells, tehua$Number.of.cells,
        diazlunana$Number.of.cells,perso$Number.of.cells,
        colli$Number.of.cells,coalco$Number.of.cells,calca$Number.of.cells,  
        finkii$Number.of.cells,conza$Number.of.cells,peri$Number.of.cells,
        names=c("Br","Lo","Ti","Cyr","Cym","Teh","Di","Pe","Col","Coa","Cal","Fin",
                "Con","Peri"),
        col=c("#899DA4","#899DA4","#899DA4","#899DA4","#899DA4","#899DA4","#899DA4","#899DA4",
              "#C93312","#C93312","#C93312","#C93312","#C93312","#C93312"))
dev.off()
#remove RayL-system, L-systemMesic and L-systemXeric
cell_lengths_subset <-subset(cell_lengths, cell_lengths$species != "RayL-system")
cell_lengths_subset <-subset(cell_lengths_subset, cell_lengths_subset$species != "L-systemMesic")
cell_lengths_subset <-subset(cell_lengths_subset, cell_lengths_subset$species != "L-systemXeric")
cell_lengths_subset$species <- factor(cell_lengths_subset$species)
lm.cellength <- lm(log10(cell_lengths_subset$Number.of.cells) ~ factor(cell_lengths_subset$species))
summary(lm.cellength)
#
lm.cellength.stdres <- rstandard(lm.cellength)
qqnorm(lm.cellength.stdres, 
         ylab="Standardized Residuals", 
         xlab="Normal Scores", 
             main="Cell lengths from fusiform initials") 
qqline(lm.cellength.stdres)
qqplot(cell_lengths$Number.of.cells)
library(emmeans)
emm1 <-emmeans(lm.cellength,specs = pairwise ~ species,adjust="tukey")
emm1$contrasts
multcomp::cld(emm1$emmeans, alpha = 0.10, Letters=LETTERS)
pdf("Figures/cell_lengths_bygroup.pdf")
boxplot(log10(cell_lengths_subset$Number.of.cells) ~ cell_lengths_subset$species, notch=T,
        col=c("#FAD510","#FAD510","#E2D200","#F2AD00","#F2AD00","#F2AD00",
             "#F98400", "#F2AD00","#F2300F","#CB2314", "#FF0000","#35274A","#35274A",
             "#35274A","#35274A","#354823"))
dev.off()
aggregate(Number.of.cells~ species, data=cell_lengths, mean)
aggregate(Number.of.cells~ species, data=cell_lengths, sd)
#Making by growth form
lm.cellength.habit <- lm(log10(cell_lengths_subset$Number.of.cells) ~ factor(cell_lengths_subset$habit))
summary(lm.cellength.habit)
emm2 <-emmeans(lm.cellength.habit,specs = pairwise ~ habit,adjust="tukey")
emm2$contrasts
multcomp::cld(emm2$emmeans, alpha = 0.10, Letters=LETTERS)
pwpm(emm2$emmeans)
pdf("Figures/cell_lengths_byhabit.pdf")
boxplot(log10(cell_lengths_subset$Number.of.cells) ~ cell_lengths_subset$habit, notch=T,
        col=c("#FAD510","#CB2314","#FF0000","#35274A"))
dev.off()

############################################################
#Graficar los tipos celulares
tipos_celulares <- read.csv("Data/wordcountsR1.csv")
colnames(tipos_celulares)[1] = "Sample"
match.id <- match(tipos_celulares$Sample,species.id$files)
tipos_celulares$species <- as.character(species.id$species[match.id])
#replace nas with 0s
tipos_celulares[is.na(tipos_celulares)] <- 0

tipos_celulares$totalcells <-rowSums(tipos_celulares[,c(2:5)])
tipos_celulares_freq <- tipos_celulares[,c(2:5,7)] / tipos_celulares[,7] 
tipos_celulares_freq$sample <- tipos_celulares$Sample
tipos_celulares_freq$species <- as.character(species.id$species[match.id])

##### Obtain a table of cells per sample and number of cells
numberofcells <- table(cell_lengths$Sample)
rownames(numberofcells)
match.id <- match(tipos_celulares$Sample,rownames(numberofcells))
tipos_celulares$numberofFiles <- as.character(numberofcells[match.id])
write.table(tipos_celulares, "Data/cell_Lengths_Celltypes.csv", row.names = FALSE)
#Obtain total of cells codified without the L-systems
tipos_celulares_withoutLS <-subset(tipos_celulares, tipos_celulares$Sample!="contextmesicLsystem_edited_cells.txt"&
         tipos_celulares$Sample != "contextxericLsystem_edited_cells.txt" & 
         tipos_celulares$Sample != "contextxericLsystem_edited_cells.txt" & 
         tipos_celulares$Sample != "probLsystem_edited_cells.txt" & 
         tipos_celulares$Sample != "probLsystembeta_edited_cells.txt" & 
         tipos_celulares$Sample != "Ray_Lsystem_edited_cells.txt")
sum(tipos_celulares_withoutLS$totalcells)
sum(tipos_celulares_withoutLS$F)/sum(tipos_celulares_withoutLS$totalcells)
sum(tipos_celulares_withoutLS$R)/sum(tipos_celulares_withoutLS$totalcells)
sum(tipos_celulares_withoutLS$P)/sum(tipos_celulares_withoutLS$totalcells)
sum(tipos_celulares_withoutLS$V)/sum(tipos_celulares_withoutLS$totalcells)

sum(as.numeric(tipos_celulares_withoutLS$numberofFiles))
tipos_celulares_freq$numberofFiles <- as.character(numberofcells[match.id])
tipos_celulares_freq$numberofFiles <- as.character(numberofcells[match.id])

#
#Reshape dat to plot
tipos_celulares <- reshape(data=tipos_celulares, idvar="x", varying = c("V","F","P","R"),
        times=c("V","F","P","R"),v.name=c("Count"),direction="long")
#
ggplot(tipos_celulares) + 
  geom_bar(mapping=aes(x = Sample, y = Count, fill= time), stat = "identity",
           position = "fill") +coord_flip()
pdf("Figures/count_freqs.pdf")
ggplot(tipos_celulares) + 
  geom_bar(mapping=aes(x = species , y = Count, fill= time), stat = "identity",
           position = "fill") +coord_flip()
dev.off()
##
tipos_celulares_freq <- reshape(data=tipos_celulares_freq, idvar="x", varying = c("V","F","P","R"),
                           times=c("V","F","P","R"),v.name=c("Count"),direction="long")
#
#### Make plots removing cell rays #####
tipos_fusiformes <- subset(tipos_celulares, time != "R")
tipos_fusiformes$countfreq <-tipos_fusiformes$Count/tipos_fusiformes$totalcells
#tipos_fusiformes <-order(tipos_fusiformes$time)
tipos_fusiformes<-tipos_fusiformes[with(tipos_fusiformes, order(time,countfreq)), ]
ggplot(tipos_fusiformes) + 
  geom_bar(mapping=aes(x =species, y = countfreq, fill= time), stat = "identity",
           position = "fill") + coord_flip()
ggplot(tipos_fusiformes,aes(x =species, y = countfreq, fill= time))+ 
  geom_col(position = position_fill(reverse = TRUE)) + coord_flip()

pdf("Figures/count_freqs_fusiform.pdf")
ggplot(tipos_fusiformes) + 
  geom_bar(mapping=aes(x =species, y = Count, fill= time), 
           stat = "identity", position = "fill") + coord_flip()
dev.off()
##
tipos_celulares_freq <- tipos_celulares[,c(2:5,7)] / tipos_celulares[,7] 
tipos_celulares_freq$sample <- tipos_celulares$Sample
tipos_celulares_freq$species <- as.character(species.id$species[match.id])

#####Make statistical test #####
tipos_celulares <- read.csv("Data/wordcountsR1.csv", row.names = 1)
tipos_celulares <-subset(tipos_celulares, select= -R)
contingency_table <-addmargins(as.matrix(tipos_celulares))
#
tipos_celulares_chi<-chisq.test(contingency_table)
#####
tipos_celulares_chi$expected
#Gtest
library(DescTools)
tipos_celulares.gi <- GTest(contingency_table)
tipos_celulares_chi$stdres
#
library(vcd) # useful package for graphics, beyond scope of this book

png("Figures/contingency_table.png")
mosaic(as.matrix(tipos_celulares), gp=shading_Friendly, residuals=tipos_celulares_chi$stdres,
       residuals_type="Std\nresiduals", labeling=labeling_residuals)
dev.off()
tipos_celulares_chi
#library(RVAideMemoire)
#chisq.multcomp(contingency_table, p.method = "none")

####
lvls <- names(sort(tapply(tipos_fusiformes$time == "Fibras", tipos_fusiformes$especies_por_muestra, mean)))
lvls
ggplot(tipos_fusiformes)+ geom_bar(aes(factor(especies_por_muestra, levels = lvls),
           stat="identity")) 
  #geom_bar(position = "fill") + scale_y_continuous(labels = count)

#ggplot(data=dat1, aes(x=time, y=total_bill, fill=sex)) +
#  geom_bar(stat="identity")

##########Graficar los indices de homogeneídad ##########
#
homogenity_index <- read.csv("Data/homogentiy_index.csv")
#Take files as factor  
files <- levels(as.factor(homogenity_index$X0))
#Add species to the homogenity index
##Make data frame for files
species.id <- as.data.frame(cbind(files,species,habit))
match.id <- match(homogenity_index$X0,species.id$files)
#
homogenity_index$sp <- as.character(species.id$species[match.id])
#    
homogenity_index$habit <- as.character(species.id$habit[match.id])
####
#Calc the length proportion standaridze per sample...
prueba<-homogenity_index %>% group_by(X0) %>%
  top_n(1, length) 
prueba <- prueba[!duplicated(prueba[ , c("X0", "length")]), ]  # Delete rows
matcher <- match(homogenity_index$X0, prueba$X0)
homogenity_index$maxlength<- prueba$length[matcher]
homogenity_index$lengthstd <- (homogenity_index$length)/(homogenity_index$maxlength)

ggplot(homogenity_index, aes(x=support, y=conductivity, colour=sp))+
  geom_point() #+ scale_color_manual(values=my_palette)
ggplot(homogenity_index, aes(x=conductivity, y=support, colour=habit))+
  geom_point() #+ scale_color_manual(values=my_palette)

ggplot(homogenity_index, aes(x=lengthstd, y=conductivity, colour=habit))+
  geom_point() #+ scale_color_manual(values=my_palette)

pdf("Figures/storage_conductivity.pdf") # Para guardar en PDF
ggplot(homogenity_index, aes(x=storage, y=conductivity, colour=habit))+
  geom_point() + scale_color_manual(values=c(pal3,pal2))
dev.off()
#
pdf("Figures/storage_support.pdf") # Para guardar en PDF
ggplot(homogenity_index, aes(x=storage, y=support, colour=habit))+
  geom_point() + scale_color_manual(values=c(pal3,pal4))
dev.off()
#
species_index <- subset(homogenity_index, habit == "xeric" | habit == "mesic" | habit=="probetaL-system")
pdf("Figures/support_conductivity.pdf") # Para guardar en PDF
ggplot(species_index, aes(x=conductivity, y=support, colour=sp))+
  geom_point() +  scale_color_manual(values=c(my_palette))
dev.off()

pdf("Figures/support_conductivity1.pdf")
ggplot(species_index, aes(x=conductivity, y=support, colour=habit))+
  geom_point() +  scale_color_manual(values=c(my_palette))
dev.off()
######
pdf("Figures/3dscatter.pdf")
scatter3D(homogenity_index$storage,homogenity_index$conductivity,
          homogenity_index$support, 
          col.var = as.integer(as.factor(homogenity_index$habit)), 
          col = pal3, 
          pch = 4, ticktype = "detailed",
          colkey = list(at = c(2,3,4,5,6)), side = 1,addlines = TRUE, length = 0.5, width = 0.5,
          labels = c("mesic","xeric","L-system","L-systemCSX","L-systemCSM"),
          #id=list(method = "mahal", n = length(homogenity_index$habit), 
           #       labels = homogenity_index$habit),
          phi = 0, bty ="g",
          main = "Homogenity index", xlab = "Storage index",
          ylab ="Conductivity index", zlab = "Support index")
dev.off()

pal3 <-pal3[as.numeric(as.factor(homogenity_index$habit))]
pdf("Figures/morphospace1.pdf")
scatterplot3d(homogenity_index$storage,homogenity_index$conductivity,
          homogenity_index$support, angle=55, pch=16, color=pal3,
          main="3D Scatter Plot",xlab = "Storage index",ylab = "Conductivity index",
          zlab = "Support index") 
legend("bottom", legend = c("mesic","xeric","L-system","L-systemCSX","L-systemCSM"),
       col =levels(as.factor(pal3)) , pch = 16, inset =-0.15,xpd=TRUE,horiz=TRUE)
dev.off()
levels(as.factor(homogenity_index$habit))

pdf("Figures/morphospace2.pdf")
scatterplot3d(homogenity_index$storage,homogenity_index$support,
              homogenity_index$conductivity,
               angle=55, pch=16, color=pal3,
              main="3D Scatter Plot",xlab = "Storage index",ylab = "Support index",
              zlab = "Conductivity index") 
legend("bottom", legend = levels(as.factor(homogenity_index$habit)),
       col =levels(as.factor(pal3)) , pch = 16, inset =-0.15,xpd=TRUE,horiz=TRUE)
dev.off()

scatter3D(homogenity_index$storage,homogenity_index$conductivity,
          homogenity_index$support, phi = 0, bty ="g",
          id=list(method = "mahal", n = length(as.factor(homogenity_index$habit)),
                  labels = as.factor(homogenity_index$habit)))
#
plot3d(geometry[,1],geometry[,2],geometry[,3])
text3d(geometry[,1],geometry[,2],geometry[,3],rownames(geometry))
points3d(geometry[,1],geometry[,2],geometry[,3], size = 5)

#scale_color_manual(values=c("#c80966","#84de66","#8369e1","#d89816",
 #                             "#0152a1","#005e18","#abb2ff"))+
  #theme_bw()

homogenity_index <- subset(homogenity_index, homogenity_index$habit != "L-systemCSM" & 
         homogenity_index$habit != "L-systemCSX" & homogenity_index$habit != "RayL-system" ) 
###Plot using length
pdf("Figures/support_filelength.pdf")
plot(NULL, xlim=c(0,700), ylim=c(-0.7,1), ylab="y label", xlab="x lablel")
points(homogenity_index$support[homogenity_index$habit=="probL-system"] ~
      homogenity_index$length[homogenity_index$habit=="probL-system"], 
       col=alpha("#CB2314",0.9),pch=19)
points(homogenity_index$support[homogenity_index$habit=="probetaL-system"] ~
         homogenity_index$length[homogenity_index$habit=="probetaL-system"], 
       col=alpha("#CB2314",0.9),pch=19)
points(homogenity_index$support[homogenity_index$habit=="mesic"] ~
       homogenity_index$length[homogenity_index$habit=="mesic"], 
       col=alpha("#FAD510", 0.6),pch=19)
points(homogenity_index$support[homogenity_index$habit=="xeric"] ~
         homogenity_index$length[homogenity_index$habit=="xeric"], 
       col=alpha("#273046",0.4), pch=19)
dev.off()
##
pdf("Figures/sotrage_filelength_separate.pdf")
par(mfrow=c(3,1))
plot(NULL, xlim=c(0,700), ylim=c(-0.7,1), ylab="y label", xlab="x lablel")
points(homogenity_index$support[homogenity_index$habit=="mesic"] ~
         homogenity_index$length[homogenity_index$habit=="mesic"], col="#FAD510",pch=19)
plot(NULL, xlim=c(0,700), ylim=c(-0.6,1), ylab="y label", xlab="x lablel")
points(homogenity_index$support[homogenity_index$habit=="xeric"] ~
         homogenity_index$length[homogenity_index$habit=="xeric"], col="#273046",pch=19)
plot(NULL, xlim=c(0,700), ylim=c(-0.6,1), ylab="y label", xlab="x lablel")
points(homogenity_index$support[homogenity_index$habit=="probL-system"] ~
         homogenity_index$length[homogenity_index$habit=="probL-system"], col="#CB2314",pch=19)
points(homogenity_index$support[homogenity_index$habit=="probetaL-system"] ~
         homogenity_index$length[homogenity_index$habit=="probetaL-system"], col="#CB2314",pch=19)
dev.off()
#
pdf("Figures/conductivity_filelength.pdf")
plot(NULL, xlim=c(0,700), ylim=c(-0.3,1), ylab="y label", xlab="x lablel")
points(homogenity_index$conductivity[homogenity_index$habit=="probL-system"] ~
         homogenity_index$length[homogenity_index$habit=="probL-system"], 
       col=alpha("#CB2314",0.9),pch=19)
points(homogenity_index$conductivity[homogenity_index$habit=="probetaL-system"] ~
         homogenity_index$length[homogenity_index$habit=="probetaL-system"],
       col=alpha("#CB2314",0.9),pch=19)
points(homogenity_index$conductivity[homogenity_index$habit=="mesic"] ~
         homogenity_index$length[homogenity_index$habit=="mesic"],  
       col=alpha("#FAD510",0.6),pch=19)
points(homogenity_index$conductivity[homogenity_index$habit=="xeric"] ~
         homogenity_index$length[homogenity_index$habit=="xeric"], 
       col=alpha("#273046",0.4),pch=19)
dev.off()
#
pdf("Figures/conductivity_filelength_separate.pdf")
par(mfrow=c(3,1))
plot(NULL, xlim=c(0,700), ylim=c(-0.3,1), ylab="y label", xlab="x lablel")
points(homogenity_index$conductivity[homogenity_index$habit=="mesic"] ~
         homogenity_index$length[homogenity_index$habit=="mesic"],  col="#FAD510",pch=19)
plot(NULL, xlim=c(0,700), ylim=c(-0.3,1), ylab="y label", xlab="x lablel")
points(homogenity_index$conductivity[homogenity_index$habit=="xeric"] ~
         homogenity_index$length[homogenity_index$habit=="xeric"], col="#273046",pch=19)
plot(NULL, xlim=c(0,700), ylim=c(-0.3,1), ylab="y label", xlab="x lablel")
points(homogenity_index$conductivity[homogenity_index$habit=="probL-system"] ~
         homogenity_index$length[homogenity_index$habit=="probL-system"], col="#CB2314",pch=19)
points(homogenity_index$conductivity[homogenity_index$habit=="probetaL-system"] ~
         homogenity_index$length[homogenity_index$habit=="probetaL-system"], col="#CB2314",pch=19)
dev.off()
#
#
homogenity_index_subset <- subset(homogenity_index,homogenity_index$sp != "L-systemXeric" & 
         homogenity_index$sp != "L-systemMesic" & homogenity_index$sp != "RayL-system")

table(homogenity_index_subset$conductivity> 0.0)
homogenity_index_justpedilanthus <- subset(homogenity_index_subset,homogenity_index$sp != "probL-system" & 
                                           homogenity_index_subset$sp !="probetaL-system")
table(homogenity_index_justpedilanthus$conductivity> 0.0)
table(homogenity_index_justpedilanthus$conductivity== 1)
temporal<-subset(homogenity_index_justpedilanthus,homogenity_index_justpedilanthus$conductivity == 1)
table(temporal$conductivityfreq==0)
temporal <- subset(temporal, temporal$conductivityfreq==0)
table(temporal$habit)
#
table(homogenity_index_justpedilanthus$conductivity <=0)
temporal<-subset(homogenity_index_justpedilanthus ,homogenity_index_justpedilanthus$conductivity <=0)
table(temporal$habit)
temporal<-subset(homogenity_index_justpedilanthus,homogenity_index_justpedilanthus$support < 0)
temporal<-subset(homogenity_index_justpedilanthus,homogenity_index_justpedilanthus$conductivity > 0.5 &
                 homogenity_index_justpedilanthus$conductivity < 0.85)
table(homogenity_index_justpedilanthus$support <1)
table(temporal$habit=="mesic")
temporal<-subset(temporal,temporal$habit=="mesic")
table(temporal$support < 0)
temporal<-subset(homogenity_index_justpedilanthus,homogenity_index_justpedilanthus$conductivity < 0.4)
table(temporal$habit=="xeric")
temporal<-subset(temporal,temporal$habit=="xeric")
table(temporal$support>0)
temporal<-subset(homogenity_index_justpedilanthus,homogenity_index_justpedilanthus$habit=="xeric")
temporal<-subset(temporal, temporal$conductivity > 0.5)
temporal1<-subset(temporal, temporal$support >= 0.0)
1279/1804
table(temporal1$support)
#Calc some values of the l-systems
homogeneity_index_lsystems<- subset(homogenity_index_subset, homogenity_index_subset$habit != "mesic" &
                               homogenity_index_subset$habit != "xeric")

table(homogeneity_index_lsystems$conductivity > 0.7)
table(homogeneity_index_lsystems$conductivity > 0.1 & 
        homogeneity_index_lsystems$conductivity < 0.7)
table(homogeneity_index_lsystems$support > 0.0 & 
        homogeneity_index_lsystems$support < 0.5)

#
pal4 <-my_palette[as.numeric(as.factor(homogenity_index_subset$habit))]
pdf("Figures/morphospace_length.pdf")
scatterplot3d(homogenity_index_subset$length,homogenity_index_subset$conductivity,
              homogenity_index_subset$support,
              angle=35, pch=16, color=pal4,
              main="3D Scatter Plot",xlab = "Cell length",ylab = "Support index",
              zlab = "Conductivity index") 
legend("bottom", legend = levels(as.factor(homogenity_index_subset$habit)),
       col =c("#273046", "#CB2314", "#FAD510") , pch = 16, inset =-0.15,xpd=TRUE,horiz=TRUE)
dev.off()
scatterplot3d(homogenity_index$support,homogenity_index$conductivity,
              homogenity_index$length,
              angle=40, pch=16, color=pal3,
              main="3D Scatter Plot",xlab = "Cell length",ylab = "Support index",
              zlab = "Conductivity index") 
legend("bottom", legend = levels(as.factor(homogenity_index$habit)),
       col =levels(as.factor(pal3)) , pch = 16, inset =-0.15,xpd=TRUE,horiz=TRUE)
#
scatter3D(homogenity_index$length,homogenity_index$conductivity,
          homogenity_index$support, 
          col.var = as.integer(as.factor(homogenity_index$habit)), 
          col = pal3, 
          pch = 4, ticktype = "detailed",
          colkey = list(at = c(2,3,4,5,6)), side = 1,addlines = TRUE, length = 0.5, width = 0.5,
          labels = c("mesic","xeric","L-system","L-systemCSX","L-systemCSM"),
          id=list(method = "mahal", n = length(homogenity_index$habit), 
                 labels = homogenity_index$habit),
          phi = 0, bty ="g",
          main = "Homogenity index", xlab = "Cell length
          ",
          ylab ="Conductivity index", zlab = "Support index")
#####Load file with info about cell lengthHIco
