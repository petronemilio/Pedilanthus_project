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
#Graficar la distribución de las longitudes de las series.
cell_lengths <- read.csv("Data/cell_lengths.csv")
files <- levels(as.factor(cell_lengths$Sample))
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
#####
match.id <- match(cell_lengths$Sample,species.id$files)
cell_lengths$species <- as.character(species.id$species[match.id])
#
cell_lengths$habit <- as.character(species.id$habit[match.id])
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
cellfilemean<-aggregate(cell_lengths$Number.of.cells, list(cell_lengths$species), FUN=mean)

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
lm.cellength <- lm(cell_lengths_subset$Number.of.cells ~ factor(cell_lengths_subset$species))
summary(lm.cellength)
library(emmeans)
emm1 <-emmeans(lm.cellength,specs = pairwise ~ species,adjust="tukey")
emm1$contrasts
multcomp::cld(emm1$emmeans, alpha = 0.10, Letters=LETTERS)
pdf("Figures/cell_lengths_bygroup.pdf")
boxplot(cell_lengths_subset$Number.of.cells ~ cell_lengths_subset$species, notch=T,
        col=c("#FAD510","#FAD510","#E2D200","#E2D200","#F2AD00","#F98400",
             "#F2300F","#CB2314", "#FF0000","#00A08A","#35274A","#35274A",
             "#35274A","#35274A","#354823","#1E1E1E"))
dev.off()
aggregate(Number.of.cells~ species, data=cell_lengths, mean)
aggregate(Number.of.cells~ species, data=cell_lengths, sd)

#cell_lengths  %>% group_by(species) %>%
 # summarise(sd = sd(Number.of.cells))

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
write.table(tipos_celulares, "Data/cell_Lengths_Celltypes.csv")
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
#tipos_fusiformes <-order(tipos_fusiformes$time)
ggplot(tipos_fusiformes) + 
  geom_bar(mapping=aes(x = Sample, y = Count, fill= time), stat = "identity",
           position = "fill") + coord_flip()
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
library(RVAideMemoire)
chisq.multcomp(contingency_table, p.method = "none")

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
  
ggplot(homogenity_index, aes(x=support, y=conductivity, colour=sp))+
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

#
sp.974 <- subset(homogenity_index, homogenity_index$X0 == '974')  
sp.883 <- subset(homogenity_index, homogenity_index$X0 == '883')  
sp.probLsystem <- subset(homogenity_index, homogenity_index$X0 == 'probLsystem') 
sp.contextmesicLsystem <- subset(homogenity_index, 
                                 homogenity_index$X0 == 'contextmesicLsystem')


ggplot(sp.974, aes(x=storage, y=conductivity, colour=X0))+
  geom_point()

ggplot(sp.883, aes(x=storage, y=conductivity, colour=X0))+
  geom_point()

ggplot(sp.probLsystem, aes(x=storage, y=conductivity, colour=X0))+
  geom_point()


ggplot(sp.contextmesicLsystem, aes(x=storage, y=conductivity, colour=X0))+
  geom_point()
#scale_color_manual(values=c("#c80966","#84de66","#8369e1","#d89816",
 #                             "#0152a1","#005e18","#abb2ff"))+
  #theme_bw()
scale_

