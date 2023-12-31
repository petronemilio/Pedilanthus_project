#Add libraries needed for the script!
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
##cHECK WORKING DIRECTORY
getwd()
########Save palettes for graphs
pal1 <- wes_palette("BottleRocket2")
pal2 <- wes_palette("Rushmore1")
pal3 <- wes_palette("Darjeeling1")
pal4 <- wes_palette("FantasticFox1")
my_palette <- c(pal1,pal2,pal3,pal4)
############################################
# Part to plot cell length distribution
# including rays
############################################
cell_lengths <- read.csv("Data/cell_lengths_notConverge.csv")
files <- levels(as.factor(cell_lengths$Sample)) #create vector of file info
#create vector with species to bind to files 
species<- c("E. peritropoides","E. peritropoides","E. bracteata","E. lomelli","E. colligata","E. coalcomanensis",
            "E. coalcomanensis","E. calcarata","E. calcarata","E. finkii","E. calcarata","E. tithymaloides",
            "E. conzattii","E. cyri","E. peritropoides", "E. cymbifera","E. tehuacana","L-systemMesic",
            "L-systemXeric","E. diazlunana","E. diazlunana","E. diazlunana","E. lomelli",
            "E. cyri", "E. cymbifera", "E. tithymaloides","E. tithymaloides","E. personata",
            "E. personata","probL-system","probetaL-system","RayL-system") 
#create vector of growth form
habit <- c("mesic","mesic","xeric","xeric","mesic","mesic",
           "mesic","mesic","mesic","mesic","mesic","xeric",
           "mesic","xeric","mesic","xeric","xeric","L-systemCSM",
           "L-systemCSX","xeric","xeric","xeric","xeric",
           "xeric","xeric","xeric","xeric","xeric",
           "xeric","probL-system","probetaL-system","RayL-system")
species.id <- as.data.frame(cbind(files,species,habit))
#####
match.id <- match(cell_lengths$Sample,species.id$files)
cell_lengths$species <- as.character(species.id$species[match.id])
#
cell_lengths$habit <- as.character(species.id$habit[match.id]) #add growth form to cell lengths
agg <- aggregate(Number.of.cells ~ species, cell_lengths, function(x){
  qq <- quantile(x, probs = c(1, 3)/4)
  iqr <- diff(qq)
  lo <- qq[1] - 1.5*iqr
  hi <- qq[2] + 1.5*iqr
  c(Mean = mean(x), IQR = unname(iqr), lower = lo, high = hi)
}) 
#####
p  <- ggplot(cell_lengths, aes(Number.of.cells, colour=species, fill=species))
p  <- p + geom_density(alpha=0.2)
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
#Make two boxplots in horizontal and vertical axes
boxplot(cell_lengths$Number.of.cells~species, cell_lengths,horizontal=TRUE,
        col=new_pal)
boxplot(cell_lengths$Number.of.cells~cell_lengths$species,
        col=new_pal)
#Remove cell lengths minor to 3 cells
cell_lengths <- subset(cell_lengths, cell_lengths$Number.of.cells > 3)
cell_lengths$species <- with(cell_lengths, reorder(species,Number.of.cells,mean))
##
boxplot(cell_lengths$Number.of.cells~cell_lengths$habit,
        col=new_pal)
##
summary(cell_lengths$Number.of.cells)
cellfilemean<-aggregate(cell_lengths$Number.of.cells, list(cell_lengths$species), FUN=mean)
cellfilemean[with(cellfilemean, order(x)),]
numbercells<-aggregate(cell_lengths$Number.of.cells, list(cell_lengths$Sample), FUN=mean)
##
q  <- ggplot(cell_lengths, aes(Number.of.cells, colour=habit, fill=habit))
q  <- q + geom_density(alpha=0.2)
q
#Check type of length distribution with descdist
descdist(cell_lengths$Number.of.cells[cell_lengths$habit=="xeric"], discrete = FALSE)
descdist(cell_lengths$Number.of.cells[cell_lengths$habit=="mesic"], discrete = FALSE)
#Chek violin plot of cell lengths
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
#Plot cell frequencies ####### 
tipos_celulares <- read.csv("Data/word_counts_all/wordcounts_all.csv")
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

###########Graficar los indices de homogeneídad ##########
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
#####
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

####
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
