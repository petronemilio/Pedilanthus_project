#Para poner tablas en formato wide agregar bibliotecas 
library(tidyr)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(ggstance)
library(vioplot)
library(plot3D)
library(scatterplot3d)
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
p_brac_long<-read.csv("../Data/Pedilanthus/P_bracteatus/bracteatus_length_filecells.csv",
                        header = TRUE)
#
p_calc_long<-read.csv("../Data/Pedilanthus/P_calcaratus/calcaratus_length_filecells.csv",
                               header = TRUE)
#E. colli
p_col_long<-read.csv("../Data/Pedilanthus/P_colligata/colligata_length_filecells.csv",
                         header = TRUE)
#E. coalcomanensis
p_coalco_long<-read.csv("../Data/Pedilanthus/P_coalcomanensis/coalcomanensis_length_filecells.csv",
                        header = TRUE)
#E. cymbifera
p_cymbi_long<-read.csv("../Data/Pedilanthus/P_cymbiferus/cymbiferus_length_filecells.csv",
                       header = TRUE)

######
p_diazluna_long<-read.csv("../Data/Pedilanthus/P_diazluna/diazluna_length_filecells.csv",
                             header = TRUE)

#
p_finkii_long<-read.csv("../Data/Pedilanthus/P_finkii/finkii_length_filecells.csv",
                        header = TRUE)
#
p_macrocarpus_long<-read.csv("../Data/Pedilanthus/P_macrocarpus/macrocarpus_length_filecells.csv",
                             header = TRUE)
#
p_peritropoides_long <-read.csv("../Data/Pedilanthus/P_peritropoides/peritropoides_length_filecells.csv",
                               header = TRUE)
#
p_personata_long <- read.csv("../Data/Pedilanthus/P_personata/personata_length_filecells.csv",
                             header = TRUE)

p_pulchellus_long<-read.csv("../Data/Pedilanthus/P_pulchellus/pulchellus_length_filecells.csv",
                            header = TRUE)

#####
p_tehuacanus_long<-read.csv("../Data/Pedilanthus/P_tehuacanus/tehuacanus_length_filecells.csv",
                            header = TRUE)
#####
p_tithy_long<-read.csv("../Data/Pedilanthus/P_tithymaloides/tithymaloides_length_filecells.csv",
                         header = TRUE)

######
p_tomentellus_long<-read.csv("../Data/Pedilanthus/P_tomentellus/tomentellus_length_filecells.csv",
                             header = TRUE)


longitud_series<-rbind(p_brac_long,p_brac_long, p_calc_long, p_coalco_long,
                       p_col_long, p_cymbi_long, p_diazluna_long,p_finkii_long,
                       p_macrocarpus_long, p_peritropoides_long, p_personata_long,
                       p_pulchellus_long,p_tehuacanus_long, p_tithy_long,
                       p_tomentellus_long)

p  <- ggplot(longitud_series, aes(Number.of.cells, colour=Species, fill=Species))
p  <- p + geom_density(alpha=0.2)
p

pal<-c(wes_palette("Cavalcanti1"),wes_palette("GrandBudapest1"),wes_palette("GrandBudapest2"),
       wes_palette("FantasticFox1"))
par(mar=c(7,5,1,1))
new_pal<-c("#cb57aa","#6db744","#8760cf","#c2af45","#7d7fc5","#dd8d4a","#45b0cf","#cf483c",
           "#5cc08c","#c56179","#3b824e","#976530","#76853a")

boxplot(Number.of.cells~Species, longitud_series,horizontal=TRUE,
        col=new_pal)

longitud_series$Species<-with(longitud_series, reorder(Species, Number.of.cells, mean))

boxplot(Number.of.cells~Species, longitud_series,horizontal = TRUE, las=1.5,
        col=c(pal), xlab="Longitud de las filas",ylab = NULL, cex=0.4,
        names=c(expression(italic("E. lomelli")),expression(italic("E. cymbifera")),
                expression(italic("E. tithymaloides")), expression(italic("E. diazlunana")),
                expression(italic("E. personata")), expression(italic("E. finkii")),
                expression(italic("E. bracteata")),expression(italic("E. coalcomanensis")),
                expression(italic("E. peritropoides")),expression(italic("E. colligata")),
                expression(italic("E. cyri")),expression(italic("E. tehuacana")),
                expression(italic("E. calcarata")),expression(italic("E. conzattii"))),
        cex.names=1.1, cex.axis=0.8, cex.lab=1.2)

boxplot(Number.of.cells~Species, longitud_series)

# Draw the boxplot using this new order
boxplot(data$note ~ new_order , ylab="sickness" , col="#69b3a2", boxwex=0.4 , main="")


ggplot(longitud_series, aes(x =Species,y=Number.of.cells, fill = Species)) +
  geom_violin()+ coord_flip()+
  scale_y_continuous(breaks = c(0,50,100,150,200,250,300,350,400,500,600))+
  scale_fill_manual(values=c(pal))+
  theme_classic()+  theme(legend.position = "none") 
 

# Draw the plot
calca <- subset(longitud_series, Especie=="E. calcarata", select=c(Especie, Longitud))
tithy <- subset(longitud_series, Especie=="E. tithymaloides", select=c(Especie, Longitud))
cyri <- subset(longitud_series, Especie=="E. cyri", select=c(Especie, Longitud))
brac <- subset(longitud_series, Especie=="E. bracteata", select=c(Especie, Longitud))
finkii<- subset(longitud_series,Especie=="E. finkii", select=c(Especie, Longitud))
colli<- subset(longitud_series,Especie=="E. colligata", select=c(Especie, Longitud))
diazlunana<- subset(longitud_series,Especie=="E. diazlunana", select=c(Especie, Longitud))
#
vioplot(brac$Longitud,calca$Longitud, colli$Longitud,cyri$Longitud,
        finkii$Longitud,diazlunana$Longitud,tithy$Longitud, 
        names=c("a","b","c","d","e","f","g"),
        col=c("#899DA4","#C93312","#FAEFD1","#DC863B", "#F1BB7B", "#FD6467","#5B1A18"))

############################################################
#Graficar los tipos celulares
tipos_celulares <- read.csv("../Data/Pedilanthus/cell_counts.csv")

especies <- levels(longitud_series$Species)
especies
especies_por_muestra <- c(especies[1], rep(especies[2],2), especies[3:5], rep(especies[6],3), 
  especies[7:8],rep(especies[10],2),especies[9],especies[11:12], rep(especies[13],2), especies[14])

tipos_celulares <- cbind(especies_por_muestra,tipos_celulares)

tipos_celulares$totalcells <-rowSums(tipos_celulares[,c(3:6)])

tipos_celulares_freq <- tipos_celulares[,c(3:6)] / tipos_celulares[,7] 
tipos_celulares_freq$especie <- tipos_celulares$especies_por_muestra

tipos_celulares <- reshape(data=tipos_celulares, idvar="x", varying = c("Vasos","Fibras","Parénquima","Radios"),
        times=c("Vasos","Fibras","Parénquima","Radios"),v.name=c("Count"),direction="long")


ggplot(tipos_celulares) + 
  geom_bar(mapping=aes(x = especies_por_muestra, y = Count, fill= time), stat = "identity",
           position = "fill") +coord_flip()

  ggplot(tipos_celulares) + 
  geom_bar(mapping=aes(x = especies_por_muestra, y = Count, fill= time), stat = "identity")

  
######################
library(data.table)
#setDT converts to a data.table and then you calculate the fraction of each expr
#grouping by the transcript_id
#setDT(tipos_celulares)[, frac := Count / sum(Count), by=X]

#("gridExtra")

tipos_fusiformes <- subset(tipos_celulares, time != "Radios")
ggplot(tipos_fusiformes) + 
  geom_bar(mapping=aes(x = especies_por_muestra, y = Count, fill= time), stat = "identity",
           position = "fill")

ggplot(tipos_fusiformes) + 
  geom_bar(mapping=aes(x = especies_por_muestra, y = Count, fill= time), stat = "identity")

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
##Make data frame for files
species.id <- as.data.frame(cbind(files,species,habit))
match.id <- match(homogenity_index$X0,species.id$files)
#
homogenity_index$sp <- as.character(species.id$species[match.id])
#    
homogenity_index$habit <- as.character(species.id$habit[match.id])
  
ggplot(homogenity_index, aes(x=storage, y=conductivity, colour=sp))+
  geom_point() #+ scale_color_manual(values=my_palette)

pdf("Figures/storage_conductivity.pdf") # Para guardar en PDF
ggplot(homogenity_index, aes(x=storage, y=conductivity, colour=habit))+
  geom_point() + scale_color_manual(values=pal3)
dev.off()
#
pdf("Figures/storage_support.pdf") # Para guardar en PDF
ggplot(homogenity_index, aes(x=storage, y=support, colour=habit))+
  geom_point() + scale_color_manual(values=pal3)
dev.off()
#
pdf("Figures/support_conductivity.pdf") # Para guardar en PDF
ggplot(homogenity_index, aes(x=conductivity, y=support, colour=habit))+
  geom_point() +  scale_color_manual(values=pal3)
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

