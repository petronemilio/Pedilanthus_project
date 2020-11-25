#Para poner tablas en formato wide agregar bibliotecas 
library(tidyr)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(ggstance)
library(vioplot)
##sET WORKING DIRECTORY
setwd("Doctorado/scripts/")
###
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
           stat="identity") 
  geom_bar(position = "fill") + scale_y_continuous(labels = count)

#ggplot(data=dat1, aes(x=time, y=total_bill, fill=sex)) +
#  geom_bar(stat="identity")

##########Graficar los indices de homogeneídad
##########
p_bracteatus<-read.csv("../Data/Pedilanthus/P_bracteatus/cell_homogeneity.txt", sep = " ")
p_bracteatus<-p_bracteatus[complete.cases(p_bracteatus), ]
#
p_calcarata_896<-read.csv("../Data/Pedilanthus/P_calcaratus/896_cell_homogeneity.txt",sep=" ")
p_calcarata_892<-read.csv("../Data/Pedilanthus/P_calcaratus/892_cell_homogeneity.txt",sep=" ")
p_calcarata_896<-p_calcarata_896[complete.cases(p_calcarata_896), ]
p_calcarata_892<-p_calcarata_892[complete.cases(p_calcarata_892), ]
#
p_colligata_867<- read.csv("../Data/Pedilanthus/P_colligata/867_cell_homogeneity.txt",sep=" ")
p_colligata_867<-p_colligata_867[complete.cases(p_colligata_867), ]
#
p_diazluna_epm10<-read.csv("../Data/Pedilanthus/P_diazluna/EPM10_cell_homogeneity.txt", sep=" ")
p_diazluna_epm11<-read.csv("../Data/Pedilanthus/P_diazluna/EPM11_cell_homogeneity.txt", sep=" ")
p_diazluna_epm12<-read.csv("../Data/Pedilanthus/P_diazluna/EPM12_cell_homogeneity.txt", sep=" ")
p_diazluna_epm10<-p_diazluna_epm10[complete.cases(p_diazluna_epm10), ]
p_diazluna_epm11<-p_diazluna_epm11[complete.cases(p_diazluna_epm11), ]
p_diazluna_epm12<-p_diazluna_epm12[complete.cases(p_diazluna_epm12), ]
#
p_finkii917<-read.csv("../Data/Pedilanthus/P_finkii/917_cell_homogeneity.txt", sep=" ")
p_finkii917<-p_finkii917[complete.cases(p_finkii917), ]
#
p_tithy_5<-read.csv("../Data/Pedilanthus/P_tithymaloides/EPM5_cell_homogeneity.txt", sep=" ")
p_tithy_6<-read.csv("../Data/Pedilanthus/P_tithymaloides/EPM6_cell_homogeneity.txt", sep=" ")
p_tithy_5<-p_tithy_5[complete.cases(p_tithy_5), ]
p_tithy_6<-p_tithy_6[complete.cases(p_tithy_6), ]
#
p_tomentellus<-read.csv("../Data/Pedilanthus/P_tomentellus/973_cell_homogeneity.txt",sep=" ")
p_tomentellus<-p_tomentellus[complete.cases(p_tomentellus),]

########
p_bracteatus<-cbind(rep("E. bracteata",nrow(p_bracteatus)), rep("845", nrow(p_bracteatus)),p_bracteatus)
#
p_calcarata_892<-cbind(rep("E. calcarata",nrow(p_calcarata_892)), rep("892", nrow(p_calcarata_892)),p_calcarata_892)
p_calcarata_896<-cbind(rep("E. calcarata",nrow(p_calcarata_896)), rep("896", nrow(p_calcarata_896)),p_calcarata_896)
#
p_colligata_867<-cbind(rep("E. colligata",nrow(p_colligata_867)), rep("867", nrow(p_colligata_867)),p_colligata_867)
#
p_diazluna_epm10<-cbind(rep("E. diazlunana",nrow(p_diazluna_epm10)),rep("10",nrow(p_diazluna_epm10)),p_diazluna_epm10)
p_diazluna_epm11<-cbind(rep("E. diazlunana",nrow(p_diazluna_epm11)),rep("11",nrow(p_diazluna_epm11)),p_diazluna_epm11)
p_diazluna_epm12<-cbind(rep("E. diazlunana",nrow(p_diazluna_epm12)),rep("12",nrow(p_diazluna_epm12)),p_diazluna_epm12)
#
p_finkii917<-cbind(rep("E. finkii",nrow(p_finkii917)),rep("917",nrow(p_finkii917)),p_finkii917)
#
p_tithy_5<-cbind(rep("E. tithymaloides",nrow(p_tithy_5)),rep("5",nrow(p_tithy_5)),p_tithy_5)
p_tithy_6<-cbind(rep("E. tithymaloides",nrow(p_tithy_6)),rep("6",nrow(p_tithy_6)),p_tithy_6)
#
p_tomentellus<-cbind(rep("E. cyri",nrow(p_tomentellus)),rep("5",nrow(p_tomentellus)),p_tomentellus)
#
###############
names(p_bracteatus) <- c("Especie","Individuo","Living/Death", "Conductive/Non-conductive")
names(p_calcarata_896) <- c("Especie","Individuo","Living/Death", "Conductive/Non-conductive")
names(p_calcarata_892) <- c("Especie","Individuo","Living/Death", "Conductive/Non-conductive")
names(p_colligata_867) <- c("Especie","Individuo","Living/Death", "Conductive/Non-conductive")
names(p_diazluna_epm10) <- c("Especie","Individuo","Living/Death", "Conductive/Non-conductive")
names(p_diazluna_epm11) <- c("Especie","Individuo","Living/Death", "Conductive/Non-conductive")
names(p_diazluna_epm12) <- c("Especie","Individuo","Living/Death", "Conductive/Non-conductive")
names(p_finkii917) <- c("Especie","Individuo","Living/Death", "Conductive/Non-conductive")
names(p_tithy_5) <- c("Especie","Individuo","Living/Death", "Conductive/Non-conductive")
names(p_tithy_6) <- c("Especie","Individuo","Living/Death", "Conductive/Non-conductive")
names(p_tomentellus) <- c("Especie","Individuo","Living/Death", "Conductive/Non-conductive")


pedilanthus<-rbind(p_calcarata_892,p_calcarata_896,p_bracteatus,
                   p_colligata_867,p_diazluna_epm10,p_diazluna_epm11,
                   p_diazluna_epm12, p_finkii917, p_tithy_5,p_tithy_6,
                   p_tomentellus)

ggplot(pedilanthus, aes(x=`Living/Death`, y=`Conductive/Non-conductive`, colour=`Especie`))+
  geom_point()+ 
  scale_color_manual(values=c("#c80966","#84de66","#8369e1","#d89816",
                              "#0152a1","#005e18","#abb2ff"))+
  theme_bw()
scale_

