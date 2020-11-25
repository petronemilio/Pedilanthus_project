#Load useful libraries
library(tidyr)
library(ggplot2)
#Set working directory 
setwd("/home/emiliopetrone/Doctorado/scripts")

#Comparar la segmentación manual contra los diferentes filtros en 10X
manualcounts15 <- read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/EPM6_015_manual_counts.csv",header = TRUE)

manualcounts15 <- manualcounts15[,c(2:5)]
manualcountstotal15 <- rowSums(manualcounts15)

#Cargar la base de datoscon centroides definidos en los diferentes métodos de segmentación.
segmentation10x <- read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/prueba/parameters_segmentation.csv",header = TRUE)

segmentation10x <- separate(segmentation10x, col = Segmentation, into = c("filter","image"), sep = "_")

segmentation15 <- subset(segmentation10x, image == 15)

table(segmentation15$filter)

#Graficar perímetros, areas y demás de las células
centroids<-read.csv("../Data/Pedilanthus/P_tithymaloides/EPM6_S2_centroids.csv", header = FALSE)
area<-read.csv("../Data/Pedilanthus/P_tithymaloides/EPM6_S2_areas.csv", header = FALSE)
perimeter<-read.csv("../Data/Pedilanthus/P_tithymaloides/EPM6_S2_perimeters.csv", header = FALSE)

centroids2<-read.csv("../Data/Pedilanthus/P_tithymaloides/EPM6_S2_smooth_centroids.csv", header = FALSE)
area2<-read.csv("../Data/Pedilanthus/P_tithymaloides/EPM6_S2_smooth_areas.csv", header = FALSE)
perimeter2<-read.csv("../Data/Pedilanthus/P_tithymaloides/EPM6_S2_smooth_perimeters.csv", header = FALSE)


boxplot(area)
boxplot(area2)
boxplot(perimeter)
boxplot(perimeter2)
hist(area$V1)
hist(area2$V1)
hist(perimeter$V1)
hist(perimeter2$V1)

#Graficar perímetros, areas y demás de las células_othercell
centroids15<-read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/EPM6_015_centroids.csv", header = FALSE)
area15<-read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/EPM6_015_areas.csv", header = FALSE)
perimeter15<-read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/EPM6_015_perimeters.csv", header = FALSE)

centroids15_2<-read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/EPM6_015_skel_centroids.csv", header = FALSE)
area15_2<-read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/EPM6_015_skel_areas.csv", header = FALSE)
perimeter15_2<-read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/EPM6_015_skel_perimeters.csv", header = FALSE)

boxplot(area)
boxplot(area2)
boxplot(perimeter)
boxplot(perimeter2)
hist(area$V1)
hist(area2$V1)
hist(perimeter$V1)
hist(perimeter2$V1)


#Load the csv file with parameters for different segm. methods.

segmentation <-read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/40X/prueba/segmentation_DF_40x001.csv")

manualcounts <- read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/40X/prueba/manual_count40x_001.csv")

manualcounts <- manualcounts[,c(2:6)]
rowSums(manualcounts)

str(segmentation)
colSums(segmentation$Segmentation == "distancefilt_binarythreshold")
table(segmentation$Segmentation)

segmentation_without0<- subset(segmentation, segmentation$areas > 50)
table(segmentation_without0$Segmentation)

boxplot(segmentation_without0$areas~segmentation_without0$Segmentation)

###############################################################################
###############################################################################
####  Cargar resultados sobre first series. El análisis de partículas de imageJ, 
####  de vecindad y de las series codificadas.

morfo11 <- read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/Results_fiji_labcolor_epm6_11.csv")

plot(morfo11$X, morfo11$Y)
plot(morfo11$Area, morfo11$Perim)
hist(morfo11$Round)
#Cargar las filas contadas 
cell_files11<-read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_011_cell_column.csv", 
         sep= "," ,header=TRUE)

matchid11 <- match(cell_files11$Cell_id, morfo11$X.1)

cell_files11 <- cbind(cell_files11, morfo11[matchid11, c(2,6:17)])
#
plot(cell_files11$ID_cell, cell_files11$Area)
plot(cell_files11$ID_cell, cell_files11$Circ.)
#
plot(cell_files11$X, cell_files11$Y, col= cell_files11$File)
###
cell_file1 <- subset(cell_files11, cell_files11$File == 1)

cell_file2 <- subset(cell_files11, cell_files11$File == 2)
hist(cell_file2$Circ.)
##
##########################Delauny triangulation
#morfo11 <- read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/Results_imagej_labcolor_epm6_11.csv")
plot(morfo11$X, morfo11$Y)

#Cargar datos delauny.
delauny11 <-read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_011_labcolor_editedimagej_delauny.csv",
                    header = TRUE)

delauny11 <- delauny11[,c(2:5)]
delauny11 <- cbind(delauny11,paste(delauny11$x1,delauny11$y1, sep = ","),paste(delauny11$x2,delauny11$y2, sep = ","))

colnames(delauny11)[5] <- "x1y1"
colnames(delauny11)[6] <- "x2y2"
##
roundx<-round(morfo11$X)
roundy<-round(morfo11$Y)

morfo11 <- cbind(morfo11, paste(roundx, roundy, sep= ","))

#rbind(delauny11$x1y1, delauny11$x2y2)
vecinitycells11<-as.data.frame(table(c(as.vector(delauny11$x1y1), as.vector(delauny11$x2y2))))
hist(vecinitycells11$Freq)
qplot(vecinitycells11$Freq, geom = "histogram", binwidth=1)
# Compute a histogram of `chol$AGE`
ggplot(data=vecinitycells11, aes(Freq)) + 
  geom_histogram(bins = 9, fill="#02401B", col="black")
# col="#972D15"

###Cargar 
morfo8 <- read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/Results_imagej_labcolor_epm6_08.csv")
plot(morfo8$X, morfo8$Y)
#Cargar las filas contadas 
cell_files8<-read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_008_cell_column.csv", 
                       sep= "," ,header=TRUE)
#
matchid8 <- match(cell_files8$Cell_id, morfo8$X.1)

cell_files8 <- cbind(cell_files8, morfo8[matchid8, c(2,6:17)])
#
plot(cell_files8$ID_cell, cell_files8$Area)
plot(cell_files8$ID_cell, cell_files8$Circ.)
#
plot(cell_files8$X, cell_files8$Y, col= cell_files8$File)

########
#Cargar datos delauny.
delauny8 <-read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_008_labcolor_editedimagej_delauny.csv",
                   header = TRUE)

delauny8 <- delauny8[,c(2:5)]
delauny8 <- cbind(delauny8,paste(delauny8$x1,delauny8$y1, sep = ","),paste(delauny8$x2,delauny8$y2, sep = ","))

colnames(delauny8)[5] <- "x1y1"
colnames(delauny8)[6] <- "x2y2"
##
roundx<-round(morfo8$X)
roundy<-round(morfo8$Y)

morfo8 <- cbind(morfo8, paste(roundx, roundy, sep= ","))
##
vecinitycells8<-as.data.frame(table(c(as.vector(delauny8$x1y1), as.vector(delauny8$x2y2))))

hist(vecinitycells8$Freq)
#outer(morfo8$roundx, delauny8$x1, "==") & outer(morfo8$roundy, delauny8$y1, "==")
ggplot(data=vecinitycells8, aes(Freq)) + 
  geom_histogram(bins = 9, fill="#976530", col="black")
ggplot(data=vecinitycells11, aes(Freq)) + 
  geom_histogram(bins = 9, fill="#02401B", col="black")

#Agregar identificador de la zona de la que proviene la foto
vecinity11 <- cbind(rep("Madera interna",nrow(vecinitycells11)), vecinitycells11)
vecinity8 <- cbind(rep("Madera externa",nrow(vecinitycells8)), vecinitycells8)

colnames(vecinity11)[1] <- "Región"
colnames(vecinity8)[1] <- "Región"

vecinity_all <- rbind(vecinity11,vecinity8)

ggplot(data=vecinity_all, aes(Freq, by=Región)) + 
  geom_histogram(bins = 9, col="black")

ggplot(vecinity_all, aes(x=Freq)) +
  geom_histogram(fill="#02401B", col="black", breaks=(seq(2,11, by=1)))

to_frequency <- as.data.frame(table(vecinity_all$Freq))
#table(vecinity_all$Freq, by(vecinity_all$Región))
to_frequency <- cbind(to_frequency, to_frequency$Freq/sum(to_frequency$Freq))
colnames(to_frequency)[3] <- "Relative Frequency"

plot(to_frequency$Var1, to_frequency$`Relative Frequency`,type="l", col="green")


ggplot(to_frequency, aes(x=to_frequency$Var1,y=to_frequency$`Relative Frequency`, group=1))+
       geom_point(size=1.5)+geom_line(color="#02401B", size=1.5)+
       xlab("Número de vecinas") + ylab("Frecuencia Relativa") 
      
######################################
###Cargar 
morfo7 <- read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/Results_fiji_labcolor_epm6_7.csv")
plot(morfo7$X, morfo7$Y)

########
#Cargar datos delauny.
delauny7 <-read.csv("../Data/Pedilanthus/P_tithymaloides/E_tithymaloides_fotos/first_series/EPM6_007_labcolor_editedimagej_delauny.csv",
                    header = TRUE)

delauny7 <- delauny7[,c(2:5)]
delauny7 <- cbind(delauny7,paste(delauny7$x1,delauny7$y1, sep = ","),paste(delauny7$x2,delauny7$y2, sep = ","))

colnames(delauny7)[5] <- "x1y1"
colnames(delauny7)[6] <- "x2y2"
##
roundx<-round(morfo7$X)
roundy<-round(morfo7$Y)

morfo7 <- cbind(morfo7, paste(roundx, roundy, sep= ","))
##
vecinitycells7 <-as.data.frame(table(c(as.vector(delauny7$x1y1), as.vector(delauny7$x2y2))))

ggplot(data=vecinitycells7, aes(Freq)) + 
  geom_histogram(bins = 9, fill="#976530", col="black")

####Frequency of the two zones
to_frequency7 <- as.data.frame(table(vecinitycells7$Freq))
#table(vecinity_all$Freq, by(vecinity_all$Región))
to_frequency7 <- cbind(to_frequency7, to_frequency7$Freq/sum(to_frequency7$Freq))
colnames(to_frequency7)[3] <- "Relative Frequency"
###
to_frequency11 <- as.data.frame(table(vecinitycells11$Freq))
to_frequency11 <- cbind(to_frequency11, to_frequency11$Freq/sum(to_frequency11$Freq))
#Agregar identificador de la zona de la que proviene la foto
to_frequency7 <- cbind(rep("Madera externa",nrow(to_frequency7)), to_frequency7)
to_frequency11 <-  cbind(rep("Madera interna",nrow(to_frequency11)), to_frequency11)
#Change col names
colnames(to_frequency11)[4] <- "Relative Frequency"
colnames(to_frequency11)[3] <- "Freq"

colnames(to_frequency7)[1] <- "Zone"
colnames(to_frequency11)[1] <- "Zone"

to_frequency <- rbind(to_frequency11,to_frequency7) 


plot(to_frequency$Var1, to_frequency$`Relative Frequency`,type="l", col="green")

ggplot(to_frequency, aes(x=to_frequency$Var1,y=to_frequency$`Relative Frequency`, color=Zone, group=Zone))+
  geom_line(size=1.5)+
  scale_linetype_manual(values=c("#02401B","#C93312"))+
  xlab("Número de vecinas") + ylab("Frecuencia Relativa") 

##################
#Ver como graficar las líneas codificadas.

par(mfrow=c(1,1))
plot(cell_files8$X, cell_files8$Y, col= cell_files8$File)
plot(cell_files11$X, cell_files11$Y, col= cell_files11$File)

plot(cell_files8$Area)
plot(cell_files8$Circ., type="l")
plot(cell_files11$Circ., type = "l")

#
?wes_palette

palette(c("#00287b","#82b730","#6a3fb1","#5bd167",
          "#9e0683","#00ae64","#e2398b","#fdce56",
          "#0159a4","#e1672b","#8f77bc","#b95100",
          "#9e0060","#e4d580","#dc2a64","#5f5600",
          "#c3103d","#ff907a","#9a1b00","#ff7b64")) 

plot(cell_files11$Width, cell_files11$Height, col = cell_files11$File, pch=19)
plot(cell_files8$Width, cell_files8$Height, col = cell_files8$File, pch= 19 )

plot(cell_files11$X, cell_files11$Y, col = cell_files11$File, pch=19)
plot(cell_files8$X, cell_files8$Y, col = cell_files8$File)

plot(cell_files8$Y, cell_files8$Area, col = cell_files8$File)
plot(cell_files11$Y, cell_files11$Area, col= cell_files11$File)
plot(cell_files8$Area, cell_files8$Perim., col=cell_files8$File)


cell_files8_sub <- subset(cell_files8, File > 20 & File < 40)

ggplot(cell_files8_sub, aes(x=ID_cell, y=Area, colour=as.factor(File))) + geom_point()+
  scale_color_manual(values=palette)



ggplot(cell_files8, aes(x=ID_cell, y=Area, colour=as.factor(File))) + geom_line()
ggplot(cell_files11, aes(x=ID_cell, y=Area, colour=as.factor(File))) + geom_line()

