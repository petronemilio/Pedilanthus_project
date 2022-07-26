#######Loading packages required for ecologic analysis################
library(letsR)
library(dismo)
library(BIEN)
library(raster)
library(ggplot2)
library(plotbiomes)
library(kgc)
library(gridExtra)
library(dplyr)
library(tidyr)
library(factoextra)
library(sf)
library(maps)
library(RColorBrewer)
#######Extract ocurrences for all species of Carlquist database
ocurrences_dataframe <-data.frame()

species <- c("Euphorbia bracteata","Euphorbia colligata","Euphorbia diazlunana",
  "Euphorbia calcarata", "Euphorbial coalcomanensis", "Euphorbia conzattii",
  "Euphorbia lomelii", "Euphorbia tehuacana", "Euphorbia cyri","Euphorbia peritropoides",
  "Euphorbia personata", "Euphorbia tithymaloides","Euphorbia cymbifera","Euphorbia finkii")
for (i in species){
  temporal <- BIEN_occurrence_species(i, political.boundaries = T)
  ocurrences_dataframe <- rbind(ocurrences_dataframe, temporal)
}
####Associate climate variables
length(unique(ocurrences_dataframe$scrubbed_species_binomial))
bioclim.data <- getData(name = "worldclim",var = "bio",
                        res = 0.5, lon=-100,lat=20)
bioclim.data2 <- getData(name = "worldclim",var = "bio",
                     res = 0.5, lon=-90,lat=15)
class(bioclim.data)
#
r <- bioclim.data[[c(1,2,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19)]]
names(r) <- c("Temp","MeanDiurnalRange","Temp-Seasonality","MaxTempWarmestMonth",
              "MinTempColdestMonth","Temp-AnnualRange","MeanTempWettestQrt","MeanTempDryestQrt",
              "MeanTempWarmestQrt","MeanTempColdestQrt","Ann-Prec","PrecWetMonth",
              "PrecDriestMonth","PrecWetestQuarter","PrecDriestQuarter",
              "PrecWarmestQuarter","PrecColdestQuarter")
#######Not latitude and longitude. Import file from meta
##Use points only from not cultivated places
pedilanthus_notcultivated <- read.csv("meta/geo_samples.csv")
coords <- data.frame(x = as.numeric(pedilanthus_notcultivated$Long_dec), 
                     y =  as.numeric(pedilanthus_notcultivated$Lat_dec))
points <- SpatialPoints(coords, proj4string = r@crs)
values <- raster::extract(r,points)
##
pedilanthus_notcultivated <- cbind.data.frame(pedilanthus_notcultivated,values)
####trying to use the bioclim2 points
r <- bioclim.data2[[c(1,2,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19)]]
names(r) <- c("Temp","MeanDiurnalRange","Temp-Seasonality","MaxTempWarmestMonth",
              "MinTempColdestMonth","Temp-AnnualRange","MeanTempWettestQrt","MeanTempDryestQrt",
              "MeanTempWarmestQrt","MeanTempColdestQrt","Ann-Prec","PrecWetMonth",
              "PrecDriestMonth","PrecWetestQuarter","PrecDriestQuarter",
              "PrecWarmestQuarter","PrecColdestQuarter")
points <- SpatialPoints(coords, proj4string = r@crs)
values <- raster::extract(r,points)
values<-na.omit(values)
#add values to personata
pedilanthus_notcultivated[24,c(9:25)]<-values
########
#ocurrences_dataframe_mediandist_sp$Temp <- (ocurrences_dataframe_mediandist_sp$Temp)/10
pedilanthus_notcultivated <- pedilanthus_notcultivated %>%
  mutate_at(vars(Temp:PrecColdestQuarter), funs(./10))
write.csv(pedilanthus_notcultivated,"meta/PedilanthusWorldClim.csv")
#############################################
###### Plotting biomes Withaker biomes ######
paleta <- c("#d39e47","#b05ac8","#56b348","#d14184","#52a676","#d54340","#49afcf","#ca6d3c",
             "#6c79cc","#a2b147","#af69a6","#77732e","#e7879e","#a84e52")
wit <- whittaker_base_plot() + aes(alpha=0.5)
witreduced <- wit + theme(legend.position = "none")
all <- wit + geom_point(data = pedilanthus_notcultivated, 
                        aes(Temp, Ann.Prec),alpha=0.3)
IDentifier <- wit + geom_point(data = pedilanthus_notcultivated, 
              aes(Temp, Ann.Prec, colour = factor(Especie),shape=factor(Habit)),
              alpha=1)
IDentifier
ggsave("Figures/Whittaker.pdf",IDentifier, device = "pdf")
all
######
#plot whitaker
IDentifier <- wit + geom_point(data = pedilanthus_notcultivated, 
                               aes(Temp, Ann.Prec, colour = factor(Habit)))
IDentifier
##### Plotting maps ####
# show map with Latitude 200 as center
mex<-maps::map('world',xlim = c(-118,-85),ylim=c(14,33))
# add axes
mex
points(coords$x, coords$y)
####Reading shape files extracted and dowinloaded from conabio page
clim_mex <- read_sf("meta/clima1mgw/clima1mgw.shp")
temp_mex <- read_sf("meta/isotm1mgw/isotm1mgw.shp")
msnm_mex <- read_sf("meta/hipso4mgw/hipso4mgw.shp")
prec_mex <- read_sf("meta/preci4mgw/preci4mgw.shp")
#####Trying with ggplot 
nb.cols <- 12
mycolors <- colorRampPalette(brewer.pal(8, "Greens"))(nb.cols)
msnm_mex$RANGO <- factor(msnm_mex$RANGO, 
                levels = c("0 a 200","200 a 500","500 a 1000","1000 a 1500",
                          "1500 a 2000","2000 a 2500","2500 a 3000","3000 a 3500",
                          "3500 a 4000","4000 a 4500","4500 a 5000", "> 5000"))
altitud <- ggplot() + geom_sf(data = msnm_mex, aes(fill = RANGO),lwd=0.1) + theme_bw()+
  scale_fill_manual(values=mycolors) + theme(text = element_text(size=7)) +
  geom_point(data = pedilanthus_notcultivated, 
             aes(Long_dec, Lat_dec,colour=factor(Especie),
                 shape=factor(Habit)))
#
nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Reds"))(nb.cols)
temp_mex$TA_RANGO <- factor(temp_mex$TA_RANGO, levels = c("MENOR DE -2","DE -2 A 5", "DE 5 A 6", "DE 6 A 8",
           "DE 8 A 10","DE 10 A 12","DE 12 A 14","DE 14 A 16","DE 16 A 18","DE 18 A 20",
           "DE 20 A 22","DE 22 A 24","DE 24 A 26","DE 26 A 28","MAYOR DE 28"))
#
temperatura <- ggplot() + geom_sf(data = temp_mex, aes(fill = TA_RANGO),lwd=0)+ theme_bw()+
  scale_fill_manual(values=mycolors)+
  theme(text = element_text(size=7)) + geom_point(data = pedilanthus_notcultivated, 
             aes(Long_dec, Lat_dec,colour=factor(Especie),
                 shape=factor(Habit)))
#
nb.cols<- 10
mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(nb.cols)
prec_mex$RANGOS <- factor(prec_mex$RANGOS, levels = c("0 a 125 mm","125 a 400 mm","400 a 600 mm",
                          "600 a 800 mm","800 a 1200 mm", "mas de 4000 mm","1200 a 1500 mm",
                          "1500 a 2000 mm","2000 a 2500 mm","2500 a 4000 mm"))
precipitacion <- ggplot()+ geom_sf(data = prec_mex, aes(fill=RANGOS),lwd=0.01)+ 
  scale_fill_manual(values=mycolors) + theme(text = element_text(size=7))+
  geom_point(data = pedilanthus_notcultivated, 
             aes(Long_dec, Lat_dec,colour=factor(Especie),
                 shape=factor(Habit)))
####
ggsave("Figures/Pedilanthus_altitude.pdf",altitud, width = 20,height = 20,units = "cm")
ggsave("Figures/Pedilanthus_temperature.pdf",temperatura, width = 20,height = 20,units = "cm")
ggsave("Figures/Pedilanthus_precipitation.pdf",precipitacion, width = 20,height = 20,units = "cm")

##########PCA #########
complete.cases(pedilanthus_notcultivated)
boxplot(pedilanthus_notcultivated$Temp, 
        pedilanthus_notcultivated$MeanTempWarmestQrt,pedilanthus_notcultivated$MinTempColdestMonth,
        pedilanthus_notcultivated$Ann.Prec,pedilanthus_notcultivated$Temp.AnnualRange,
        pedilanthus_notcultivated$PrecWetMonth,pedilanthus_notcultivated$PrecDriestMonth)
rownames(pedilanthus_notcultivated)<- paste0(substr(pedilanthus_notcultivated$Especie, 11,14),pedilanthus_notcultivated$Colecta)
#Check variables useful for pca analysis
pairs(pedilanthus_notcultivated[,c(9,14,17,18,19,20,21)])
pairs(scale(pedilanthus_notcultivated[,c(9,14,17,18,19,20,21)]))
var(pedilanthus_notcultivated[,c(9,14,17,18,19,20,21)])
cor(pedilanthus_notcultivated[,c(9,14,17,18,19,20,21)])
cov(pedilanthus_notcultivated[,c(9,14,17,18,19,20,21)])
pca.pedilanthus <- prcomp(pedilanthus_notcultivated[,c(9,14,17,18,19,20,21)], scale=TRUE,center=TRUE,retx=TRUE)
summary(pca.pedilanthus)
pca.pedilanthus$rotation #eigenvectores
rownames(pca.pedilanthus$x)
pdf("Figures/PCA_climate.pdf")
ggbiplot(pca.pedilanthus,ellipse=TRUE,  labels=rownames(pca.pedilanthus), 
         groups=pedilanthus_notcultivated$Habit)
dev.off()
png("Figures/PCA_climate.png")
ggbiplot(pca.pedilanthus,ellipse=TRUE,  labels=rownames(pca.pedilanthus), 
         groups=pedilanthus_notcultivated$Habit)
dev.off()
####
traits_dataframe <- c()
for (i in species){
  temporal <- BIEN_trait_species(i, all.taxonomy=TRUE, political.boundaries = TRUE)
  traits_dataframe <- rbind(traits_dataframe, temporal)
}
####PithXylem
pedilanthus_samples <- read.csv("meta/samples.csv")
matcher <- match(pedilanthus_samples$Species, pedilanthus_notcultivated$Especie)
pedilanthus_samples$Habit<-pedilanthus_notcultivated$Habit[matcher]
#select onoe column
pedilanthus_samples$stemdiamcons <- with(pedilanthus_samples, 
                ifelse(is.na(stemdiametertotal),stemdiameter_reportedCacho,
                       stemdiametertotal))
stdiam_height.lm<- lm(log10(pedilanthus_samples$stemdiamcons)~
                        log10(pedilanthus_samples$Height.cm.))
summary(stdiam_height.lm)
plot(log10(pedilanthus_samples$stemdiamcons)~
       log10(pedilanthus_samples$Height.cm.))
points(log10(pedilanthus_samples$stemdiamcons[pedilanthus_samples$Habit=="Xeric"])~
         log10(pedilanthus_samples$Height.cm.[pedilanthus_samples$Habit=="Xeric"]),
       col="red")
points(log10(pedilanthus_samples$stemdiamcons[pedilanthus_samples$Habit=="Mesic"])~
         log10(pedilanthus_samples$Height.cm.[pedilanthus_samples$Habit=="Mesic"]),
       col="blue")
plot(log10(pedilanthus_samples$stemdiameter.xylem)~
       log10(pedilanthus_samples$Height.cm.))
pith.height<- lm(log10(pedilanthus_samples$pithdiameter)~ log10(pedilanthus_samples$Height.cm.))
summary(pith.height)
plot(log10(pedilanthus_samples$pithdiameter)~
       log10(pedilanthus_samples$Height.cm.))
points(log10(pedilanthus_samples$pithdiameter[pedilanthus_samples$Habit=="Xeric"])~
         log10(pedilanthus_samples$Height.cm.[pedilanthus_samples$Habit=="Xeric"]),
       col="red")
points(log10(pedilanthus_samples$pithdiameter[pedilanthus_samples$Habit=="Mesic"])~
         log10(pedilanthus_samples$Height.cm.[pedilanthus_samples$Habit=="Mesic"]),
       col="blue")

stdiamxylempith.lm<-lm(log10(pedilanthus_samples$stemdiameter.xylem)~
                         log10(pedilanthus_samples$Height.cm.))
summary(stdiamxylempith.lm)
plot(log10(pedilanthus_samples$stemdiameter.xylem)~
       log10(pedilanthus_samples$Height.cm.))
points(log10(pedilanthus_samples$stemdiameter.xylem[pedilanthus_samples$Habit=="Xeric"])~
         log10(pedilanthus_samples$Height.cm.[pedilanthus_samples$Habit=="Xeric"]),
       col="red")
points(log10(pedilanthus_samples$stemdiameter.xylem[pedilanthus_samples$Habit=="Mesic"])~
         log10(pedilanthus_samples$Height.cm.[pedilanthus_samples$Habit=="Mesic"]),
       col="blue")
abline(stdiam_height.lm)
#
pedilanthus_samples$xylempith.ratio <- pedilanthus_samples$stemdiameter.xylem/pedilanthus_samples$pithdiameter
xpratio.lm<-lm(log10(pedilanthus_samples$xylempith.ratio)~log10(pedilanthus_samples$Height.cm.))
summary(xpratio.lm)
plot(log10(pedilanthus_samples$xylempith.ratio)~
       log10(pedilanthus_samples$Height.cm.))
points(log10(pedilanthus_samples$xylempith.ratio[pedilanthus_samples$Habit=="Xeric"])~
         log10(pedilanthus_samples$Height.cm.[pedilanthus_samples$Habit=="Xeric"]),
       col="red")
points(log10(pedilanthus_samples$xylempith.ratio[pedilanthus_samples$Habit=="Mesic"])~
         log10(pedilanthus_samples$Height.cm.[pedilanthus_samples$Habit=="Mesic"]),
       col="blue")
boxplot(pedilanthus_samples$xylempith.ratio ~ pedilanthus_samples$Habit)
pedilanthus_samples$stemdiametertotal-pedilanthus_samples$pithdiameter
####
plot(log10(pedilanthus_samples$stemdiameter.xylem)~ log10(pedilanthus_samples$Height.cm.))
plot(log10(pedilanthus_samples$pithdiameter)~ log10(pedilanthus_samples$Height.cm.))
pedilanthus_samples$stemdiameter.xylem - pedilanthus_samples$pithdiameter
