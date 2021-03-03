##Load libraries
library(dplyr)
library(ggplot2)

#Load file with coordinates
chaos_representation <- read.csv("Data/chaos_representation.csv", row.names = 1)

files <- levels(as.factor(chaos_representation$Sample))

pdf("Figures/chaos_representation.pdf")
for(i in files) {
  temp <- subset(chaos_representation, chaos_representation$Sample == i)
  plot(temp$X,temp$Y,main= i)
}
dev.off()
##First do it with subset
e.peritropoides.974 <- subset(chaos_representation, Sample == "974_edited_cells.txt")

e.perirtropoides.974 %>%
  group_by(count)

counter <- c()
y <- 0
x <- 1
for (i in 1:nrow(e.peritropoides.974)){
  y <- e.peritropoides.974$count[i]
  if (i == nrow(e.peritropoides.974))
  { 
    break
  } else if ((y + e.peritropoides.974$count[i+1]) > y)
  {
    counter <- c(counter, x)
  } else if((y + e.peritropoides.974$count[i+1]) == y)
  {  x <- x + 1
  }
}
tail(counter)
e.peritropoides.974$count[1] < e.peritropoides.974$count[1+1]

#Add species to the chaos data frame
species<- c("E. bracteata","E. lomelli","E. colligata","E. coalcomanensis",
            "E. calcarata","E. calcarata","E. finkii","E. conzattii","E. cyri",
            "E. peritropoides", "E. cymbifera","E. tehuacana","L-systemMesic",
            "L-systemXeric","E. diazlunana","E. diazlunana","E. diazlunana",
            "E. tithymaloides","E. tithymaloides","E. personata","E. personata",
            "L-system","RayL-system")
habit <- c("xeric","xeric","mesic","mesic",
           "mesic","mesic","mesic","mesic","xeric",
           "mesic","xeric","xeric","L-systemCSM","L-systemCSX","xeric",
           "xeric","xeric","xeric","xeric","xeric","xeric","L-system","RayL-system")

##Make data frame for files
species.id <- as.data.frame(cbind(files,species,habit))
#Create mathc between chaos and species
match.id <- match(chaos_representation$Sample,species.id$files)
#
chaos_representation$sp <- as.character(species.id$species[match.id])
#    
chaos_representation$habit <- as.character(species.id$habit[match.id])
#Plot with species
species.factor <- levels(as.factor(species))
pdf("Figures/chaos_representation_bysp.pdf")
for(i in species.factor) {
  temp <- subset(chaos_representation, chaos_representation$sp == i)
  plot(temp$X,temp$Y,main= i,col=as.factor(chaos_representation$Sample))
}
dev.off()

#Plot by habit
habit.factor <- levels(as.factor(habit))
pdf("Figures/chaos_representation_byhabit.pdf")
for(i in habit.factor) {
  temp <- subset(chaos_representation, chaos_representation$habit == i)
  plot(temp$X,temp$Y,main= i,col=as.factor(chaos_representation$Sample))
}
dev.off()
