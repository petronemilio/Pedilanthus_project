#Load file with coordinates
chaos_representation <- read.csv("Data/chaos_representation.csv", row.names = 1)


e.tithy.947 <- subset(chaos_representation, chaos_representation$Sample == '974_edited_cells.txt')
plot(e.tithy.947$X,e.tithy.947$Y)

files <- levels(as.factor(chaos_representation$Sample))

pdf("Figures/chaos_representation.pdf")
for(i in files) {
  temp <- subset(chaos_representation, chaos_representation$Sample == i)
  plot(temp$X,temp$Y,main= i)
}
dev.off()
