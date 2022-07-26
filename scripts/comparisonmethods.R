###Script to compare both methods
comparison <- read.csv("Data/Images/P_tithymaloides/EPM6/ComparisonsDF.csv")

comparison$lineagelength <- nchar(comparison$CellNumber)
t.test(comparison$lineagelength ~ comparison$Method)
pdf("Figures/length_methods.pdf")
png("Figures/length_methods.png")
boxplot(comparison$lineagelength ~ comparison$Method, xlab = "Method",
        ylab="Lineage length")
dev.off()

live <- subset(comparison, comparison$Method == "Live")
image <- subset(comparison, comparison$Method == "Image")
tempmatcher <- match(image$Lineage,live$Lineage)
image$lineageLive <- live$CellNumber[tempmatcher]
image$lineageLiveLength <- live$lineagelength[tempmatcher]

image$lengthdiference <- abs(image$lineagelength - image$lineageLiveLength)
table(image$lengthdiference)
table((image$lengthdiference > 5))
library(dplyr)

write.csv(image$CellNumber,"Data/Images/P_tithymaloides/EPM6/comparison/image.txt",row.names = FALSE)
write.csv(image$lineageLive, "Data/Images/P_tithymaloides/EPM6/comparison/live.txt",row.names = FALSE)
