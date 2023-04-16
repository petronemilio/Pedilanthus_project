##########
#this script generates the plots and models of the 
# maximum number of words and the similitude based on word frequencies
#plots and tables are made
# to analize different measures of word diversity
# and complexity
#Load color blind pallete generated from https://medialab.github.io/iwanthue/
##### Load useful libraries #####
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(corrplot)
library(gplots)
library(wesanderson)
library(pheatmap)
library(vegan)
library(scales)
##### Load files
## Loading species id file
species.id <- read.csv("meta/speciesID.csv")
insert <- "_NotConverge"
species.id$filesconverge <- sub("(?<=\\w)(?=\\.)", insert, species.id$files, perl=TRUE)
## => A1.A1c<=7 A1.A1c<=8
#
#####Load word counts file to plot as a function of word length #####
wc <- read.csv("Data/word_counts_all/wordcounts_all.csv",stringsAsFactors = FALSE)
class(wc)
wc$X <-rep(c(1:34),each=32)
plot(wc$NumberOfWords~ wc$X)
###
nwords <- subset(wc, wc$file != "contextmesicLsystem_edited_cells_NotConverge.txt" &
                   wc$file != "contextxericLsystem_edited_cells_NotConverge.txt" &
                   wc$file != "Ray_Lsystem_edited_cells_NotConverge.txt")
matcher <- match(nwords$file, species.id$files)
nwords$species <- species.id$species[matcher]
nwords$habit <- species.id$habit[matcher]
#
write.csv(nwords,"Data/wc_resumen.csv")#define the correct folder!Maybe move it to meta!!!
plot(nwords$X,nwords$NumberOfWords)
nwordslog <- ggplot(data= nwords)+
             aes(x= X, y= log10(NumberOfWords),group = file,  color= habit) + 
  geom_line(aes(linetype=habit),size = 1)+
  scale_color_manual(values = c("#CB2314","#FAD510","#273046"))
#
nwordsnormal <- ggplot(data= nwords)+
  aes(x= X, y= NumberOfWords,group = file,  color= habit) + 
  geom_line(aes(linetype=habit),size = 1)+
  scale_color_manual(values = c("#CB2314","#FAD510","#273046"))
nwordsnormal
nwordslog
ggsave("Figures/nwords_log.pdf",nwordslog, device = "pdf")
ggsave("Figures/nwords_notlog.pdf",nwordsnormal, device = "pdf")
ggsave("Figures/nwords_log.png",nwordslog, device = "png")
ggsave("Figures/nwords_notlog.png",nwordsnormal, device = "png")
#
#remove specific factors from plot
nwords.temp <- subset(nwords, nwords$habit != "L-system")
#
nwords_sp <- ggplot(data= nwords.temp) + 
  aes(x= X, y= NumberOfWords, group = file, color= species) + 
  geom_line(aes(linetype=species),size = 1)

nwords_sphabit <- ggplot(data= nwords.temp)+
   aes(x= X, y= NumberOfWords, group = file, color= habit) + 
  geom_line(aes(linetype=habit),size = 1) +
   scale_color_manual(values = c("#FAD510","#273046"))

ggsave("Figures/nwords_sp.pdf",nwords_sp, device = "pdf")
ggsave("Figures/nwords_sphabit.pdf",nwords_sphabit, device = "pdf")
ggsave("Figures/nwords_sp.png",nwords_sp, device = "png")
ggsave("Figures/nwords_sphabit.png",nwords_sphabit, device = "png")
####Do the same for the words appearing more than once####
wcm1 <- read.csv("Data/word_counts_morethanone/wordcounts_all.csv",stringsAsFactors = FALSE)
wcm1$X <-rep(c(1:34),each=32)
plot(wcm1$NumberOfWords~ wcm1$X)
###
nwordsm1 <- subset(wcm1, wcm1$file != "contextmesicLsystem_edited_cells_NotConverge.txt" &
                   wcm1$file != "contextxericLsystem_edited_cells_NotConverge.txt" &
                   wcm1$file != "Ray_Lsystem_edited_cells_NotConverge.txt")
matcher <- match(nwordsm1$file, species.id$files)
nwordsm1$species <- species.id$species[matcher]
nwordsm1$habit <- species.id$habit[matcher]
#write.csv(nwords,"Data/wc_resumen.csv")#define the correct folder
plot(nwordsm1$X,nwordsm1$NumberOfWords)
plot(log10(nwordsm1$X), log10(nwordsm1$NumberOfWords))
nwordslog <- ggplot(data= nwordsm1) + 
  aes(x= X, y= log10(NumberOfWords),group = file,color= habit) + geom_line(aes(linetype=habit),size = 1) +
  scale_color_manual(values = c("#CB2314","#FAD510","#273046"))
nwordsnormal <- ggplot(data= nwordsm1) + 
  aes(x= X, y= NumberOfWords, group = file, color= habit) +  geom_line(aes(linetype=habit),size = 1) +
  scale_color_manual(values = c("#CB2314","#FAD510","#273046"))
####
ggsave("Figures/nwordsmoreone_log.pdf",nwordslog, device = "pdf")
ggsave("Figures/nwordsmoreone.pdf",nwordsnormal, device = "pdf")
ggsave("Figures/nwordsmoreone_log.png",nwordslog, device = "png")
ggsave("Figures/nwordsmoreone.png",nwordsnormal, device = "png")
####Table to see max number of words
maxwords <- nwords %>% group_by(file) %>% top_n(1, NumberOfWords)
maxwords <- maxwords[order(maxwords$file, decreasing = TRUE), ]      # Order data
maxwords <- maxwords[!duplicated(maxwords$file), ]   # Unique rows of ordered data
maxwordsm1 <- nwordsm1 %>% group_by(file) %>% top_n(1, NumberOfWords)
maxwordsm1 <- maxwordsm1[order(maxwordsm1$file, decreasing = TRUE), ]      # Order data
maxwordsm1 <- maxwordsm1[!duplicated(maxwordsm1$file), ]   # Unique rows of ordered data
#Subset l-systems
boxplot(maxwords$X~maxwords$habit)
boxplot(maxwordsm1$X~maxwordsm1$habit)
plot(log10(maxwordsm1$NumberOfWords) ~ maxwordsm1$X)
mean(maxwords$X)
mean(maxwordsm1$X)
#
maxwordsm1<-subset(maxwordsm1, 
       maxwordsm1$file != "probLsystembeta_edited_cells_NotConverge.txt" &
      maxwordsm1$file != "probLsystem_edited_cells_NotConverge.txt")
summary(maxwordsm1$X)
####Check mean coded row length ####
cell_lengths <- read.csv("Data/cell_lengths_notConverge.csv")
cell_lengths_sinR <- read.csv("Data/cell_lengths_withoutR.csv")
#
cell_lengths_mean <- aggregate(Number.of.cells ~ Sample, data = cell_lengths, 
                          FUN = mean, na.rm = TRUE)
cell_lengths_sinR_mean <- aggregate(Number.of.cells ~ Sample, data = cell_lengths_sinR, 
                               FUN = mean, na.rm = TRUE)
cell_lengths_sum <- aggregate(Number.of.cells ~ Sample, data = cell_lengths, 
                               FUN = sum, na.rm = TRUE)
cell_lengths_sinR_sum <- aggregate(Number.of.cells ~ Sample, data = cell_lengths_sinR, 
                                    FUN = sum, na.rm = TRUE)
cell_lengths<-cbind(cell_lengths_mean,cell_lengths_sinR_mean$Number.of.cells,
      cell_lengths_sum$Number.of.cells, cell_lengths_sinR_sum$Number.of.cells)
colnames(cell_lengths)[c(2,3,4,5)]<- c("MeanCellLength","MeanCellLengthWithoutRays",
                                       "TotalCodedCells","TotalCodedCellsWithoutRays")
#load stem diameter data
stemdiameter <- read.csv("meta/samples.csv")
#Load total of coded cells
#Graficar mean cell length contra stem diameter y #of words
matcher<-match(maxwords$file,cell_lengths$Sample)
maxwords$MeanCellLength <- cell_lengths$MeanCellLength[matcher]
maxwords$MeanCellLengthWithoutRays <- cell_lengths$MeanCellLengthWithoutRays[matcher]
maxwords$TotalCodedCellsWithoutRays<- cell_lengths$TotalCodedCellsWithoutRays[matcher]
#Add to m1 
matcher <- match(maxwordsm1$file,cell_lengths$Sample)
maxwordsm1$MeanCellLength <- cell_lengths$MeanCellLength[matcher]
maxwordsm1$MeanCellLengthWithoutRays <- cell_lengths$MeanCellLengthWithoutRays[matcher]
maxwordsm1$TotalCodedCellsWithoutRays <- cell_lengths$TotalCodedCellsWithoutRays[matcher]
#Ajustar un modelo 
hist(log10(maxwordsm1$NumberOfWords))
hist(log10(maxwords$MeanCellLengthWithoutRays))
#Delete l-systems to adjust model
withoutLsystems <- subset(maxwordsm1,
                          maxwordsm1$file != "probLsystembeta_edited_cells_NotConverge.txt" &
                          maxwordsm1$file != "probLsystem_edited_cells_NotConverge.txt")
lm.maxwords_meancelllength <- lm(log10(withoutLsystems$NumberOfWords)~ 
                                   log10(withoutLsystems$MeanCellLengthWithoutRays))
summary(lm.maxwords_meancelllength)                        
lm.maxwords_meancelllengthhabit <- lm(log10(withoutLsystems$NumberOfWords)~ 
                                   log10(withoutLsystems$MeanCellLengthWithoutRays)+
                                   withoutLsystems$habit)
summary(lm.maxwords_meancelllengthhabit)
interceptmesic <- lm.maxwords_meancelllengthhabit$coefficients[1]
interceptxeric <- lm.maxwords_meancelllengthhabit$coefficients[1]+lm.maxwords_meancelllengthhabit$coefficients[3]
slope <- lm.maxwords_meancelllengthhabit$coefficients[2]
pdf("Figures/meancellengthnumberwordssinR.pdf", width = 7, height = 5)
plot(maxwordsm1$NumberOfWords~ maxwordsm1$MeanCellLengthWithoutRays,log="xy",
     xlab="Mean cell length (number of cells in row)",ylab="Total number of words")
points(maxwordsm1$NumberOfWords[maxwordsm1$habit=="xeric"] ~ 
       maxwordsm1$MeanCellLengthWithoutRays[maxwordsm1$habit=="xeric"],col="#273046",pch=16)
points(maxwordsm1$NumberOfWords[maxwords$habit=="mesic"] ~ 
         maxwordsm1$MeanCellLengthWithoutRays[maxwordsm1$habit=="mesic"],col="#FAD510",pch=22)
points(maxwordsm1$NumberOfWords[maxwords$habit=="L-system"] ~ 
         maxwordsm1$MeanCellLengthWithoutRays[maxwordsm1$habit=="L-system"],col="#CB2314",pch=16)
abline(interceptmesic,slope,col="#FAD510")
abline(interceptxeric,slope,col="#273046")
dev.off()
lm.maxwords_totalcodedcells <- lm(log10(withoutLsystems$NumberOfWords)~ 
                                 log10(withoutLsystems$TotalCodedCellsWithoutRays))
summary(lm.maxwords_totalcodedcells)                        
lm.maxwords_totalcodedcellshabit <- lm(log10(withoutLsystems$NumberOfWords)~ 
                                    log10(withoutLsystems$TotalCodedCellsWithoutRays)*
                                    withoutLsystems$habit)
summary(lm.maxwords_totalcodedcellshabit)                        
lm.maxwords_totalcodedcellshabit <- lm(log10(withoutLsystems$NumberOfWords)~ 
                                         log10(withoutLsystems$TotalCodedCellsWithoutRays)+
                                         withoutLsystems$habit)
summary(lm.maxwords_totalcodedcellshabit)                        
confint(lm.maxwords_meancelllengthhabit)
interceptxeric<-lm.maxwords_totalcodedcellshabit$coefficients[1]+lm.maxwords_totalcodedcellshabit$coefficients[3]
interceptmesic<-lm.maxwords_totalcodedcellshabit$coefficients[1] 
slope <-lm.maxwords_totalcodedcellshabit$coefficients[2] 
pdf("Figures/codedcellsnumberwordssinR.pdf", width = 7, height = 5)
plot(maxwordsm1$NumberOfWords~ maxwordsm1$TotalCodedCellsWithoutRays,log="xy",
     xlab="Total Wood cells coded", ylab = "Total number of words")
points(maxwordsm1$NumberOfWords[maxwordsm1$habit=="xeric"] ~ 
         maxwordsm1$TotalCodedCellsWithoutRays[maxwordsm1$habit=="xeric"],col="#273046",pch=16)
points(maxwordsm1$NumberOfWords[maxwordsm1$habit=="mesic"] ~ 
         maxwordsm1$TotalCodedCellsWithoutRays[maxwordsm1$habit=="mesic"],col="#FAD510",pch=22)
points(maxwordsm1$NumberOfWords[maxwordsm1$habit=="L-system"] ~ 
         maxwordsm1$TotalCodedCellsWithoutRays[maxwordsm1$habit=="L-system"],col="#CB2314",pch=16)
abline(interceptmesic,slope,col="#FAD510")
abline(interceptxeric,slope,col="#273046")
dev.off()
#
plot(maxwords$X~ maxwords$MeanCellLengthWithoutRays,log="xy")
#Also check stem diameter
#Problem: there are some measures in one variable and other in more than one....
stemdiameter$stemdiameter.xylem
#
insert <- "_edited_cells_NotConverge.txt"
stemdiameter$filesconverge <- paste0(stemdiameter$Sample,insert)

matcher <-match(maxwords$file, stemdiameter$filesconverge)
maxwords$xylemdiameter <- stemdiameter$stemdiameter.xylem[matcher]
maxwords$pithdiameter <- stemdiameter$pithdiameter[matcher]
#maxwords$stemdiameter <- cell_lengths$MeanCellLengthWithoutRays[matcher]
#maxwords$TotalCodedCellsWithoutRays<- cell_lengths$TotalCodedCellsWithoutRays[matcher]
plot(maxwords$NumberOfWords ~ maxwords$xylemdiameter)
#
plot(maxwords$NumberOfWords ~ maxwords$pithdiameter)

#Make loop to get the total possible words at different lengths:
nwords$universe <- 3**nwords$X
nwords$devpot<-nwords$NumberOfWords/nwords$universe
plot(nwords$devpot ~nwords$X)
boxplot(nwords$devpot~nwords$species)
boxplot(nwords$devpot~nwords$habit)
#
ggplot(data= nwords, aes(x= X, y=devpot, group = file, 
       color= habit)) + geom_line(size = 1)

#####Using wc.n files to make a matrix and determine ######
# euclidean or bray-kurtis distances #####Only wordcounts at 2, 3, 7, 19 
wc.2 <- as.matrix(read.csv("Data/word_counts_morethanone/wordcounts2.csv",row.names=1))
wc.3 <- as.matrix(read.csv("Data/word_counts_morethanone/wordcounts3.csv", row.names = 1))
wc.3[is.na(wc.3)] <- 0
wc.7 <- as.matrix(read.csv("Data/word_counts_morethanone/wordcounts7.csv", row.names = 1))
wc.7[is.na(wc.7)] <- 0
wc.19 <- as.matrix(read.csv("Data/word_counts_morethanone/wordcounts19.csv", row.names = 1))
wc.19[is.na(wc.19)] <- 0
#
euc.dist.r.2 <- dist(wc.2)
euc.dist.r.7 <- dist(wc.7)
#Get proportion of counts
wc.2.freq <- wc.2/rowSums(wc.2)
euc.dist.r2.freq <- dist(wc.2.freq)
plot(hclust(euc.dist.r2.freq))
#
#Euclidean distances calculated with r are the same as the ones in python
#Try to solve problem of different total of cells
boxplot(rowSums(wc.2))
hist(rowSums(wc.3))
hist(rowSums(wc.7))
#Try to count cell number for each and calc the counts per 30,000.Remove L-systems
#remove l-systems
remove <- rownames(wc.2)[c(16,25,30)]
wc.2.sinL<- wc.2[!rownames(wc.2) %in% remove, ]  # ! is logical negation
porcadamiles <- rowSums(wc.2.sinL)/30000
wc.2.portreintamil <- round((wc.2.sinL/porcadamiles),0)
#wc.2.portreintamil.freq <- wc.2.portreintamil/rowSums(wc.2.portreintamil)
euc.dist.2.std <- dist(wc.2.portreintamil)
plot(hclust(euc.dist.2.std))
####para 3
remove <- rownames(wc.3)[c(16,19,21,25,30)]
wc.3.sinL<- wc.3[!rownames(wc.3) %in% remove, ]  # ! is logical negation
porcadamiles <- rowSums(wc.3.sinL)/29000
wc.3.portreintamil <- round((wc.3.sinL/porcadamiles),0)
#wc.2.portreintamil.freq <- wc.2.portreintamil/rowSums(wc.2.portreintamil)
euc.dist.3.std <- dist(wc.3.portreintamil)
plot(hclust(euc.dist.3.std))
#
###for 7
remove <- rownames(wc.7)[c(18,20,24,29,32)]
wc.7.sinL<- wc.7[!rownames(wc.7) %in% remove, ]  # ! is logical negation
#Determine sample with more word counted. 
head(sort(rowSums(wc.7.sinL),decreasing = TRUE))
porcadamiles <- rowSums(wc.7.sinL)/28000
wc.7.portreintamil <- round((wc.7.sinL/porcadamiles),0)
euc.dist.7.std <- dist(wc.7.portreintamil)
plot(hclust(euc.dist.7.std))
#Calc the bray-curtis distance and create heatmap with thosee values
bray.dist.7 <- vegdist(wc.7.portreintamil, method="bray")
#Applying hclust to the not normalized
#######NMDS
ord7 <- metaMDS(wc.7.portreintamil)#
species.id_onlysp <-subset(species.id, species.id$habit != "L-system"  & species.id$habit != "RayL-system" &
                             species.id$habit != "L-systemCSM" & species.id$habit != "L-systemCSX")
# Match the id of files from datafrmae with more info with the id files 
# of the matrix counts
matcher<-match(rownames(wc.7.portreintamil), species.id_onlysp$files)
ordiplot(ord7, type = "n", main = "ellipses")
#Plot the if of the files distingushing between habits.
orditorp(ord7, display = "sites", labels = F, 
         pch = c(16, 18) [as.numeric(as.factor(species.id_onlysp$habit[matcher]))], 
         col = c("#FAD510","#273046") [as.numeric(as.factor(species.id_onlysp$habit[matcher]))], 
         cex = 1)
#Test to get the significant words!
# Check words that are driving the speccies distribution pattern
ord7.dist.fit <- envfit(ord7, wc.7.portreintamil, permutations = 9999)
head(ord7.dist.fit)
pdf("Figures/NDMS_Kmer7.pdf", width = 7, height = 5)
ordiplot(ord7, type = "n", main = "ellipses")
#Plot the if of the files distingushing between habits.
plot(ord7.dist.fit, p.max = 0.0001, col = alpha("black",0.6), cex=0.6)  
orditorp(ord7, display = "sites", labels = F, 
         pch = c(16, 18) [as.numeric(as.factor(species.id_onlysp$habit[matcher]))], 
         col = c("#FAD510","#273046") [as.numeric(as.factor(species.id_onlysp$habit[matcher]))], 
         cex = 1)
dev.off()
#
ordiellipse(ord7, groups = as.factor(species.id$habit[matcher]), 
            draw = "polygon", lty = 1, col = "grey90")
#
par(mfrow = c(1, 2))
stressplot(ord7, main = "Shepard plot")
gof <- goodness(ord7)
plot(ord7, type = "t", main = "Goodness of fit")
points(ord7, display = "sites", cex = gof * 300)
#######Determine counts based on growth from############################
#First do it with word counts not scaled
wc7.df<-as.data.frame(cbind(rownames(wc.7.sinL),wc.7.sinL))#transform matrix to dataframe
colnames(wc7.df)[1]<- "Sample" 
dim(wc7.df) #check number of columns
wc7.df <- wc7.df %>% pivot_longer(cols=c(2:2188),
                                    names_to='Words',  
                                    values_to='Numberofcounts') #Make dataframe in long format
wc7.df$Numberofcounts <- as.numeric(wc7.df$Numberofcounts) #transform counts to numeric
matcher <- match(wc7.df$Sample, species.id$files) #add species info based on species.id
wc7.df$species <- species.id$species[matcher]
wc7.df$habit <- species.id$habit[matcher]
wc7.df.species<-aggregate(Numberofcounts ~ 
                            species+Words, wc7.df, sum) #summarize counts at species level
matcher <- match(wc7.df.species$species, species.id$species) #add growth form to df
wc7.df.species$habit <- species.id$habit[matcher]
#
wc7.df.habit <- aggregate(Numberofcounts ~ 
                            habit+Words, wc7.df.species, sum) #summarize based on growth form
sum(wc7.df.habit$Numberofcounts==0) #count words with zero
sum(wc7.df.habit$Numberofcounts==0 & wc7.df.habit$habit=="xeric") #count words with 0s and habit
sum(wc7.df.habit$Numberofcounts==0 & wc7.df.habit$habit=="mesic")
####Do it with counts scaled
wc7.df.scaled<-as.data.frame(cbind(rownames(wc.7.portreintamil),wc.7.portreintamil))
colnames(wc7.df.scaled)[1]<- "Sample"
dim(wc7.df.scaled)
wc7.df.scaled <- wc7.df.scaled %>% pivot_longer(cols=c(2:2188),
                                  names_to='Words',
                                  values_to='Numberofcounts')
wc7.df.scaled$Numberofcounts <- as.numeric(wc7.df.scaled$Numberofcounts)
#
matcher <- match(wc7.df.scaled$Sample, species.id$files)
wc7.df.scaled$species <- species.id$species[matcher]
wc7.df.scaled.species <- aggregate(Numberofcounts ~ species+Words, 
                                   wc7.df.scaled, sum)
matcher <- match(wc7.df.scaled.species$species, species.id$species)
wc7.df.scaled.species$habit <- species.id$habit[matcher]
#
wc7.df.scaled.habit<-aggregate(Numberofcounts ~ habit+Words, wc7.df.scaled.species, sum)
sum(wc7.df.scaled.habit$Numberofcounts<10)
sum(wc7.df.scaled.habit$Numberofcounts<10 & wc7.df.scaled.habit$habit=="xeric")
sum(wc7.df.scaled.habit$Numberofcounts<10 & wc7.df.scaled.habit$habit=="mesic")

######For 19
remove <- rownames(wc.19)[c(22,27,30,31)]
wc.19.sinL<- wc.19[!rownames(wc.19) %in% remove, ]  # ! is logical negation
head(sort(rowSums(wc.19.sinL),decreasing = TRUE))
porcadamiles <- rowSums(wc.19.sinL)/17000
wc.19.porcadadiecisiete <- round((wc.19.sinL/porcadamiles),0)
euc.dist.19.std <- dist(wc.19.porcadadiecisiete)
plot(hclust(euc.dist.19.std))
bray.dist.19 <- vegdist(wc.19.porcadadiecisiete, method="bray")
#Create ndms and plot
ord19 <- metaMDS(wc.19.porcadadiecisiete)#
matcher<-match(rownames(wc.19.porcadadiecisiete), species.id_onlysp$files)
#
ord19.dist.fit <- envfit(ord19, wc.19.porcadadiecisiete, permutations = 9999)
head(ord19.dist.fit)
saveRDS(ord19.dist.fit, file="../Meta/ord19.dist.envfit.RData")

pdf("Figures/NDMS_Kmer19.pdf", width = 7, height = 5)
ordiplot(ord19, type = "n", main = "ellipses")
plot(ord19.dist.fit, p.max = 0.0001, col = alpha("black",0.4), cex=0.6)  
orditorp(ord19, display = "sites", labels = F, 
         pch = c(16, 18) [as.numeric(as.factor(species.id_onlysp$habit[matcher]))], 
         col = c("#FAD510","#273046") [as.numeric(as.factor(species.id_onlysp$habit[matcher]))], 
         cex = 1)
dev.off()
ordiellipse(ord19, groups = as.factor(species.id$habit[matcher]), 
            draw = "polygon", lty = 1, col = "grey90")
#
stressplot(ord19, main = "Shepard plot")
gof <- goodness(ord19)
plot(ord19, type = "t", main = "Goodness of fit")
points(ord, display = "sites", cex = gof * 300)
#####
pheatmap2 <- pheatmap(as.matrix(euc.dist.2.std), display_numbers = F)
pdf("Figures/heatmaplengthStandarized_2.pdf") # Para guardar en PDF
pheatmap2
dev.off()
#3
pheatmap3 <- pheatmap(as.matrix(euc.dist.3.std), display_numbers = F)
pdf("Figures/heatmaplengthStandarized_3.pdf") # Para guardar en PDF
pheatmap3
dev.off()
#7
pheatmap7 <- pheatmap(as.matrix(euc.dist.7.std), display_numbers = F)
pdf("Figures/heatmaplengthStandarized_7.pdf") # Para guardar en PDF
pheatmap7
dev.off()
#
pheatmap7 <- pheatmap(as.matrix(bray.dist.7), display_numbers = F)
pdf("Figures/heatmaplengthbray_7.pdf") # Para guardar en PDF
pheatmap7
dev.off()
#
pheatmap19 <- pheatmap(as.matrix(euc.dist.19.std), display_numbers = F)
pdf("Figures/heatmaplengthStandarized_19.pdf") # Para guardar en PDF
pheatmap19
dev.off()
#
pheatmap19 <- pheatmap(as.matrix(bray.dist.19), display_numbers = F)
pdf("Figures/heatmaplengthbray_19.pdf") # Para guardar en PDF
pheatmap19
dev.off()
####
wc.7.prueba <- vegdist(wc.7.portreintamil)

wc.7.pcoa <- cmdscale(wc.7.prueba, k = (nrow(wc.7.portreintamil) - 1), eig = TRUE)
# Plot of the sites
ordiplot(scores(wc.7.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with species weighted averages")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
#####
wc.19.prueba <- vegdist(wc.19.porcadadiecisiete)

wc.19.pcoa <- cmdscale(wc.19.prueba, k = (nrow(wc.19.porcadadiecisiete) - 1), eig = TRUE)
# Plot of the sites
ordiplot(scores(wc.19.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with species weighted averages")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

########Make a tree
library(phangorn)
tree2<-upgma(euc.dist.2.std)
plot(tree2)
tree7<-upgma(euc.dist.7.std)
plot(tree7)
tree19<-upgma(euc.dist.19.std)
plot(tree19)
####
## Para guardar en Png
write.tree(tree7,"Figures/upgmaStandarized7")
png("Figures/upgmaStandarized_7.png")
plot(tree7)
dev.off()
write.tree(tree19,"Figures/upgmaStandarized19")
png("Figures/upgmaStandarized_19.png") # Para guardar en PDF
plot(tree19)
dev.off()
#################################################################
###Make the same for word counts with appearing also one time###
summary(maxwords$X)
wc.12 <- as.matrix(read.csv("Data/word_counts_all/wordcounts12.csv",row.names=1))
wc.12[is.na(wc.12)] <- 0
#wc.30 <- as.matrix(read.csv("Data/word_counts_all/wordcounts30.csv", row.names = 1))
#wc.30[is.na(wc.30)] <- 0
#
euc.dist.r.12 <- dist(wc.12)
#Get proportion of counts
wc.12.freq <- wc.12/rowSums(wc.12)
euc.dist.r12.freq <- dist(wc.12.freq)
plot(hclust(euc.dist.r12.freq))
#
#Euclidean distances calculated with r are the same as the ones in python
#Try to solve problem of different total of cells
boxplot(rowSums(wc.12))
hist(rowSums(wc.12))
#Try to count cell number for each and calc the counts per 30,000.Remove L-systems
#remove l-systems
remove <- rownames(wc.12)[c(18,23,28,31,32)]
wc.12.sinL<- wc.12[!rownames(wc.12) %in% remove, ]  # ! is logical negation
head(sort(rowSums(wc.12.sinL),decreasing = TRUE))
porcadamiles <- rowSums(wc.12.sinL)/27292
wc.12.portreintamil <- round((wc.12.sinL/porcadamiles),0)
euc.dist.12.std <- dist(wc.12.portreintamil)
plot(hclust(euc.dist.12.std))
#####
#heatmap.2(as.matrix(euc.dist.r.2),key=TRUE,scale = "row",
#            margins = c(10,12), cexRow=0.5)
pheatmap12 <- pheatmap(as.matrix(euc.dist.12.std), display_numbers = F)
pdf("Figures/heatmaplengthStandarized_12.pdf") # Para guardar en PDF
pheatmap12
dev.off()
