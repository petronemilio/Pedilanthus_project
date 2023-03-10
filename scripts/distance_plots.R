##### Distance measures #####
# In this script, plots and tables are made
# to analize different measures of word diversity
# and complexity
pal1<-c("#a5c533", "#9658d8", "#59d261", "#9c3db2", "#75bb38", "#5f77f3", "#bebe35", "#3659cb",
        "#40af43", "#ca37a3", "#39c477", "#ee68d2", "#4f901e", "#cf6fe4", "#8fa837", "#8163d6",
        "#e6ac35", "#614bae", "#7dbe68", "#a33791", "#4dc38d", "#d32f80", "#49852d", "#ac84e6",
        "#bca53a", "#507fe3", "#e47f2e", "#4a9de3", "#e55632", "#3ec7da", "#b82e25", "#59ceb8",
        "#de3953", "#308949", "#e861ac", "#3d733d", "#ea93e1", "#606b16", "#82479d", "#90be7e",
        "#6b4fa0", "#cb8e35", "#345fa9", "#b3511f", "#58b8e3", "#c33361", "#54a27a", "#eb658b",
        "#287e63", "#c972bd", "#7c9451", "#5259a9", "#928432", "#8897ec", "#7e5f16", "#8775be",
        "#bcb66e", "#a15195", "#1a6447", "#b04a79", "#33a29e", "#c24f4e", "#2e73a9", "#eb8665",
        "#5e99c8", "#af6938", "#788acb", "#8e4b27", "#c0a8e8", "#5b672f", "#775696", "#d69e6c",
        "#576196", "#8f6f3d", "#b17db7", "#9e4d51", "#e18db3", "#92465f", "#dc8484", "#8b4c76")
##### Load useful libraries #####
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(corrplot)
library(gplots)
library(wesanderson)
library(pheatmap)
##### Load files
## Loading species id file
species.id <- read.csv("meta/speciesID.csv")
insert <- "_NotConverge"
species.id$filesconverge <- sub("(?<=\\w)(?=\\.)", insert, species.id$files, perl=TRUE)
## => A1.A1c<=7 A1.A1c<=8
#
euc.dist2<-read.csv("euclidean_distance2.csv")
euc.dist2[c("File1", "File2")] <- lapply(euc.dist2[c("File1", "File2")], factor,
                                         levels=unique(unlist(euc.dist2[c("File1", "File2")])))
#####Transform to matrix format#####
matrix_euc2 <- xtabs(Euc_dist ~ File1+ File2,euc.dist2)
corrplot(matrix_euc2, is.corr = FALSE,tl.pos='n')
##
matrix_euc3 <- xtabs(Euc_dist ~ File1+ File2,euc_dist3)
corrplot(matrix_euc3, is.corr = FALSE,type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05,
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=0.5, outline=T)
#####Load word counts to plot as a function of word length #####
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
write.csv(nwords,"Data/wc_resumen.csv")#define the correct folder
plot(nwords$X,nwords$NumberOfWords)
nwordslog <- ggplot(data= nwords,aes(x= X, y= log10(NumberOfWords),group = file, 
                                     color= habit)) + geom_line(size = 1)
nwordsnormal <- ggplot(data= nwords, aes(x= X, y= NumberOfWords, group = file, color= habit)) +
  geom_line(size = 1)
nwordsnormal
nwordslog
ggsave("Figures/nwords_log.pdf",nwordslog, device = "pdf")
ggsave("Figures/nwords_notlog.pdf",nwordsnormal, device = "pdf")
ggsave("Figures/nwords_log.png",nwordslog, device = "png")
ggsave("Figures/nwords_notlog.png",nwordsnormal, device = "png")
#
#remove specific factors from plot
nwords.temp <- subset(nwords, nwords$habit != "probL-system" &
                        nwords$habit != "probetaL-system")
#
nwords_sp <- ggplot(data= nwords.temp, aes(x= X, y= NumberOfWords, group = file, 
                                           color= species)) + geom_line(size = 1)
nwords_sphabit <- ggplot(data= nwords.temp,aes(x= X, y= NumberOfWords, group = file, 
                         color= habit)) + geom_line(size = 1)
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
nwordslog <- ggplot(data= nwordsm1,aes(x= X, y= log10(NumberOfWords),group = file, 
                                     color= habit)) + geom_line(size = 1)
nwordsnormal <- ggplot(data= nwordsm1, aes(x= X, y= NumberOfWords, group = file, color= habit)) +
  geom_line(size = 1)
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

boxplot(maxwords$X~maxwords$habit)
boxplot(maxwordsm1$X~maxwordsm1$habit)
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
       maxwordsm1$MeanCellLengthWithoutRays[maxwordsm1$habit=="xeric"],col="red",pch=16)
points(maxwordsm1$NumberOfWords[maxwords$habit=="mesic"] ~ 
         maxwordsm1$MeanCellLengthWithoutRays[maxwordsm1$habit=="mesic"],col="blue",pch=16)
abline(interceptmesic,slope,col="blue")
abline(interceptxeric,slope,col="red")
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
interceptxeric<-lm.maxwords_totalcodedcellshabit$coefficients[1]+lm.maxwords_totalcodedcellshabit$coefficients[3]
interceptmesic<-lm.maxwords_totalcodedcellshabit$coefficients[1] 
slope <-lm.maxwords_totalcodedcellshabit$coefficients[2] 
pdf("Figures/codedcellsnumberwordssinR.pdf", width = 7, height = 5)
plot(maxwordsm1$NumberOfWords~ maxwordsm1$TotalCodedCellsWithoutRays,log="xy",
     xlab="Total Wood cells coded", ylab = "Total number of words")
points(maxwordsm1$NumberOfWords[maxwordsm1$habit=="xeric"] ~ 
         maxwordsm1$TotalCodedCellsWithoutRays[maxwordsm1$habit=="xeric"],col="red",pch=16)
points(maxwordsm1$NumberOfWords[maxwordsm1$habit=="mesic"] ~ 
         maxwordsm1$TotalCodedCellsWithoutRays[maxwordsm1$habit=="mesic"],col="blue",pch=16)
abline(interceptmesic,slope,col="blue")
abline(interceptxeric,slope,col="red")
dev.off()
#
plot(maxwords$X~ maxwords$MeanCellLengthWithoutRays,log="xy")
#Also check stem diameter

#Make loop to get the total possible words at different lengths:
nwords$universe <- 3**nwords$X
nwords$devpot<-nwords$NumberOfWords/nwords$universe
plot(nwords$devpot ~nwords$X)
boxplot(nwords$devpot~nwords$species)
boxplot(nwords$devpot~nwords$habit)
#
ggplot(data= nwords, aes(x= X, y=devpot, group = file, 
       color= habit)) + geom_line(size = 1)

#####Trying to use wc.n files to make a matrix and determine ######
# euclidean distances #####
wc.2 <- as.matrix(read.csv("Data/word_counts_morethanone/wordcounts2.csv",row.names=1))
wc.3 <- as.matrix(read.csv("Data/word_counts_morethanone/wordcounts3.csv", row.names = 1))
wc.3[is.na(wc.3)] <- 0
wc.5 <- as.matrix(read.csv("Data/word_counts_morethanone/wordcounts5.csv", row.names = 1))
wc.5[is.na(wc.5)] <- 0
wc.7 <- as.matrix(read.csv("Data/word_counts_morethanone/wordcounts7.csv", row.names = 1))
wc.7[is.na(wc.7)] <- 0
wc.12 <- as.matrix(read.csv("Data/word_counts_morethanone/wordcounts12.csv", row.names = 1))
wc.12[is.na(wc.12)] <- 0
wc.15 <- as.matrix(read.csv("Data/word_counts_morethanone/wordcounts15.csv", row.names = 1))
wc.15[is.na(wc.15)] <- 0
wc.18 <- as.matrix(read.csv("Data/word_counts_morethanone/wordcounts18.csv", row.names = 1))
wc.18[is.na(wc.18)] <- 0
wc.19 <- as.matrix(read.csv("Data/word_counts_morethanone/wordcounts19.csv", row.names = 1))
wc.19[is.na(wc.19)] <- 0
wc.22 <- as.matrix(read.csv("Data/word_counts_morethanone/wordcounts22.csv", row.names = 1))
wc.22[is.na(wc.22)] <- 0
#
euc.dist.r.2 <- dist(wc.2)
euc.dist.r.8 <- dist(wc.8)
#Get proportion of counts
wc.2.freq <- wc.2/rowSums(wc.2)
euc.dist.r2.freq <- dist(wc.2.freq)
plot(hclust(euc.dist.r2.freq))
#
#Euclidean distances calculated with r are the same as the ones in python
#Try to solve problem of different total of cells
boxplot(rowSums(wc.2))
hist(rowSums(wc.3))
hist(rowSums(wc.8))
rowSums(wc.12)
#Try to count cell number for each and calc the counts per 30,000.Remove L-systems
#remove l-systems
remove <- rownames(wc.2)[c(16,19,21,25,30)]
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
####para 5
remove <- rownames(wc.5)[c(16,19,21,25,30)]
wc.5.sinL<- wc.5[!rownames(wc.5) %in% remove, ]  # ! is logical negation
head(sort(rowSums(wc.5.sinL),decreasing = TRUE))
porcadamiles <- rowSums(wc.5.sinL)/29000
wc.5.porveintinuevemil <- round((wc.5.sinL/porcadamiles),0)
#wc.2.portreintamil.freq <- wc.2.portreintamil/rowSums(wc.2.portreintamil)
euc.dist.5.std <- dist(wc.5.porveintinuevemil)
plot(hclust(euc.dist.5.std))
###para 8
remove <- rownames(wc.7)[c(18,20,24,29,32)]
wc.7.sinL<- wc.7[!rownames(wc.7) %in% remove, ]  # ! is logical negation
#Determine sample with more word counted. 
head(sort(rowSums(wc.7.sinL),decreasing = TRUE))
porcadamiles <- rowSums(wc.7.sinL)/28000
wc.7.portreintamil <- round((wc.7.sinL/porcadamiles),0)
#wc.2.portreintamil.freq <- wc.2.portreintamil/rowSums(wc.2.portreintamil)
euc.dist.7.std <- dist(wc.7.portreintamil)
plot(hclust(euc.dist.7.std))
#para 15
remove <- rownames(wc.15)[c(22,27,30:32)]
wc.15.sinL<- wc.15[!rownames(wc.15) %in% remove, ]  # ! is logical negation
head(sort(rowSums(wc.15.sinL),decreasing = TRUE))
porcadamiles <- rowSums(wc.15.sinL)/23000
wc.15.porveintitresmil <- round((wc.15.sinL/porcadamiles),0)
#wc.2.portreintamil.freq <- wc.2.portreintamil/rowSums(wc.2.portreintamil)
euc.dist.15.std <- dist(wc.15.porveintitresmil)
plot(hclust(euc.dist.15.std))
#Para 19
remove <- rownames(wc.19)[c(22,27,30,31)]
wc.19.sinL<- wc.19[!rownames(wc.19) %in% remove, ]  # ! is logical negation
head(sort(rowSums(wc.19.sinL),decreasing = TRUE))
porcadamiles <- rowSums(wc.19.sinL)/17000
wc.19.porcadadiecisiete <- round((wc.19.sinL/porcadamiles),0)
#wc.2.portreintamil.freq <- wc.2.portreintamil/rowSums(wc.2.portreintamil)
euc.dist.19.std <- dist(wc.19.porcadadiecisiete)
plot(hclust(euc.dist.19.std))
#Para 22
remove <- rownames(wc.22)[c(22,27,30)]
wc.22.sinL<- wc.22[!rownames(wc.22) %in% remove, ]  # ! is logical negation
head(sort(rowSums(wc.22.sinL),decreasing = TRUE))
porcadamiles <- rowSums(wc.22.sinL)/12000
wc.22.porcadadocemil <- round((wc.22.sinL/porcadamiles),0)
#wc.2.portreintamil.freq <- wc.2.portreintamil/rowSums(wc.2.portreintamil)
euc.dist.22.std <- dist(wc.22.porcadadocemil)
plot(hclust(euc.dist.22.std))
####
#heatmap.2(as.matrix(euc.dist.r.2),key=TRUE,scale = "row",
#            margins = c(10,12), cexRow=0.5)
pheatmap(as.matrix(euc.dist.2.std), display_numbers = F)
pdf("Figures/heatmaplengthStandarized_2.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.2.std), display_numbers = F)
dev.off()
#3
pdf("Figures/heatmaplengthStandarized_3.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.3.std), display_numbers = F)
dev.off()
#5
pdf("Figures/heatmaplengthStandarized_5.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.5.std), display_numbers = F)
dev.off()
#
pdf("Figures/heatmaplengthStandarized_7.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.7.std), display_numbers = F)
dev.off()
#
pdf("Figures/heatmaplengthStandarized_15.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.15.std), display_numbers = F)
dev.off()
#
pdf("Figures/heatmaplengthStandarized_19.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.19.std), display_numbers = F)
dev.off()
#
pdf("Figures/heatmaplengthStandarized_22.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.22.std), display_numbers = F)
dev.off()
#
matcher <- match(wc.2$sample, species.id$files)
wc.2$habit <- species.id$habit[matcher]
#3
matcher <- match(wc.3$sample, species.id$files)
wc.3$habit <- species.id$habit[matcher]
#4
matcher <- match(wc.4$sample, species.id$files)
wc.4$habit <- species.id$habit[matcher]
######Determine just 
wc.2 <- subset(wc.2, wc.2$value != 0)
wc.2.gf <- as.data.frame(table(wc.2$habit, wc.2$variable))
wc.2.gf <- subset(wc.2.gf, wc.2.gf$Freq != 0)
nwords2 <- wc.2.gf %>% group_by(Var1) %>%
  tally()
##
wc.3 <- subset(wc.3, wc.3$value != 0)
wc.3.gf <- as.data.frame(table(wc.3$habit,wc.3$variable))
wc.3.gf <- subset(wc.3.gf, wc.3.gf$Freq != 0)
nwords3 <- wc.3.gf %>% group_by(Var1) %>%
  tally()
##
wc.4 <- subset(wc.4, wc.4$value != 0)
wc.4.gf <- as.data.frame(table(wc.4$habit, wc.4$variable))
wc.4.gf <- subset(wc.4.gf, wc.4.gf$Freq != 0)
nwords4 <- wc.4.gf %>% group_by(Var1) %>%
  tally()
##
nwords <- rbind(cbind(nwords2,wordLength =rep(2,nrow(nwords2))),
                cbind(nwords3,wordLength =rep(3,nrow(nwords3))),
                cbind(nwords4,wordLength =rep(4,nrow(nwords4))))
####
m1 <-lm(nwords$n~nwords$wordLength*nwords$habit)
summary(m1)
m2 <-lm(log10(nwords$n)~ log10(nwords$wordLength) *nwords$habit)
summary(m2)

#remove l-systems
remove <- rownames(wc.10)[c(18,23,28,31,32)]
wc.sinL<- wc.10[!rownames(wc.10) %in% remove, ]  # ! is logical negation
diezporcadatreintamil <- rowSums(wc.sinL[!rownames(wc.sinL) %in% remove, ])/30000
wc.sinL <- round((wc.sinL/diezporcadatreintamil),0)
#
wc.10.portreintamil.freq <- wc.10.portreintamil/rowSums(wc.10.portreintamil)
euc.dist.r.10.freq.std <- dist(wc.10.portreintamil.freq)
plot(hclust(euc.dist.r.10.freq.std))

#####
euc.dist.r.2.std <- dist(wc.2.portreintamil)
euc.dist.r.3.std <- dist(wc.3.portreintamil)
euc.dist.r.4.std <- dist(wc.4.portreintamil)
euc.dist.r.10.std <- dist(wc.10.portreintamil)
#
euc.dist.r.10.std <- dist(wc.sinL)
library(phangorn)
tree<-upgma(euc.dist.r.10.std)
plot(tree)
####
plot(hclust(euc.dist.r.5.std))



library(wordcloud)
#pdf("Figures/WordCloudCoal886.pdf")
#wordcloud(coal886$variable,coal886$value, max.words =100,min.freq=3,scale=c(4,.5),
#          random.order = FALSE)
#dev.off()
#pdf("Figures/WordCloudDia11.pdf")
#wordcloud(dia11$variable,dia11$value, max.words =100,min.freq=3,scale=c(4,.5),
#          random.order = FALSE)
#dev.off()
#pdf("Figures/WordCloudPeri843.pdf")
#wordcloud(peri843$variable,peri843$value, max.words =100,min.freq=3,scale=c(4,.5),
#          random.order = FALSE)
#dev.off()
#pdf("Figures/WordCloudFin917.pdf")
#wordcloud(finki917$variable,finki917$value, max.words =100,min.freq=3,scale=c(4,.5),
#          random.order = FALSE)
#dev.off()
#pdf("Figures/WordCloudCym979.pdf")
#wordcloud(cym979$variable,cym979$value, max.words =100,min.freq=3,scale=c(4,.5),
#          random.order = FALSE)
#dev.off()
#pdf("Figures/WordCloudTyth2.pdf")
#wordcloud(tyth2$variable,tyth2$value, max.words =100,min.freq=3,scale=c(4,.5),
#          random.order = FALSE)
#dev.off()
#pdf("Figures/WordCloudCyr14.pdf")
#wordcloud(cyr14$variable,cyr14$value, max.words =100,min.freq=3,scale=c(4,.5),
#          random.order = FALSE)
#dev.off()
#pdf("Figures/WordCloudTeh981.pdf")
#wordcloud(teh981$variable,teh981$value, max.words =100,min.freq=3,scale=c(4,.5),
#          random.order = FALSE)
#dev.off()
#pdf("Figures/problsys.pdf")
#wordcloud(problsys$variable,problsys$value, max.words =100,min.freq=3,scale=c(4,.5),
#          random.order = FALSE)
#dev.off()
######
#plot(epm11_wc2$value ~ epm11_wc2$variable)
#plot(log10(nwords10$n)~ as.factor(nwords10$sample))
#plot(nwords5$n ~ as.factor(nwords5$sample))

######
library(vegan)

as.matrix(problsys)
div2 <- as.data.frame(diversity(wc.2.portreintamil))
div3 <- as.data.frame(diversity(wc.3.portreintamil))
div4 <- as.data.frame(diversity(wc.4.portreintamil))
div5 <- as.data.frame(diversity(wc.5.portreintamil))
div6 <- as.data.frame(diversity(wc.6.portreintamil))
div7 <- as.data.frame(diversity(wc.7.portreintamil))
div8 <- as.data.frame(diversity(wc.8.portreintamil))
div9 <- as.data.frame(diversity(wc.9.portreintamil))
div10 <- as.data.frame(diversity(wc.10.portreintamil))
diversity(wc.10.portreintamil)
?diversity
######
rownames(div2)
#
div3$sample <- rownames(div3)
matcher <- match(div3$sample,species.id$files)
div3$species <- species.id$species[matcher]
div3$habit <- species.id$habit[matcher]
boxplot(div3$`diversity(wc.3.portreintamil)` ~ div3$species)
boxplot(div3$`diversity(wc.3.portreintamil)` ~ div3$habit)
#
div5$sample <- rownames(div5)
matcher <- match(div5$sample,species.id$files)
div5$species <- species.id$species[matcher]
div5$habit <- species.id$habit[matcher]
boxplot(div5$`diversity(wc.5.portreintamil)` ~ div5$species)
boxplot(div5$`diversity(wc.5.portreintamil)` ~ div5$habit)
#
matcher <- match(rownames(div8),species.id$files)
div8$species <- species.id$species[matcher]
div8$habit <- species.id$habit[matcher]
boxplot(div8$`diversity(wc.8.portreintamil)` ~ div8$species)
boxplot(div8$`diversity(wc.8.portreintamil)` ~ div8$habit)
lm8 <- lm(div8$`diversity(wc.8.portreintamil)`~ div8$habit)
#
matcher <- match(rownames(div10),species.id$files)
div10$species <- species.id$species[matcher]
div10$habit <- species.id$habit[matcher]
boxplot(div10$`diversity(wc.10.portreintamil)` ~ div10$species)
boxplot(div10$`diversity(wc.10.portreintamil)` ~ div10$habit)
lm10 <- lm(div10$`diversity(wc.10.portreintamil)`~ div10$habit)
summary(lm10)
pdf("Figures/wordsdiversity.pdf")
par(mfrow=c(2, 2))
boxplot(div3$`diversity(wc.3.portreintamil)` ~ div3$habit)
boxplot(div5$`diversity(wc.5.portreintamil)` ~ div5$habit)
boxplot(div8$`diversity(wc.8.portreintamil)` ~ div8$habit)
boxplot(div10$`diversity(wc.10.portreintamil)` ~ div10$habit)
dev.off()
###Rarefraction curve
spAbund <- rowSums(wc.10.portreintamil)  #gives the number of individuals found in each plot
raremin <- min(rowSums(wc.10.portreintamil))  #rarefaction uses the smallest number of observations per sample to extrapolate the expected number if all other samples only had that number of observations
sRare <- rarefy(wc.10.portreintamil, raremin) # now use function rarefy
rarecurve(wc.5.portreintamil, col = "blue")##### Distance measures #####
# In this script, plots and tables are made
# to analize different measures of word diversity
