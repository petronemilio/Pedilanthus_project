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


##Maybe include radial files
##### Load files of lempelziv measurments #####
lempel.by.cell <- read.csv("Data/lemplzivbyfile.csv", row.names=1)
#See a exploratory boxplot
boxplot(lempel.by.cell$Value ~ lempel.by.cell$Name)
###Add species identifiers to files
match.id <- match(lempel.by.cell$Name,species.id$files)
lempel.by.cell$species <- as.character(species.id$species[match.id])
match.id <- match(lempel.by.cell$Name,species.id$files)
lempel.by.cell$habit <- as.character(species.id$habit[match.id])
###boxplot with species
boxplot(lempel.by.cell$Value ~ lempel.by.cell$species)
#Remove Lsystems:
lempel.by.cell <- subset(lempel.by.cell, lempel.by.cell$Name != "contextmesicLsystem_edited_cells_NotConverge.txt" &
         lempel.by.cell$Name != "contextxericLsystem_edited_cells_NotConverge.txt" &
         lempel.by.cell$Name != "Ray_Lsystem_edited_cells_NotConverge.txt")
hist(lempel.by.cell$Value)
hist(log10(lempel.by.cell$Value))

lempel.aov <- aov(log10(lempel.by.cell$Value+1)~ lempel.by.cell$habit)
lempel.species.aov <- aov(log10(lempel.by.cell$Value+1)~ lempel.by.cell$species)
summary(lempel.aov)
plot(lempel.aov)
TukeyHSD(lempel.aov)
summary(lempel.species.aov)
plot(lempel.aov)
lempel.species.tukey <- TukeyHSD(lempel.species.aov)
plot(lempel.species.tukey, las = 1)
#Make palete
colores<- c("#71982d","#7464d7","#5aba46","#d369d3","#44be7a","#a2409a","#b8b436","#696fba",
            "#d68937","#5f9ed7","#c94c34","#3fc1bf","#da477c","#67b88c","#a14761","#91b869",
            "#c987c4","#4a772f","#dc8279","#327e58","#c2a864","#846a2b","#71982d")

pdf("Figures/lempelzivFile.pdf") # Para guardar en PDF
boxplot(lempel.by.cell$Value ~ lempel.by.cell$Name, col=colores,
        las=1.5, xlab="Cell-file Lempel-Ziv algorithm",ylab = "Score", cex=0.4,
        cex.names=0.5, cex.axis = 0.3, cex.lab=1.2)
dev.off()
pdf("Figures/lempelzivsp.pdf")
boxplot(lempel.by.cell$Value ~ lempel.by.cell$species, col=colores,
        las=1.5, xlab = "Cell-file Lempel-Ziv algorithm",ylab = "Score", cex=0.4,
        cex.names=0.5, cex.axis=0.6, cex.lab=1.2)
dev.off()
pdf("Figures/lempelzivhabit.pdf")
boxplot(lempel.by.cell$Value ~ lempel.by.cell$habit, 
        col=c("#FAD510","#CB2314","#CB2314","#273046"),
        xlab = "Cell-file Lempel-Ziv algorithm")
dev.off()
pdf("Figures/lempelzivhabit_log10.pdf")
boxplot((log10(lempel.by.cell$Value+1)) ~ lempel.by.cell$habit, 
        col=c("#FAD510","#CB2314","#CB2314","#273046"),
        xlab = "Cell-file Lempel-Ziv algorithm")
dev.off()
#
#lempel.cellhabit.lm <-lm(lempel.by.cell$Value~lempel.by.cell$habit)
#summary(lempel.cellhabit.lm)
#plot(lempel.cellhabit.lm)
library(emmeans)
library(multcompView)
emm1 <- emmeans(lempel.species.aov, specs = pairwise ~ species,adjust="tukey")
emm1$contrasts
multcomp::cld(emm1$emmeans, alpha = 0.10, Letters=LETTERS)
#
#lempel.cellsp.lm <-lm(lempel.by.cell$Value~lempel.by.cell$species)
#summary(lempel.cellsp.lm)
##
emm.habit.lempel <- emmeans(lempel.aov, specs = pairwise ~ habit, adjust="tukey")
emm.habit.lempel$contrasts
multcomp::cld(emm.habit.lempel$emmeans, alpha = 0.10, Letters=LETTERS)
#
pdf("Figures/lempelzivsp_reduced.pdf")
boxplot(lempel.by.cell$Value ~ lempel.by.cell$species,
        col=c("#273046", "#FAD510","#FAD510","#FAD510","#FAD510","#273046","#273046","#273046",
                       "#FAD510","#273046","#FAD510","#273046","#273046","#273046","#CB2314","#CB2314"),
        las=1.5, xlab = "Cell-file Lempel-Ziv algorithm",ylab = "Score", cex=0.4,
        cex.names=0.5, cex.axis=0.6, cex.lab=1.2)
dev.off()
pdf("Figures/lempelzivsp_reduced_log.pdf")
boxplot(log10(lempel.by.cell$Value) ~ lempel.by.cell$species,
        col=c("#273046", "#FAD510","#FAD510","#FAD510","#FAD510","#273046","#273046","#273046",
                       "#FAD510","#273046","#FAD510","#273046","#273046","#273046","#CB2314","#CB2314"),
                       las=1.5, xlab = "Cell-file Lempel-Ziv algorithm",ylab = "Score", cex=0.4,
        cex.names=0.5, cex.axis=0.6, cex.lab=1.2)
dev.off()

ggplot(lempel.by.cell, aes(x = Name, y = Value, color = species)) + geom_boxplot()
ggplot(lempel.by.cell, aes(x = species, y = Value, color = species)) + geom_boxplot()
#
species.id<-species.id[order(species.id$species),]
myColors <- ifelse(levels(as.factor(species.id$habit))=="mesic", "#71982d", 
             ifelse(levels(as.factor(species.id$habit))=="xeric","#c94c34",
                    "grey90"))
boxplot(lempel.by.cell$Value ~ lempel.by.cell$species , 
        col=myColors)
#Now plot distribution of compression
bract<-subset(lempel.by.cell, lempel.by.cell$species == "E. bracteata")
#Make loop to plot each file
plot(lempel.by.cell$Value[lempel.by.cell$species=="E. tithymaloides"], type ="l")
#####Make palete
pal1 <- wes_palette("BottleRocket2")
pal2 <- wes_palette("Rushmore1")
pal3 <- wes_palette("Darjeeling1")
pal4 <- wes_palette("FantasticFox1")
pal5 <- wes_palette("IsleofDogs1")
pal6 <- wes_palette("Cavalcanti1")
my_palette <- c(pal1,pal2,pal3,pal4,pal5,pal6)
#
pdf("Figures/lempelzivbyfile.pdf")
par(mfrow = c(2,1))
z=1
for(i in files) {
  temp <- subset(lempel.by.cell, lempel.by.cell$Name == i)
  plot(temp$Value, type ="l", xlab = "Cell file position", ylab="Lempel-ziv value",
       main= i, col= my_palette[z])
  z = z +1
}
dev.off()
##
plot(lempel.by.cell$Value[lempel.by.cell$species=="probL-system"], type ="b", xlab="Compression distance",
     col=my_palette[1], lwd=2)
plot(lempel.by.cell$Value[lempel.by.cell$species=="probetaL-system"], type ="b", xlab="Compression distance",
     col=my_palette[1], lwd=2)

z=1
for(i in files) {
  plot(lempel.by.cell$Value[lempel.by.cell$Name==i], type ="b", xlab="Compression distance",
       col=my_palette[1], lwd=2,  main= i)
  lines(lempel.by.cell$Value[lempel.by.cell$Name=="probLsystem_edited_cells.txt"], type ="b", xlab="Compression distance",
        col = my_palette[z], lwd = 2)
  z=z+1
}
#
pdf("Figures/lempelziv_bycell.pdf")
plot(lempel.by.cell$Value[lempel.by.cell$Name=="974_edited_cells.txt"], type ="b", xlab="Compression distance",
     col=my_palette[1], lwd=2,  main= "Lempel-ziv")
lines(lempel.by.cell$Value[lempel.by.cell$Name=="probLsystem_edited_cells.txt"], type ="b", xlab="Compression distance",
      col = my_palette[2], lwd = 2)
lines(lempel.by.cell$Value[lempel.by.cell$Name=="EPM10_edited_cells.txt"], type ="b", xlab="Compression distance",
      col = my_palette[3], lwd = 2)
legend("bottom", legend = c("E. peritropoides","L-system","E. diazlunana"),
       col =levels(as.factor(my_palette[1:3])) , pch = 16, inset =-0.15,xpd=TRUE,horiz=TRUE)
dev.off()
#
par(mfrow = c(2,1))
plot(lempel.by.cell$Value[lempel.by.cell$species=="probL-system"], type ="l", xlab="Compression distance")
plot(lempel.by.cell$Value[lempel.by.cell$Name=="EPM5_edited_cells.txt"], type ="l", xlab="Compression distance")
dev.off()
#
lempel.by.cell <- transform(lempel.by.cell, Sequence=ave(seq_along(Name), Name, FUN=seq_along))
ggplot(lempel.by.cell, aes(x = Sequence,y = Value,color = species)) +
  geom_line()
#Add identifier of habit
ggplot(lempel.by.cell, aes(x = Sequence,y = Value,color = habit)) +
  geom_line()

######## Shannon entropy #######
shannonentropy.by.cell <- read.csv("Data/shannonentropy.csv", row.names=1)
#Add species and habit
match.id <- match(shannonentropy.by.cell$Name,species.id$files)
shannonentropy.by.cell$species <- as.character(species.id$species[match.id])
shannonentropy.by.cell$habit <- as.character(species.id$habit[match.id])
#
pdf("Figures/shannonbyfile.pdf") # Para guardar en PDF
boxplot(shannonentropy.by.cell$Value ~ shannonentropy.by.cell$Name, col=colores,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        cex.names=0.5, cex.axis=0.3, cex.lab=1.2)
dev.off()
pdf("Figures/shanonbysp.pdf")
boxplot(shannonentropy.by.cell$Value ~ shannonentropy.by.cell$species, col=colores,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        names=c( expression(italic("E. bracteata")),expression(italic("E. calcarata")) ,
                 expression(italic("E. coalcomanensis")),expression(italic("E. colligata")),
                 expression(italic("E. conzattii")),expression(italic("E. cymbifera")),
                 expression(italic("E. cyri")), expression(italic("E. diazlunana")),
                 expression(italic("E. finkii")),expression(italic("E. lomelli")),
                 expression(italic("E. peritropoides")),  expression(italic("E. personata")),
                 expression(italic("E. tehuacana")),expression(italic("E. tithymaloides")),
                 expression("L-systemM"),expression("L-systemX"), expression("probBetaLsystem"),
                 expression("probL-system"),expression("RayL-system")),
        cex.names=0.5, cex.axis=0.8, cex.lab=1.2)
dev.off()
#
pdf("Figures/entropybyfile.pdf")
par(mfrow = c(2,1))
z=1
for(i in files) {
  temp <- subset(shannonentropy.by.cell, shannonentropy.by.cell$Name == i)
  plot(temp$Value, type ="l", xlab = "Cell file position", ylab="Shannon-Entropy value",
       main= i, col= my_palette[z], lwd = 2)
  z = z +1
}
dev.off()
######
#Small subset
shannon_subset<-subset(shannonentropy.by.cell, Name == "probLsystem_edited_cells.txt" | Name=="974_edited_cells.txt" |
                         Name == "EPM10_edited_cells.txt")
shannon_subset <- transform(shannon_subset, Sequence=ave(seq_along(species),
                                                         species, FUN=seq_along))
#Shannon
pdf("Figures/shannon_bycell.pdf")
ggplot(shannon_subset, aes(x = Sequence,y = Value,color = species)) +
  geom_line()
dev.off()
#########Remove rlsystem, and context sensitive lsystems
shannonentropy.by.cell_subset <-subset(shannonentropy.by.cell, shannonentropy.by.cell$Name != "contextmesicLsystem_edited_cells_NotConverge.txt" &
                                         shannonentropy.by.cell$Name != "contextxericLsystem_edited_cells_NotConverge.txt" &
                                         shannonentropy.by.cell$Name != "Ray_Lsystem_edited_cells_NotConverge.txt" )
#
hist(shannonentropy.by.cell_subset$Value)
shannon_morethanzerolessonefive <- subset(shannonentropy.by.cell_subset, shannonentropy.by.cell_subset$Value > 0.1 & 
                                 shannonentropy.by.cell_subset$Value < 1.5)
hist(shannon_morethanzerolessonefive$Value)
shannon.cellsp.aov <- aov(shannonentropy.by.cell_subset$Value ~
                         factor(shannonentropy.by.cell_subset$species))
plot(shannon.cellsp.aov)
shannon.cellsp.filtered.aov <- aov(shannon_morethanzerolessonefive$Value ~
                          shannon_morethanzerolessonefive$species)
plot(shannon.cellsp.filtered.aov)
shannon.cellsp.filtered.aov
shannon.cellsp.filtered.kruskal <- kruskal.test(shannonentropy.by.cell_subset$Value ~
  shannonentropy.by.cell_subset$species)
shannon.cellhabit.filtered.kruskal <- kruskal.test(shannonentropy.by.cell_subset$Value ~
                                                  shannonentropy.by.cell_subset$habit)
#Check if dunn test should follow.
#########################################
emm.shannon <- emmeans(shannon.cellsp.aov, specs = pairwise ~ species,adjust="tukey")
emm.shannon$contrasts
multcomp::cld(emm.shannon$emmeans, alpha = 0.10, Letters=LETTERS)
emm.shannon.filt <- emmeans(shannon.cellsp.filtered.aov, specs = pairwise ~ species,adjust="tukey")
emm.shannon.filt$contrasts
multcomp::cld(emm.shannon.filt$emmeans, alpha = 0.10, Letters=LETTERS)

#
pdf("Figures/shanonbysp_subset.pdf")
boxplot(shannonentropy.by.cell_subset$Value ~ shannonentropy.by.cell_subset$species,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        names=c( expression(italic("E. bracteata")),expression(italic("E. calcarata")) ,
                 expression(italic("E. coalcomanensis")),expression(italic("E. colligata")),
                 expression(italic("E. conzattii")),expression(italic("E. cymbifera")),
                 expression(italic("E. cyri")), expression(italic("E. diazlunana")),
                 expression(italic("E. finkii")),expression(italic("E. lomelli")),
                 expression(italic("E. peritropoides")),  expression(italic("E. personata")),
                 expression(italic("E. tehuacana")),expression(italic("E. tithymaloides")),
                 expression("probBetaLsystem"),expression("probL-system")),
        cex.names=0.5, cex.axis=0.8, cex.lab=1.2)
dev.off()

#                
pdf("Figures/shanonbysp_subset_filtered.pdf")
boxplot(shannon_morethanzerolessonefive$Value ~ shannon_morethanzerolessonefive$species,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        col=c("#273046", "#FAD510","#FAD510","#FAD510","#FAD510","#273046","#273046","#273046",
              "#FAD510","#273046","#FAD510","#273046","#273046","#273046","#CB2314","#CB2314"),
        names=c( expression(italic("E. bracteata")),expression(italic("E. calcarata")) ,
                 expression(italic("E. coalcomanensis")),expression(italic("E. colligata")),
                 expression(italic("E. conzattii")),expression(italic("E. cymbifera")),
                 expression(italic("E. cyri")), expression(italic("E. diazlunana")),
                 expression(italic("E. finkii")),expression(italic("E. lomelli")),
                 expression(italic("E. peritropoides")),  expression(italic("E. personata")),
                 expression(italic("E. tehuacana")),expression(italic("E. tithymaloides")),
                 expression("probBetaLsystem"),expression("probL-system")),
        cex.names=0.5, cex.axis=0.8, cex.lab=1.2)
dev.off()
#by category
pdf("Figures/shanonbyhabit.pdf") # Para guardar en PDF
boxplot(shannon_morethanzerolessonefive$Value ~ shannon_morethanzerolessonefive$habit,
        col=c("#FAD510","#CB2314","#CB2314","#273046"),
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4)
dev.off()    
boxplot(shannonentropy.by.cell_subset$Value ~ shannonentropy.by.cell_subset$habit,
        col=c("#FAD510","#CB2314","#CB2314","#273046"),
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4)

#######Checkout files with ray cells########
shannonentropy.by.cell.with.ray <- read.csv("Data/shannonentropywithR.csv", row.names=1)

#Add species and habit
match.id <- match(shannonentropy.by.cell.with.ray$Name,species.id$files)
shannonentropy.by.cell.with.ray$species <- as.character(species.id$species[match.id])
shannonentropy.by.cell.with.ray$habit <- as.character(species.id$habit[match.id])

boxplot(shannonentropy.by.cell.with.ray$Value ~ shannonentropy.by.cell.with.ray$Name,
        col=colores,las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        cex.names=0.5, cex.axis=0.3, cex.lab=1.2)

plot(shannonentropy.by.cell.with.ray$Value[shannonentropy.by.cell.with.ray$Name=="974_edited_cells.txt"],
     type ="l", xlab="Compression distance")
lines(shannonentropy.by.cell$Value[shannonentropy.by.cell$Name=="974_edited_cells.txt"],
      type ="l", col =my_palette[2])
######Check out Shanon by window #####
shannonentropy.by.window <- read.csv("Data/shannonentropy_window.csv", row.names=1)
library(dplyr)
shannonentropy.by.window <- shannonentropy.by.window %>%
  group_by(Sample,Lineage) %>% mutate(id = row_number())

write.csv(shannonentropy.by.window, "Data/shannonentropy_window_id.csv")
window947 <- subset(shannonentropy.by.window, Sample == "974_edited_cells.txt")
window883 <- subset(shannonentropy.by.window, Sample == "883_edited_cells.txt")
window883 <- subset(shannonentropy.by.window, Sample == "883_edited_cells.txt")
windowepm10 <- subset(shannonentropy.by.window, Sample == "EPM10_edited_cells.txt")
prbeta <-subset(shannonentropy.by.window, Sample == "probLsystembeta_edited_cells.txt")
#PLot every sample
samples<-factor(shannonentropy.by.window$Sample)
for (i in samples){
  plot(shannonentropy.by.window$id[shannonentropy.by.window==i],
       shannonentropy.by.window$Value[shannonentropy.by.window==i],type = "l",lwd = 3)
}
plot(window947$id,window947$Value,type = "l")
plot(window947$id[window947$Lineage==26],window947$Value[window947$Lineage==26],type = "l",lwd = 2,
     xlab="Window position",ylab = "Shannon Value")

lineage<-factor(windowepm10$Lineage)
for (i in lineage){
  points(windowepm10$id[windowepm10$Lineage==i], windowepm10$Value[windowepm10$Lineage==i],
         col="green",type = "l",lwd = 3,
         col=pal1)
}
plot(windowepm10$id,windowepm10$Value, type = "l",col=pal1)
plot(window883$id,window883$Value,type = "l",lwd = 3)


pdf("Figures/shannon_window947.pdf")
plot(window947$id[window947$Lineage==2],window947$Value[window947$Lineage==2],type = "l",lwd = 3)
points(window947$id[window947$Lineage==3],window947$Value[window947$Lineage==3], col="green",type = "l",lwd = 3)
points(window947$id[window947$Lineage==4],window947$Value[window947$Lineage==4], col="red",type = "l",lwd = 3)
points(window947$id[window947$Lineage==5],window947$Value[window947$Lineage==5], col=pal1[1],type = "l",lwd = 3)
points(window947$id[window947$Lineage==6],window947$Value[window947$Lineage==6], col=pal1[2],type = "l",lwd = 3)
#points(window947$id[window947$Lineage==7],window947$Value[window947$Lineage==7], col=pal1[3],pch=19)
#points(window947$id[window947$Lineage==8],window947$Value[window947$Lineage==8], col=pal2[1],pch=19)
#points(window947$id[window947$Lineage==9],window947$Value[window947$Lineage==9], col=pal2[2],pch=19)
#points(window947$id[window947$Lineage==10],window947$Value[window947$Lineage==10], col=pal2[3],pch=19)
dev.off()

pdf("Figures/shannon_window883.pdf")
plot(window883$id[window883$Lineage==3],window883$Value[window883$Lineage==3],type = "l",lwd = 3)
points(window883$id[window883$Lineage==2],window883$Value[window883$Lineage==2], col="green",type = "l",lwd = 3)
points(window883$id[window883$Lineage==6],window883$Value[window883$Lineage==6], col="red",type = "l",lwd = 3)
points(window883$id[window883$Lineage==4],window883$Value[window883$Lineage==4], col=pal1[1],type = "l",lwd = 3)
points(window883$id[window883$Lineage==7],window883$Value[window883$Lineage==7], col=pal1[2],type = "l",lwd = 3)
dev.off()

pdf("Figures/probetaS.pdf")
plot(prbeta$id[prbeta$Lineage==6],prbeta$Value[prbeta$Lineage==6],type = "l",lwd = 3)
points(prbeta$id[prbeta$Lineage==2],prbeta$Value[prbeta$Lineage==2], col="green",type = "l",lwd = 3)
points(prbeta$id[prbeta$Lineage==1],prbeta$Value[prbeta$Lineage==1], col="red",type = "l",lwd = 3)
points(prbeta$id[prbeta$Lineage==4],prbeta$Value[prbeta$Lineage==4], col=pal1[1],type = "l",lwd = 3)
points(prbeta$id[prbeta$Lineage==7],prbeta$Value[prbeta$Lineage==7], col=pal1[2],type = "l",lwd = 3)
dev.off()

summary(cell_lengths)

ggplot(window947, aes(x= id, y = Value, color=Lineage)) +
  geom_point()

ggplot(lempel.by.cell, aes(x = Name, y = Value, color = species)) + geom_boxplot()


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
## Loading euclidean distance files
euc_dist2 <- read.csv("Data/euclidean_distance_morethanone/euclidean_distance2.csv")
euc_dist3 <- read.csv("Data/euclidean_distance_morethanone/euclidean_distance3.csv")
euc_dist4 <- read.csv("Data/euclidean_distance_morethanone/euclidean_distance4.csv")
euc_dist5 <- read.csv("Data/euclidean_distance_morethanone/euclidean_distance5.csv")
euc_dist10 <- read.csv("Data/euclidean_distance_morethanone/euclidean_distance10.csv")
euc_dist15 <- read.csv("Data/euclidean_distance_morethanone/euclidean_distance15.csv")
euc_dist19 <- read.csv("Data/euclidean_distance_morethanone/euclidean_distance19.csv")
euc_dist22 <- read.csv("Data/euclidean_distance_morethanone/euclidean_distance22.csv")
euc_dist25 <- read.csv("Data/euclidean_distance_morethanone/euclidean_distance25.csv")

##Euc_dist objects are not properly ordered to make matrix
euc_dist2[c("File1", "File2")] <- lapply(euc_dist2[c("File1", "File2")], factor,
                                         levels=unique(unlist(euc_dist2[c("File1", "File2")])))
#
euc_dist3[c("File1", "File2")] <- lapply(euc_dist3[c("File1", "File2")], factor,
                                         levels=unique(unlist(euc_dist3[c("File1", "File2")])))
#
euc_dist4[c("File1", "File2")] <- lapply(euc_dist4[c("File1", "File2")], factor,
                                         levels=unique(unlist(euc_dist4[c("File1", "File2")])))
#
euc_dist5[c("File1", "File2")] <- lapply(euc_dist5[c("File1", "File2")], factor,
                                         levels=unique(unlist(euc_dist5[c("File1", "File2")])))
#
euc_dist10[c("File1", "File2")] <- lapply(euc_dist10[c("File1", "File2")], factor,
                                         levels=unique(unlist(euc_dist10[c("File1", "File2")])))
#
euc_dist15[c("File1", "File2")] <- lapply(euc_dist15[c("File1", "File2")], factor,
                                         levels=unique(unlist(euc_dist15[c("File1", "File2")])))
#
euc_dist19[c("File1", "File2")] <- lapply(euc_dist19[c("File1", "File2")], factor,
                                         levels=unique(unlist(euc_dist19[c("File1", "File2")])))
#
euc_dist22[c("File1", "File2")] <- lapply(euc_dist22[c("File1", "File2")], factor,
                                         levels=unique(unlist(euc_dist22[c("File1", "File2")])))
#
euc_dist25[c("File1", "File2")] <- lapply(euc_dist25[c("File1", "File2")], factor,
                                          levels=unique(unlist(euc_dist25[c("File1", "File2")])))
#####Transform to matrix format#####
matrix_euc2 <- xtabs(Euc_dist ~ File1+ File2,euc_dist2)
corrplot(matrix_euc2, is.corr = FALSE)
##
matrix_euc3 <- xtabs(Euc_dist ~ File1+ File2,euc_dist3)
corrplot(matrix_euc3, is.corr = FALSE,type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05,
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=0.5, outline=T)
##
matrix_euc4 <- xtabs(Euc_dist ~ File1+ File2,euc_dist4)
corrplot(matrix_euc4, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05,
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
##
matrix_euc5 <- xtabs(Euc_dist ~ File1+ File2,euc_dist5)
corrplot(matrix_euc5, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05,
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
###
matrix_euc6 <- xtabs(Euc_dist ~ File1+ File2,euc_dist6)
corrplot(matrix_euc6, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05,
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
#
matrix_euc7 <- xtabs(Euc_dist ~ File1+ File2,euc_dist7)
corrplot(matrix_euc7, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05,
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
#
matrix_euc8 <- xtabs(Euc_dist ~ File1+ File2,euc_dist8)
corrplot(matrix_euc8, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05,
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
#
matrix_euc9 <- xtabs(Euc_dist ~ File1+ File2,euc_dist9)
corrplot(matrix_euc9, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05,
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)
#
matrix_euc10 <- xtabs(Euc_dist ~ File1+ File2,euc_dist10)
corrplot(matrix_euc10, is.corr = FALSE, type="upper", #order="hclust",   #p.mat = leng.leaf.matrix$P, sig.level = 0.05,
         bg="WHITE", tl.col = "black", tl.srt = 45, pch.cex=1, outline=T)

#####Load word counts to calc euc distance in R #####
wc.2 <- as.matrix(read.csv("Data/wordcounts2.csv",row.names=1))
wc.3 <- as.matrix(read.csv("Data/wordcounts3.csv", row.names = 1))
wc.3[is.na(wc.3)] <- 0
wc.4 <- as.matrix(read.csv("Data/wordcounts4.csv", row.names = 1))
wc.4[is.na(wc.4)] <- 0
wc.5 <- as.matrix(read.csv("Data/wordcounts5.csv", row.names = 1))
wc.5[is.na(wc.5)] <- 0
wc.6 <- as.matrix(read.csv("Data/wordcounts6.csv", row.names = 1))
wc.6[is.na(wc.6)] <- 0
wc.7 <- as.matrix(read.csv("Data/wordcounts7.csv", row.names = 1))
wc.7[is.na(wc.7)] <- 0
wc.8 <- as.matrix(read.csv("Data/wordcounts8.csv", row.names = 1))
wc.8[is.na(wc.8)] <- 0
wc.9 <- as.matrix(read.csv("Data/wordcounts9.csv", row.names = 1))
wc.9[is.na(wc.9)] <- 0
wc.10 <- as.matrix(read.csv("Data/wordcounts10.csv", row.names = 1))
wc.10[is.na(wc.10)] <- 0
wc.11 <- as.matrix(read.csv("Data/wordcounts11.csv", row.names = 1))
wc.11[is.na(wc.11)] <- 0
wc.12 <- as.matrix(read.csv("Data/wordcounts12.csv", row.names = 1))
wc.12[is.na(wc.12)] <- 0
wc.13 <- as.matrix(read.csv("Data/wordcounts13.csv",row.names = 1))
wc.13[is.na(wc.13)] <- 0
wc.14 <- as.matrix(read.csv("Data/wordcounts14.csv",row.names = 1))
wc.14[is.na(wc.14)] <- 0
#
euc.dist.r.2 <- dist(wc.2)
euc.dist.r.3 <- dist(wc.3)
euc.dist.r.4 <- dist(wc.4)
euc.dist.r.5 <- dist(wc.5)
euc.dist.r.6 <- dist(wc.6)
euc.dist.r.7 <- dist(wc.7)
euc.dist.r.8 <- dist(wc.8)
euc.dist.r.9 <- dist(wc.9)
euc.dist.r.10 <- dist(wc.10)
#Get proportion of counts
wc.2.freq <- wc.2/rowSums(wc.2)
euc.dist.r2.freq <- dist(wc.2.freq)
plot(hclust(euc.dist.r2.freq))
#
wc.3.freq <- wc.3/rowSums(wc.3)
euc.dist.r3.freq <- dist(wc.3.freq)
plot(hclust(euc.dist.r3.freq))
#
wc.4.freq <- wc.4/rowSums(wc.4)
euc.dist.r4.freq <- dist(wc.4.freq)
plot(hclust(euc.dist.r4.freq))
#
wc.5.freq <- wc.5/rowSums(wc.5)
euc.dist.r5.freq <- dist(wc.5.freq)
plot(hclust(euc.dist.r5.freq))
#
wc.6.freq <- wc.6/rowSums(wc.6)
euc.dist.r6.freq <- dist(wc.6.freq)
plot(hclust(euc.dist.r6.freq))
#
wc.7.freq <- wc.7/rowSums(wc.7)
euc.dist.r7.freq <- dist(wc.7.freq)
plot(hclust(euc.dist.r7.freq))
#
wc.8.freq <- wc.8/rowSums(wc.8)
euc.dist.r8.freq <- dist(wc.8.freq)
plot(hclust(euc.dist.r8.freq))
#
wc.10.freq <- wc.10/rowSums(wc.10)
euc.dist.r10.freq <- dist(wc.10.freq)
plot(hclust(euc.dist.r10.freq))
#Euclidean distances calculated with r are the same as the ones in python
#Try to solve problem of different total of cells
boxplot(rowSums(wc.2))
hist(rowSums(wc.3))
hist(rowSums(wc.3))
hist(rowSums(wc.4))
rowSums(wc.10)
#Try to count cell number for each and calc the counts per 30,000
rowSums(wc.2)
dosporcadatreintamil <- rowSums(wc.2)/30000
wc.2.portreintamil <- round((wc.2/dosporcadatreintamil),0)
wc.2.portreintamil.freq <- wc.2.portreintamil/rowSums(wc.2.portreintamil)
euc.dist.r.2.freq.std <- dist(wc.2.portreintamil.freq)
plot(hclust(euc.dist.r.2.freq.std))
#para 3
tresporcadatreintamil <- rowSums(wc.3)/30000
wc.3.portreintamil <- round((wc.3/tresporcadatreintamil),0)
wc.3.portreintamil.freq <- wc.3.portreintamil/rowSums(wc.3.portreintamil)
euc.dist.r.3.freq.std <- dist(wc.3.portreintamil.freq)
plot(hclust(euc.dist.r.3.freq.std))
#wc.3.portreintamil["probLsystem_edited_cells.txt",]
#para 4
cuatroporcadatreintamil <- rowSums(wc.4)/30000
wc.4.portreintamil <- round((wc.4/cuatroporcadatreintamil),0)
wc.4.portreintamil.freq <- wc.4.portreintamil/rowSums(wc.4.portreintamil)
euc.dist.r.4.std.pr <- dist(wc.4.portreintamil)
euc.dist.r.4.freq.std <- dist(wc.4.portreintamil.freq)
plot(hclust(euc.dist.r.4.freq.std))
plot(hclust(euc.dist.r.4.std.pr))

#para 5
cincoporcadatreintamil <- rowSums(wc.5)/30000
wc.5.portreintamil <- round((wc.5/cincoporcadatreintamil),0)
wc.5.portreintamil.freq <- wc.5.portreintamil/rowSums(wc.5.portreintamil)
euc.dist.r.5.freq.std <- dist(wc.5.portreintamil.freq)
plot(hclust(euc.dist.r.5.freq.std))
#para 6
seisporcadatreintamil <- rowSums(wc.6)/30000
wc.6.portreintamil <- round((wc.6/seisporcadatreintamil),0)
wc.6.portreintamil.freq <- wc.6.portreintamil/rowSums(wc.6.portreintamil)
euc.dist.r.6.freq.std <- dist(wc.6.portreintamil.freq)
plot(hclust(euc.dist.r.6.freq.std))
#Para 7
sieteporcadatreintamil <- rowSums(wc.7)/30000
wc.7.portreintamil <- round((wc.7/sieteporcadatreintamil),0)
wc.7.portreintamil.freq <- wc.7.portreintamil/rowSums(wc.7.portreintamil)
euc.dist.r.7.freq.std <- dist(wc.7.portreintamil.freq)
plot(hclust(euc.dist.r.7.freq.std))
#
#para 8
ochoporcadatreintamil <- rowSums(wc.8)/30000
wc.8.portreintamil <- round((wc.8/ochoporcadatreintamil),0)
#Para 9
nueveporcadatreintamil <- rowSums(wc.9)/30000
wc.9.portreintamil <- round((wc.9/nueveporcadatreintamil),0)
#para 10
diezporcadatreintamil <- rowSums(wc.10)/30000
wc.10.portreintamil <- round((wc.10/diezporcadatreintamil),0)
wc.10.portreintamil.freq <- wc.10.portreintamil/rowSums(wc.10.portreintamil)
euc.dist.r.10.freq.std <- dist(wc.10.portreintamil.freq)
plot(hclust(euc.dist.r.10.freq.std))

#####
euc.dist.r.2.std <- dist(wc.2.portreintamil)
euc.dist.r.3.std <- dist(wc.3.portreintamil)
euc.dist.r.4.std <- dist(wc.4.portreintamil)
euc.dist.r.5.std <- dist(wc.5.portreintamil)
euc.dist.r.6.std <- dist(wc.6.portreintamil)
euc.dist.r.7.std <- dist(wc.7.portreintamil)
euc.dist.r.8.std <- dist(wc.8.portreintamil)
euc.dist.r.9.std <- dist(wc.9.portreintamil)
euc.dist.r.10.std <- dist(wc.10.portreintamil)
####
plot(hclust(euc.dist.r.5.std))

pdf("Figures/eucdistwordlength2.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r2.freq))
plot(hclust(euc.dist.r.2.freq.std))
dev.off()
pdf("Figures/eucdistwordlength2abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.2))
plot(hclust(euc.dist.r.2.std))
dev.off()
pdf("Figures/eucdistwordlength3.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r3.freq))
plot(hclust(euc.dist.r.3.freq.std))
dev.off()
pdf("Figures/eucdistwordlength3abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.3))
plot(hclust(euc.dist.r.3.std))
dev.off()
pdf("Figures/eucdistwordlength4.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r4.freq))
plot(hclust(euc.dist.r.4.freq.std))
dev.off()
pdf("Figures/eucdistwordlength4abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.4))
plot(hclust(euc.dist.r.4.std))
dev.off()
pdf("Figures/eucdistwordlength5.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r5.freq))
plot(hclust(euc.dist.r.5.freq.std))
dev.off()
pdf("Figures/eucdistwordlength5abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.5))
plot(hclust(euc.dist.r.5.std))
dev.off()
pdf("Figures/eucdistwordlength6.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r6.freq))
plot(hclust(euc.dist.r.6.freq.std))
dev.off()
pdf("Figures/eucdistwordlength6abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.6))
plot(hclust(euc.dist.r.6.std))
dev.off()
pdf("Figures/eucdistwordlength7.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r7.freq))
plot(hclust(euc.dist.r.7.freq.std))
dev.off()
pdf("Figures/eucdistwordlength7abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.7))
plot(hclust(euc.dist.r.7.std))
dev.off()
pdf("Figures/eucdistwordlength10.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r10.freq))
plot(hclust(euc.dist.r.10.freq.std))
dev.off()
pdf("Figures/eucdistwordlength10abs.pdf") # Para guardar en PDF
par(mfrow= c(1,2))
plot(hclust(euc.dist.r.10))
plot(hclust(euc.dist.r.10.std))
dev.off()
####
#heatmap.2(as.matrix(euc.dist.r.2),key=TRUE,scale = "row",
#            margins = c(10,12), cexRow=0.5)

pheatmap(as.matrix(euc.dist.r.2.std), display_numbers = F)

pdf("Figures/heatmaplength2.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r2.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.2.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.2), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.2.std), display_numbers = F)
dev.off()
pdf("Figures/heatmaplength3.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r3.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.3.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.3), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.3.std), display_numbers = F)
dev.off()

pdf("Figures/heatmaplength4.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r4.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.4.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.4), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.4.std), display_numbers = F)
dev.off()
pdf("Figures/heatmaplength5.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r5.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.5.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.5), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.5.std), display_numbers = F)
dev.off()
pdf("Figures/heatmaplength6.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r6.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.6.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.6), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.6.std), display_numbers = F)
dev.off()
pdf("Figures/heatmaplength7.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r7.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.7.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.7), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.7.std), display_numbers = F)
dev.off()
pdf("Figures/heatmaplength10.pdf") # Para guardar en PDF
pheatmap(as.matrix(euc.dist.r10.freq), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.10.freq.std), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.10), display_numbers = F)
pheatmap(as.matrix(euc.dist.r.10.std), display_numbers = F)
dev.off()
##Falta hacer los análisis de words incluyendo el sistema radial

##### Load files of lempelziv measurments #####
lempel.by.cell <- read.csv("Data/lemplzivbyfile.csv", row.names=1)
#See a exploratory boxplot
boxplot(lempel.by.cell$Value ~ lempel.by.cell$Name)
###Add species identifiers to files
match.id <- match(lempel.by.cell$Name,species.id$files)
lempel.by.cell$species <- as.character(species.id$species[match.id])
match.id <- match(lempel.by.cell$Name,species.id$files)
lempel.by.cell$habit <- as.character(species.id$habit[match.id])
###boxplot with species
boxplot(lempel.by.cell$Value ~ lempel.by.cell$species)
#Make palete
colores<- c("#71982d","#7464d7","#5aba46","#d369d3","#44be7a","#a2409a","#b8b436","#696fba",
            "#d68937","#5f9ed7","#c94c34","#3fc1bf","#da477c","#67b88c","#a14761","#91b869",
            "#c987c4","#4a772f","#dc8279","#327e58","#c2a864","#846a2b","#71982d")

pdf("Figures/lempelzivFile.pdf") # Para guardar en PDF
boxplot(lempel.by.cell$Value ~ lempel.by.cell$Name, col=colores,
        las=1.5, xlab="Cell-file Lempel-Ziv algorithm",ylab = "Score", cex=0.4,
        cex.names=0.5, cex.axis = 0.3, cex.lab=1.2)
dev.off()
pdf("Figures/lempelzivsp.pdf")
boxplot(lempel.by.cell$Value ~ lempel.by.cell$species, col=colores,
        las=1.5, xlab = "Cell-file Lempel-Ziv algorithm",ylab = "Score", cex=0.4,
        cex.names=0.5, cex.axis=0.6, cex.lab=1.2)
dev.off()
pdf("Figures/lempelzivhabit.pdf")
boxplot(lempel.by.cell$Value ~ lempel.by.cell$habit, col=colores,
        xlab = "Cell-file Lempel-Ziv algorithm")
dev.off()

lempel.cellhabit.lm <-lm(lempel.by.cell$Value~lempel.by.cell$habit)
summary(lempel.cellhabit.lm)
plot(lempel.cellhabit.lm)
library(emmeans)
library(multcompView)
emm1 <- emmeans(lempel.cellhabit.lm, specs = pairwise ~ habit,adjust="tukey")
emm1$contrasts
multcomp::cld(emm1$emmeans, alpha = 0.10, Letters=LETTERS)

lempel.cellsp.lm <-lm(lempel.by.cell$Value~lempel.by.cell$species)
summary(lempel.cellsp.lm)
##Remove rlsystem, and context sensitive lsystems
lempel.by.cell_subset <-subset(lempel.by.cell, lempel.by.cell$species != "RayL-system")
lempel.by.cell_subset <-subset(lempel.by.cell_subset, lempel.by.cell_subset$species != "L-systemMesic")
lempel.by.cell_subset <-subset(lempel.by.cell_subset, lempel.by.cell_subset$species != "L-systemXeric")
lempel.by.cell_subset$species <- factor(lempel.by.cell_subset$species)
#
lempel.cellsp.lm <-lm(lempel.by.cell_subset$Value ~ factor(lempel.by.cell_subset$species))
emm.lempel <- emmeans(lempel.cellsp.lm, specs = pairwise ~ species,adjust="tukey")
emm.lempel$contrasts
multcomp::cld(emm.lempel$emmeans, alpha = 0.10, Letters=LETTERS)
#
pdf("Figures/lempelzivsp_reduced.pdf")
boxplot(lempel.by.cell_subset$Value ~ lempel.by.cell_subset$species,
        las=1.5, xlab = "Cell-file Lempel-Ziv algorithm",ylab = "Score", cex=0.4,
        cex.names=0.5, cex.axis=0.6, cex.lab=1.2)
dev.off()
ggplot(lempel.by.cell, aes(x = Name, y = Value, color = species)) + geom_boxplot()
ggplot(lempel.by.cell, aes(x = species, y = Value, color = species)) + geom_boxplot()

#Now plot distribution of compression
bract<-subset(lempel.by.cell, lempel.by.cell$species == "E. bracteata")
#Make loop to plot each file
plot(lempel.by.cell$Value[lempel.by.cell$species=="E. tithymaloides"], type ="l")

#####Make palete
pal1 <- wes_palette("BottleRocket2")
pal2 <- wes_palette("Rushmore1")
pal3 <- wes_palette("Darjeeling1")
pal4 <- wes_palette("FantasticFox1")
pal5 <- wes_palette("IsleofDogs1")
pal6 <- wes_palette("Cavalcanti1")
my_palette <- c(pal1,pal2,pal3,pal4,pal5,pal6)
#
pdf("Figures/lempelzivbyfile.pdf")
par(mfrow = c(2,1))
z=1
for(i in files) {
  temp <- subset(lempel.by.cell, lempel.by.cell$Name == i)
  plot(temp$Value, type ="l", xlab = "Cell file position", ylab="Lempel-ziv value",
       main= i, col= my_palette[z])
  z = z +1
}
dev.off()
##
plot(lempel.by.cell$Value[lempel.by.cell$species=="probL-system"], type ="b", xlab="Compression distance",
     col=my_palette[1], lwd=2)
plot(lempel.by.cell$Value[lempel.by.cell$species=="probetaL-system"], type ="b", xlab="Compression distance",
     col=my_palette[1], lwd=2)

z=1
for(i in files) {
  plot(lempel.by.cell$Value[lempel.by.cell$Name==i], type ="b", xlab="Compression distance",
       col=my_palette[1], lwd=2,  main= i)
  lines(lempel.by.cell$Value[lempel.by.cell$Name=="probLsystem_edited_cells.txt"], type ="b", xlab="Compression distance",
        col = my_palette[z], lwd = 2)
  z=z+1
}
#
pdf("Figures/lempelziv_bycell.pdf")
plot(lempel.by.cell$Value[lempel.by.cell$Name=="974_edited_cells.txt"], type ="b", xlab="Compression distance",
     col=my_palette[1], lwd=2,  main= "Lempel-ziv")
lines(lempel.by.cell$Value[lempel.by.cell$Name=="probLsystem_edited_cells.txt"], type ="b", xlab="Compression distance",
      col = my_palette[2], lwd = 2)
lines(lempel.by.cell$Value[lempel.by.cell$Name=="EPM10_edited_cells.txt"], type ="b", xlab="Compression distance",
      col = my_palette[3], lwd = 2)
legend("bottom", legend = c("E. peritropoides","L-system","E. diazlunana"),
       col =levels(as.factor(my_palette[1:3])) , pch = 16, inset =-0.15,xpd=TRUE,horiz=TRUE)
dev.off()
#


par(mfrow = c(2,1))
plot(lempel.by.cell$Value[lempel.by.cell$species=="probL-system"], type ="l", xlab="Compression distance")
plot(lempel.by.cell$Value[lempel.by.cell$Name=="EPM5_edited_cells.txt"], type ="l", xlab="Compression distance")

dev.off()
#
lempel.by.cell <- transform(lempel.by.cell, Sequence=ave(seq_along(Name), Name, FUN=seq_along))

ggplot(lempel.by.cell, aes(x = Sequence,y = Value,color = species)) +
  geom_line()
#Add identifier of habit
ggplot(lempel.by.cell, aes(x = Sequence,y = Value,color = habit)) +
  geom_line()

######## Shannon entropy #######
shannonentropy.by.cell <- read.csv("Data/shannonentropy.csv", row.names=1)

#Add species and habit
match.id <- match(shannonentropy.by.cell$Name,species.id$files)
shannonentropy.by.cell$species <- as.character(species.id$species[match.id])
shannonentropy.by.cell$habit <- as.character(species.id$habit[match.id])

#
pdf("Figures/shannonbyfile.pdf") # Para guardar en PDF
boxplot(shannonentropy.by.cell$Value ~ shannonentropy.by.cell$Name, col=colores,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        cex.names=0.5, cex.axis=0.3, cex.lab=1.2)
dev.off()
pdf("Figures/shanonbysp.pdf")
boxplot(shannonentropy.by.cell$Value ~ shannonentropy.by.cell$species, col=colores,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        names=c( expression(italic("E. bracteata")),expression(italic("E. calcarata")) ,
                 expression(italic("E. coalcomanensis")),expression(italic("E. colligata")),
                 expression(italic("E. conzattii")),expression(italic("E. cymbifera")),
                 expression(italic("E. cyri")), expression(italic("E. diazlunana")),
                 expression(italic("E. finkii")),expression(italic("E. lomelli")),
                 expression(italic("E. peritropoides")),  expression(italic("E. personata")),
                 expression(italic("E. tehuacana")),expression(italic("E. tithymaloides")),
                 expression("L-systemM"),expression("L-systemX"), expression("probBetaLsystem"),
                 expression("probL-system"),expression("RayL-system")),
        cex.names=0.5, cex.axis=0.8, cex.lab=1.2)
dev.off()
pdf("Figures/shanonbyhabit.pdf") # Para guardar en PDF
boxplot(shannonentropy.by.cell$Value ~ shannonentropy.by.cell$habit,col=colores,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4)
dev.off()

pdf("Figures/entropybyfile.pdf")
par(mfrow = c(2,1))
z=1
for(i in files) {
  temp <- subset(shannonentropy.by.cell, shannonentropy.by.cell$Name == i)
  plot(temp$Value, type ="l", xlab = "Cell file position", ylab="Shannon-Entropy value",
       main= i, col= my_palette[z], lwd = 2)
  z = z +1
}
dev.off()
######
#Small subset
shannon_subset<-subset(shannonentropy.by.cell, Name == "probLsystem_edited_cells.txt" | Name=="974_edited_cells.txt" |
                         Name == "EPM10_edited_cells.txt")
shannon_subset <- transform(shannon_subset, Sequence=ave(seq_along(species),
                                                         species, FUN=seq_along))

#Shannon
pdf("Figures/shannon_bycell.pdf")
ggplot(shannon_subset, aes(x = Sequence,y = Value,color = species)) +
  geom_line()
dev.off()
#######
##Remove rlsystem, and context sensitive lsystems
shannonentropy.by.cell_subset <-subset(shannonentropy.by.cell, shannonentropy.by.cell$species != "RayL-system")
shannonentropy.by.cell_subset <-subset(shannonentropy.by.cell_subset,
                                       shannonentropy.by.cell_subset$species != "L-systemMesic")
shannonentropy.by.cell_subset <-subset(shannonentropy.by.cell_subset,
                                       shannonentropy.by.cell_subset$species != "L-systemXeric")
shannonentropy.by.cell_subset$species <- factor(shannonentropy.by.cell_subset$species)
#
shannon.cellsp.lm <-lm(shannonentropy.by.cell_subset$Value ~
                         factor(shannonentropy.by.cell_subset$species))
emm.shannon <- emmeans(shannon.cellsp.lm, specs = pairwise ~ species,adjust="tukey")
emm.shannon$contrasts
multcomp::cld(emm.shannon$emmeans, alpha = 0.10, Letters=LETTERS)
#
pdf("Figures/shanonbysp_subset.pdf")
boxplot(shannonentropy.by.cell_subset$Value ~ shannonentropy.by.cell_subset$species,
        las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        names=c( expression(italic("E. bracteata")),expression(italic("E. calcarata")) ,
                 expression(italic("E. coalcomanensis")),expression(italic("E. colligata")),
                 expression(italic("E. conzattii")),expression(italic("E. cymbifera")),
                 expression(italic("E. cyri")), expression(italic("E. diazlunana")),
                 expression(italic("E. finkii")),expression(italic("E. lomelli")),
                 expression(italic("E. peritropoides")),  expression(italic("E. personata")),
                 expression(italic("E. tehuacana")),expression(italic("E. tithymaloides")),
                 expression("probBetaLsystem"),expression("probL-system")),
        cex.names=0.5, cex.axis=0.8, cex.lab=1.2)
dev.off()
#######Checkout files with ray cells########
shannonentropy.by.cell.with.ray <- read.csv("Data/shannonentropywithR.csv", row.names=1)

#Add species and habit
match.id <- match(shannonentropy.by.cell.with.ray$Name,species.id$files)
shannonentropy.by.cell.with.ray$species <- as.character(species.id$species[match.id])
shannonentropy.by.cell.with.ray$habit <- as.character(species.id$habit[match.id])

boxplot(shannonentropy.by.cell.with.ray$Value ~ shannonentropy.by.cell.with.ray$Name,
        col=colores,las=1.5, xlab="Average Cell-file Shannon-Entropy",ylab = NULL, cex=0.4,
        cex.names=0.5, cex.axis=0.3, cex.lab=1.2)

plot(shannonentropy.by.cell.with.ray$Value[shannonentropy.by.cell.with.ray$Name=="974_edited_cells.txt"],
     type ="l", xlab="Compression distance")
lines(shannonentropy.by.cell$Value[shannonentropy.by.cell$Name=="974_edited_cells.txt"],
      type ="l", col =my_palette[2])
######Check out Shanon by window #####
shannonentropy.by.window <- read.csv("Data/shannonentropy_window.csv", row.names=1)
library(dplyr)
shannonentropy.by.window <- shannonentropy.by.window %>%
  group_by(Sample,Lineage) %>% mutate(id = row_number())

write.csv(shannonentropy.by.window, "Data/shannonentropy_window_id.csv")
window947 <- subset(shannonentropy.by.window, Sample == "974_edited_cells.txt")
window883 <- subset(shannonentropy.by.window, Sample == "883_edited_cells.txt")
window883 <- subset(shannonentropy.by.window, Sample == "883_edited_cells.txt")
windowepm10 <- subset(shannonentropy.by.window, Sample == "EPM10_edited_cells.txt")
prbeta <-subset(shannonentropy.by.window, Sample == "probLsystembeta_edited_cells.txt")
#PLot every sample
samples<-factor(shannonentropy.by.window$Sample)
for (i in samples){
  plot(shannonentropy.by.window$id[shannonentropy.by.window==i],
       shannonentropy.by.window$Value[shannonentropy.by.window==i],type = "l",lwd = 3)
}
plot(window947$id,window947$Value,type = "l")
plot(window947$id[window947$Lineage==26],window947$Value[window947$Lineage==26],type = "l",lwd = 2,
     xlab="Window position",ylab = "Shannon Value")

lineage<-factor(windowepm10$Lineage)
for (i in lineage){
  points(windowepm10$id[windowepm10$Lineage==i], windowepm10$Value[windowepm10$Lineage==i],
         col="green",type = "l",lwd = 3,
         col=pal1)
}
plot(windowepm10$id,windowepm10$Value, type = "l",col=pal1)
plot(window883$id,window883$Value,type = "l",lwd = 3)


pdf("Figures/shannon_window947.pdf")
plot(window947$id[window947$Lineage==2],window947$Value[window947$Lineage==2],type = "l",lwd = 3)
points(window947$id[window947$Lineage==3],window947$Value[window947$Lineage==3], col="green",type = "l",lwd = 3)
points(window947$id[window947$Lineage==4],window947$Value[window947$Lineage==4], col="red",type = "l",lwd = 3)
points(window947$id[window947$Lineage==5],window947$Value[window947$Lineage==5], col=pal1[1],type = "l",lwd = 3)
points(window947$id[window947$Lineage==6],window947$Value[window947$Lineage==6], col=pal1[2],type = "l",lwd = 3)
#points(window947$id[window947$Lineage==7],window947$Value[window947$Lineage==7], col=pal1[3],pch=19)
#points(window947$id[window947$Lineage==8],window947$Value[window947$Lineage==8], col=pal2[1],pch=19)
#points(window947$id[window947$Lineage==9],window947$Value[window947$Lineage==9], col=pal2[2],pch=19)
#points(window947$id[window947$Lineage==10],window947$Value[window947$Lineage==10], col=pal2[3],pch=19)
dev.off()

pdf("Figures/shannon_window883.pdf")
plot(window883$id[window883$Lineage==3],window883$Value[window883$Lineage==3],type = "l",lwd = 3)
points(window883$id[window883$Lineage==2],window883$Value[window883$Lineage==2], col="green",type = "l",lwd = 3)
points(window883$id[window883$Lineage==6],window883$Value[window883$Lineage==6], col="red",type = "l",lwd = 3)
points(window883$id[window883$Lineage==4],window883$Value[window883$Lineage==4], col=pal1[1],type = "l",lwd = 3)
points(window883$id[window883$Lineage==7],window883$Value[window883$Lineage==7], col=pal1[2],type = "l",lwd = 3)
dev.off()

pdf("Figures/probetaS.pdf")
plot(prbeta$id[prbeta$Lineage==6],prbeta$Value[prbeta$Lineage==6],type = "l",lwd = 3)
points(prbeta$id[prbeta$Lineage==2],prbeta$Value[prbeta$Lineage==2], col="green",type = "l",lwd = 3)
points(prbeta$id[prbeta$Lineage==1],prbeta$Value[prbeta$Lineage==1], col="red",type = "l",lwd = 3)
points(prbeta$id[prbeta$Lineage==4],prbeta$Value[prbeta$Lineage==4], col=pal1[1],type = "l",lwd = 3)
points(prbeta$id[prbeta$Lineage==7],prbeta$Value[prbeta$Lineage==7], col=pal1[2],type = "l",lwd = 3)
dev.off()

summary(cell_lengths)

ggplot(window947, aes(x= id, y = Value, color=Lineage)) +
  geom_point()

ggplot(lempel.by.cell, aes(x = Name, y = Value, color = species)) + geom_boxplot()

#####Trying to get total of words per length#####
wc.2 <- as.data.frame(wc.2)
wc.2$sample <- rownames(wc.2)
wc.2<-melt(wc.2)
#Create table counting number of words per sample
wc.2 <-subset(wc.2, wc.2$value != 0) #First remove words with 0 counts
nwords2 <- wc.2 %>% group_by(sample) %>% tally()
#
wc.3 <- as.data.frame(wc.3)
wc.3$sample <- rownames(wc.3)
wc.3 <- melt(wc.3)
#remove 0 and create table
wc.3 <-subset(wc.3, wc.3$value != 0)
nwords3 <- wc.3 %>% group_by(sample) %>% tally()
#
wc.4 <- as.data.frame(wc.4)
wc.4$sample <- rownames(wc.4)
wc.4 <- melt(wc.4)
#remove 0 counts and table
wc.4 <- subset(wc.4, wc.4$value != 0)
nwords4 <- wc.4 %>% group_by(sample) %>% tally()
#
wc.5 <- as.data.frame(wc.5)
wc.5$sample <- rownames(wc.5)
wc.5 <- melt(wc.5)
#remove 0s
wc.5 <- subset(wc.5, wc.5$value != 0)
nwords5 <- wc.5 %>% group_by(sample) %>% tally()
#
wc.6 <- as.data.frame(wc.6)
wc.6$sample <- rownames(wc.6)
wc.6 <- melt(wc.6)
#remove 0s
wc.6 <- subset(wc.6, wc.6$value != 0)
nwords6 <- wc.6 %>% group_by(sample) %>% tally()
#
wc.7 <- as.data.frame(wc.7)
wc.7$sample <- rownames(wc.7)
wc.7 <- melt(wc.7)
#remove 0s
wc.7 <- subset(wc.7, wc.7$value != 0)
nwords7 <- wc.7 %>% group_by(sample) %>% tally()
#
wc.8 <- as.data.frame(wc.8)
wc.8$sample <- rownames(wc.8)
wc.8 <- melt(wc.8)
#remove 0s
wc.8 <- subset(wc.8, wc.8$value != 0)
nwords8 <- wc.8 %>% group_by(sample) %>% tally()
#
wc.9 <- as.data.frame(wc.9)
wc.9$sample <- rownames(wc.9)
wc.9 <- melt(wc.9)
#remove 0s
wc.9 <- subset(wc.9, wc.9$value != 0)
nwords9 <- wc.9 %>% group_by(sample) %>% tally()
####
wc.10 <- as.data.frame(wc.10)
wc.10$sample <- rownames(wc.10)
wc.10 <- melt(wc.10)
######
wc.10 <- subset(wc.10, wc.10$value != 0)
nwords10 <- wc.10 %>% group_by(sample) %>% tally()
#11 and 12
wc.11 <- as.data.frame(wc.11)
wc.11$sample <- rownames(wc.11)
wc.11 <- melt(wc.11)
######
wc.11 <- subset(wc.11, wc.11$value != 0)
nwords11 <- wc.11 %>% group_by(sample) %>% tally()
#12
wc.12 <- as.data.frame(wc.12)
wc.12$sample <- rownames(wc.12)
wc.12 <- melt(wc.12)
######
wc.12 <- subset(wc.12, wc.12$value != 0)
nwords12 <- wc.12 %>% group_by(sample) %>% tally()
#13
wc.13 <- as.data.frame(wc.13)
wc.13$sample <- rownames(wc.13)
wc.13 <- melt(wc.13)
######
wc.13 <- subset(wc.13, wc.13$value != 0)
nwords13 <- wc.13 %>% group_by(sample) %>% tally()
#14
wc.14 <- as.data.frame(wc.14)
wc.14$sample <- rownames(wc.14)
wc.14 <- melt(wc.14)
######
wc.14 <- subset(wc.14, wc.14$value != 0)
nwords14 <- wc.14 %>% group_by(sample) %>% tally()

##Making common data frame
nwords <- rbind(cbind(nwords2,wordLength =rep(2,nrow(nwords2))),
                cbind(nwords3,wordLength =rep(3,nrow(nwords3))),
                cbind(nwords4,wordLength =rep(4,nrow(nwords4))),
                cbind(nwords5,wordLength =rep(5,nrow(nwords5))),
                cbind(nwords6,wordLength =rep(6,nrow(nwords6))),
                cbind(nwords7,wordLength =rep(7,nrow(nwords7))),
                cbind(nwords8,wordLength =rep(8,nrow(nwords8))),
                cbind(nwords9,wordLength =rep(9,nrow(nwords9))),
                cbind(nwords10,wordLength =rep(10,nrow(nwords10))),
                cbind(nwords11,wordLength =rep(11,nrow(nwords11))),
                cbind(nwords12,wordLength =rep(12,nrow(nwords12))),
                cbind(nwords13,wordLength =rep(13,nrow(nwords13))),
                cbind(nwords14,wordLength =rep(14,nrow(nwords14))))
#
rm(nwords10,nwords2,nwords3,nwords4, nwords5,nwords6,
   nwords7, nwords8, nwords9,nwords11,nwords12,nwords13,nwords14)
#
plot(log10(nwords$n) ~ nwords$wordLength)
nwords <- subset(nwords, nwords$sample != "contextmesicLsystem_edited_cells.txt" &
                   nwords$sample != "contextxericLsystem_edited_cells.txt" &
                   nwords$sample != "Ray_Lsystem_edited_cells.txt")
matcher <- match(nwords$sample, species.id$files)
nwords$species <- species.id$species[matcher]
nwords$habit <- species.id$habit[matcher]

nwordslog <-nwords %>% ggplot(aes(x= wordLength, y= log10(n), group = sample, color= habit)) +
  geom_line(size = 1)
nwordsnormal <-nwords %>% ggplot(aes(x= wordLength, y= n, group = sample, color= habit)) +
  geom_line(size = 1)
ggsave("Figures/nwords_log.pdf",nwordslog, device = "pdf")
ggsave("Figures/nwords_notlog.pdf",nwordsnormal, device = "pdf")
ggsave("Figures/nwords_log.png",nwordslog, device = "png")
ggsave("Figures/nwords_notlog.png",nwordsnormal, device = "png")
#

#remove specific factors from plot
nwords.temp <- subset(nwords, nwords$habit != "probL-system" &
                        nwords$habit != "probetaL-system")
nwords_sp <- nwords.temp %>% ggplot(aes(x= wordLength, y= n, group = sample, color= species)) +
  geom_line(size = 1)
nwords_sphabit <- nwords.temp %>% ggplot(aes(x= wordLength, y= n, group = sample, color= habit)) +
  geom_line(size = 1)
ggsave("Figures/nwords_sp.pdf",nwords_sp, device = "pdf")
ggsave("Figures/nwords_sphabit.pdf",nwords_sphabit, device = "pdf")
ggsave("Figures/nwords_sp.png",nwords_sp, device = "png")
ggsave("Figures/nwords_sphabit.png",nwords_sphabit, device = "png")

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
rarecurve(wc.5.portreintamil, col = "blue")
####
wc.all.5 <- as.matrix(read.csv("Data/word_counts_all/wordcounts5.csv", row.names = 1))
wc.all.5[is.na(wc.all.5)] <- 0
#
wc.all.10 <- as.matrix(read.csv("Data/word_counts_all/wordcounts10.csv", row.names = 1))
wc.all.10[is.na(wc.all.10)] <- 0

