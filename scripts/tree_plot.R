library("ape")
library("Biostrings")
library("ggplot2")
library("tidyverse")
library("ggtree")
library("phytools")
## try http:// if https:// URLs are not supported

tree <- read.tree("../Data/Pedilanthus/tree.txt")
ggtree(tree,size=1) + geom_tiplab() + 
  geom_cladelabel(node=17, label="Pedilanthus clade", color="red")+
  geom_cladelabel(node=24, label="Xéricas", 
                  color="red2", offset=.8)  
  
ggtree(tree)+geom_text(aes(label=node), hjust=-.3)

geom_cladelabel(node=10, label="Some random clade", 
                color="red2", offset=.8, align=TRUE)

colors = c("#FF000060","#0000ff60","#551a8b60","#00FF0033","#FFFF0033")

plot.phylo(tree)
rect(0.09,38.5,2,51, col = colors[1],border = F) #group 1
rect(0.6,37.8,2.4,25.5,col =colors[2],border=F) #group 2
rect(0.4,25.5,2.7,18.5,col = colors[3],border = F) #group 3

tree_p<-read.tree("../Meta/prueba_tree.txt")
x<-as.matrix(read.csv("../Meta/x.csv",row.names=1))[,1]
dotTree(tree_p,x,length=10,ftype="i")

?phytools
dotTree(tree,x,length=10,ftype="i")
nodelabels("Pedilanthus clade", 17,frame = "none", bg = "tomato", font = 3)
?cladelabels
cladelabels(tree=NULL ,"Xéricas", 24)
cladelabels(tree=NULL ,"Mésica", 20)

trans<-read.csv("../Data/Pedilanthus/tabla_transiciones.csv", sep = ",")
trans_filter<-trans[c(1,2,4:9,12:17),c(1,4:12)]
f_f<-as.matrix(trans_filter, row.names=as.character(trans_filter$Especie))[,4]
?as.matrix
dotTree(tree)
phylo.heatmap(tree,f_f,standardize=TRUE)

as.numeric(as.character(trans[,c(4:12)]))
heatmap(as.numeric(trans[,c(4:12)]))

        