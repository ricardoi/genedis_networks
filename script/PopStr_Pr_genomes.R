#'@author Ricardo I Alcala
#'@title Genetic Distance - Network Analysis
#'@description 
#'

#- setting working directory
setwd("~/Box/OSU/P_ramorum/data/")

#----functions

#------ loading libraries
library(reshape2)
library(ggplot2)
library(ape)
library(igraph)

#------
gd <- read.csv("na1_gendists.txt", header = T, sep = " ")
colnames(gd) <- rownames(gd)
class(gd)
gd[1:4,1:4]

gd1 <- as.data.frame(gd$`7612-A1`)
rownames(gd1) <- rownames(gd)
gdnet <- melt(as.matrix(gd1)) 
colnames(gdnet) <- c("From", "To", "value")
head(gdnet)
gdnet$To <- rep("7612-A",dim(gdnet)[1],)
#
gdnet0 <- gdnet[!rowSums(gdnet[-c(1:2)] == 0) >= 1,] #removing 0's
gdnet <- gdnet0[gdnet0[-c(1:2)] > 85,] # setting up the treshold
head(gdnet)


#----- Network
x = gdnet
igraph<-graph_from_data_frame(x, directed = FALSE)
rbPal <- colorRampPalette(c("green", "yellow", "red")) 
counPal <- colorRampPalette(c("red", "yellow", "blue", "white", "brown"), bias = 1)
V(igraph)$size=5
# V(igraph)$xx <- as.numeric(as.factor(V(igraph))) # make the categories of x into numeric values for color ramp
# V(igraph)$color <- counPal(10)[cut(as.numeric(V(igraph)$xx),breaks = 10)]
CounColor <- unique(cbind(V(igraph)$Country, V(igraph)$color))
E(igraph)$xx <- as.numeric(unlist(E(igraph))) # make the categories of x into numeric values for color ramp
E(igraph)$color <- rbPal(10)[cut(as.numeric(E(igraph)$xx),breaks = 10)]


plot(igraph,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_kk)
# legend(x= 0.9, y= -0.9, CounColor[,1], pch=21,  col="black", pt.bg=CounColor[,2], pt.cex=2,cex=.8, bty="n", ncol=1)
