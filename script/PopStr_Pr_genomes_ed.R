#'@author Ricardo I Alcala
#'@title Genetic Distance - Network Analysis
#'@description
#'

#- setting working directory
setwd("~/Box/OSU/P_ramorum/data/")

#----functions

#------ loading libraries
library(reshape2)
library(tidyverse)
library(ape)
library(igraph)

library(poppr)
library(castor)

#------
gd <- read.csv("na1_gendists_snpsubset.txt", header = T, sep = " ") |>
    as.matrix()
colnames(gd) <- rownames(gd)
gd[1:4,1:4]

hist(as.matrix(gd))

eu1_mac_tre <- aboot(as.matrix(gdm),
                     sample = 100, showtree = T, tree = "nj")
class(eu1_mac_tre)

all.dis <- get_all_distances_to_root(eu1_mac_tre)
node.dis <- all.dis[(132+1):(132+eu1_mac_tre$Nnode)]

nd0 <- node.dis[node.dis > 0]

hist(nd0, xlab="distance to root", ylab="# nodes", prob=FALSE,
     breaks = 20, xlim = c(0,0.25))

dist.nodes(eu1_mac_tre)|>
  as_tibble()

# long format
gd1 <- as.data.frame(gd)#$`7612-A1`)
rownames(gd1) <- rownames(gd)
gdnet <- melt(as.matrix(gd1))

# comparison by 1
colnames(gdnet) <- c("From", "To", "value")
head(gdnet)
# gdnet$To <- rep("7612-A",dim(gdnet)[1],)
#
gdnet0 <- gdnet[!rowSums(gdnet[-c(1:2)] == 0) >= 1,] #removing 0's
gdnet <- gdnet0[gdnet0[-c(1:2)] < 10000,] # setting up the treshold
head(gdnet)
hist(gdnet$value, breaks= 100)

# cutoff for entire matrix
# gdnet <- gdnet0[gdnet0[-c(1:2)] < 1000,] # setting up the treshold


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

