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
gd <- read.csv("eu1_gendists_snpsubset.txt", header = T, sep = " ") |>
                as.matrix()
colnames(gd) <- rownames(gd)
gd[1:4,1:4]

gd0 = gd[gd > 0]

hist(gd0, xlim = c(0,0.4))

# long format
gd1 <- as.data.frame(gd)#$`7612-A1`)
# rownames(gd1) <- rownames(gd)
gdnet <- melt(as.matrix(gd1))
colnames(gdnet) <- c("From", "To", "value")
head(gdnet)
# gdnet$To <- rep("7612-A",dim(gdnet)[1],)
#
gdnet0 <- gdnet[!rowSums(gdnet[-c(1:2)] == 0) >= 1,] #removing 0's
gdnet <- gdnet0[gdnet0[-c(1:2)] > 0,] # setting up the treshold
head(gdnet)
hist(gdnet$value, breaks= 100)

#----- Network
x = gdnet
GRPH <- graph_from_data_frame(x, directed = FALSE)
rbPal <- colorRampPalette(c("yellow", "red"))
counPal <- colorRampPalette(c("red", "yellow", "blue", "white", "brown"), bias = 1)
V(GRPH)$size=5
# V(GRPH)$xx <- as.numeric(as.factor(V(GRPH))) # make the categories of x into numeric values for color ramp
# V(igraph)$color <- counPal(10)[cut(as.numeric(V(GRPH)$xx),breaks = 10)]
CounColor <- unique(cbind(V(GRPH)$Country, V(igraph)$color))

E(GRPH)$xx <- x$value # make the categories of x into numeric values for color ramp
# E(igraph)$color <- rbPal(10)[cut(as.numeric(E(igraph)$xx),breaks = 3)]

# Identify isolated nodes
cond  <- E(GRPH)[E(GRPH)$xx <=  0.125 ]
# Remove edges nodes
GRPH <- delete.edges(GRPH, cond)
E(GRPH)$color <- rbPal(10)[cut(as.numeric(E(GRPH)$xx),breaks = 10)]

plot(GRPH,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_kk)
# legend(x= 0.9, y= -0.9, CounColor[,1], pch=21,  col="black", pt.bg=CounColor[,2], pt.cex=2,cex=.8, bty="n", ncol=1)

