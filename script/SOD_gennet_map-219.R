#------ loading libraries
library(reshape2)
library(ggplot2)
library(ape)
library(igraph)
#
library(sf)
library(rgeos)
#
library("rnaturalearth")
library(rnaturalearthdata)
library("maps")
library("tools")

#--- setting working directory
setwd("/Users/ricardoi/git_local/genedis_networks/")

#----- loading SOD metadata
data <- read.csv("data/hazel_population_data.csv")
data.wgs <- read.csv("data/hazel_population_data.wgs.csv")
#----- loading genetic distance data
gd <- read.csv("~/Box/OSU/P_ramorum/data/na1_gendists.txt", header = T, sep = " ")
colnames(gd) <- rownames(gd)
gd[1:4,1:4]


#' subsetting data sets
dat <- data[data$ID %in% data.wgs$ID,]
dat
#' fornattign matrix table to long format
gd1 <- as.data.frame(gd)#$`7612-A1`)
rownames(gd1) <- rownames(gd)
gdnet <- melt(as.matrix(gd1))
# renaming columns
colnames(gdnet) <- c("From", "To", "value")
head(gdnet)

# data manipulation - setting up arbitrary tresholds
gdnet0 <- gdnet[!rowSums(gdnet[-c(1:2)] == 0) >= 1,] #removing 0's
gdnet <- gdnet0[gdnet0[-c(1:2)] < 1000,] # setting up the treshold
head(gdnet)
# hist(gdnet$value, breaks= 100)

#----- Network
x = gdnet
igraph<-graph_from_data_frame(x, directed = FALSE)
rbPal <- colorRampPalette(c("green", "yellow", "red"))
counPal <- colorRampPalette(c("red", "yellow", "blue", "white", "brown"), bias = 1)
V(igraph)$size=5
V(igraph)$xx <- as.numeric(as.factor(V(igraph))) # make the categories of x into numeric values for color ramp
# V(igraph)$color <- counPal(10)[cut(as.numeric(V(igraph)$xx),breaks = 10)]
CounColor <- unique(cbind(V(igraph)$Country, V(igraph)$color))
E(igraph)$xx <- as.numeric(unlist(E(igraph))) # make the categories of x into numeric values for color ramp
E(igraph)$color <- rbPal(10)[cut(as.numeric(E(igraph)$xx),breaks = 10)]

plot(igraph,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_kk)
# legend(x= 0.9, y= -0.9, CounColor[,1], pch=21,  col="black", pt.bg=CounColor[,2], pt.cex=2,cex=.8, bty="n", ncol=1)


#----- map
# # Select world data
# world <- ne_countries(scale = "medium", returnclass = "sf")
# class(world)
#
# ggplot(data = world) +
#   geom_sf() +
#   geom_point(dat = dat, aes(x = Lon, y = Lat, col = Lin), size = 2,
#              shape = 1) +
#   coord_sf(xlim = c(-124.45, -124.2), ylim = c(42.05, 42.4), expand = T)

## put it on a map
library(rworldmap)
library(rworldxtra)
worldmap <- getMap(resolution = "high") #grab the world map
NorthAmerica <- worldmap[which(worldmap$REGION == "North America"),] # grab north america

plot(igraph,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_kk)


plot(NorthAmerica, xlim = c(-124.5, -124.3), ylim = c(42.1, 42.325))
lo <- as.matrix(dat[,c(6,5)]) ## set the layout to our coordinates
plot(igraph, vertex.size = 0.5, edge.arrow.size =.05, vertex.label.cex= 0.3,
     vertex.label.color = "", edge.width = 0.2, edge.curved = TRUE, layout = lo,
     xlim = c(-124.45, -124.2), ylim = c(42.05, 42.4), rescale = FALSE, add = TRUE)
      ## very important to set layout to "lo" (or whatever you named it above), rescale = FALSE, add = TRUE

