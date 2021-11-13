#'@author Ricardo I Alcala
#'@title Genetic Distance - Network Analysis
#'@description
#'

#- setting working directory
# setwd("~/Box/OSU/P_ramorum/data/")
setwd("~/git_local/genedis_networks/data/")

#----functions

#------ loading libraries
library(reshape2)
library(tidyverse)
library(ape)
library(igraph)
library(poppr)
library(TreeTools)
library(dendextend)

#------
eu1 <- read.csv("eu1_gendists_snpsubset.txt", header = T, sep = " ") |>
                as.matrix()
eu1_names = rownames(eu1)
eu1_gd <- matrix(as.numeric(eu1), nrow = dim(eu1)[1])
colnames(eu1) <- rownames(eu1) <- names
eu1[1:4,1:4]
class(eu1)
#------
na1 <- read.csv("na1_gendists_snpsubset.txt", header = T, sep = " ") |>
  as.matrix()
na1_names = rownames(na1)
na1_gd <- matrix(as.numeric(na1), nrow = dim(na1)[1])
colnames(na1) <- rownames(na1) <- names
na1[1:4,1:4]
class(na1)
#----- loading SOD metadata
data <- read.csv("hazel_population_data.csv")
dim(data)
data.wgs <- read.csv("hazel_population_data.wgs.csv")
dim(data.wgs)
#' subsetting data sets
dat = data.wgs[data$ID %in% data.wgs$ID,]
dat

# Minimum Spaning Tree
MSTEdges(eu1, plot = T)
MSTEdges(na1, plot = T)
MSTLength(eu1)
MSTLength(na1)

# ploting distances
par(mfrow = c(2, 1))
hist(eu1[eu1 > 0], xlim = c(0,0.35), breaks = 100)
hist(na1[na1 > 0], xlim = c(0,0.35), breaks = 100)
dev.off()
eq <- function(x) (x*1)
# recalculating distances
eu1_mac_tre <- aboot(as.matrix(eu1), dist = eq ,
                     sample = 100, showtree = T, tree = "nj")
# recalculating distances
na1_mac_tre <- aboot(as.matrix(na1), dist = eq ,
                     sample = 100, showtree = T, tree = "nj")

# tree manipulation
library("phylogram")
eu1d <- as.dendrogram(eu1_mac_tre)
eu1.dend <- eu1d |>
  dendextend::get_nodes_attr("members") # node's membership
na1d <- as.dendrogram(na1_mac_tre)
na1.dend <- na1d |>
  dendextend::get_nodes_attr("members") # node's height
par(mfrow=c(2,1))
hist(eu1.dend, breaks= 100, xlim = c(0,75))
abline(v=010, col="red")
hist(na1.dend, breaks= 100, xlim = c(0,75))
abline(v=010, col="red")
dev.off()

library("treeman")

y1 <- dat |>
  subset(Year == 2016 & ID %in% eu1_mac_tre$tip.label)|>
  select(ID, Year)
y2 <- dat |>
      subset(Year == 2017 & ID %in% eu1_mac_tre$tip.label)|>
  select(ID, Year)
y3 <- dat |>
  subset(Year == 2018 & ID %in% eu1_mac_tre$tip.label)|>
  select(ID, Year)
y4 <- dat |>
  subset(Year == 2019 & ID %in% eu1_mac_tre$tip.label)|>
  select(ID, Year)

datos <- dat |>
          subset( ID %in% eu1_mac_tre$tip.label)
cmatrix = matrix(0, ncol = 5, nrow = nrow(datos))
  for (i in 1:nrow(datos)){
        cmatrix[i,1] <- datos[i,1]
        cmatrix[i,2] <- ifelse(datos[i,1] %in% y1$ID, yes = 1, no = 0)
        cmatrix[i,3] <- ifelse(datos[i,1] %in% y2$ID, yes = 1, no = 0)
        cmatrix[i,4] <- ifelse(datos[i,1] %in% y3$ID, yes = 1, no = 0)
        cmatrix[i,5] <- ifelse(datos[i,1] %in% y4$ID, yes = 1, no = 0)
  }
rownames(cmatrix) <- cmatrix[,1]
cmatrix <- cmatrix[,c(2:5)] |>
            t()
cmatrix <- apply(cmatrix, 2, FUN=as.numeric)
rownames(cmatrix) <- c(2016, 2017, 2018, 2019)
head(cmatrix)

eu1_gc <- commplot(cmatrix, eu1_mac_tre, groups=c(1:4), no.margin=FALSE)

# NA1
y1 <- dat |>
  subset(Year == 2001 & ID %in% na1_mac_tre$tip.label)|>
  select(ID, Year)
y2 <- dat |>
  subset(Year == 2002 & ID %in% na1_mac_tre$tip.label)|>
  select(ID, Year)
y3 <- dat |>
  subset(Year == 2003 & ID %in% na1_mac_tre$tip.label)|>
  select(ID, Year)
y4 <- dat |>
  subset(Year == 2004 & ID %in% na1_mac_tre$tip.label)|>
  select(ID, Year)
y5 <- dat |>
  subset(Year == 2005 & ID %in% na1_mac_tre$tip.label)|>
  select(ID, Year)

datos <- dat |>
  subset( ID %in% na1_mac_tre$tip.label)
cmatrix = matrix(0, ncol = 6, nrow = nrow(datos))
for (i in 1:nrow(datos)){
  cmatrix[i,1] <- datos[i,1]
  cmatrix[i,2] <- ifelse(datos[i,1] %in% y1$ID, yes = 1, no = 0)
  cmatrix[i,3] <- ifelse(datos[i,1] %in% y2$ID, yes = 1, no = 0)
  cmatrix[i,4] <- ifelse(datos[i,1] %in% y3$ID, yes = 1, no = 0)
  cmatrix[i,5] <- ifelse(datos[i,1] %in% y4$ID, yes = 1, no = 0)
  cmatrix[i,6] <- ifelse(datos[i,1] %in% y5$ID, yes = 1, no = 0)
}
rownames(cmatrix) <- cmatrix[,1]
cmatrix <- cmatrix[,c(2:6)] |>
  t()
cmatrix <- apply(cmatrix, 2, FUN=as.numeric)
rownames(cmatrix) <- c(2001, 2002, 2003, 2004, 2005)
head(cmatrix)
library(MoreTreeTools)

tree <- compute.brlen(na1_mac_tre, method = "Grafen")
tree$edge.length <- tree$edge.length/getSize(tree, "rtt")
plot(tree)
commplot(cmatrix, na1_mac_tre, groups=c(1:5), no.margin=TRUE)

tree_tm <- as(na1_mac_tre, 'TreeMan')
c1_ids <- colnames(cmatrix)[cmatrix[1, ] > 0]
c2_ids <- colnames(cmatrix)[cmatrix[2, ] > 0]
c3_ids <- colnames(cmatrix)[cmatrix[3, ] > 0]
c4_ids <- colnames(cmatrix)[cmatrix[4, ] > 0]
c5_ids <- colnames(cmatrix)[cmatrix[5, ] > 0]

obs_ovrlp <- calcOvrlp(tree_tm, c1_ids, c2_ids, c3_ids, c4_ids, c5_ids)
iterations <- 99
null <- rep(NA, iterations)

for(i in 1:iterations) {
  # cat('....[', i, ']\n', sep= "")
  null_tips <- sample(tree_tm['tips'], length(c5_ids))
  null[i] <- calcOvrlp(tree_tm, c5_ids, null_tips)
}
p_value <- sum(obs_ovrlp >= null)/iterations
hist(null, main="", xlab="", ylab="", xlim=c(0,1), ylim=c(0,20))
abline(v=obs_ovrlp, col="red")
mtext(paste0("P-value: ", signif(p_value, 3)))
cat("P-value: ", signif(p_value, 3), "\n", sep="")


# genetic distance from-to network
eu1_gdnet <- as.data.frame(eu1_mac_tre$edge)
eu1_gdnet$length <- as.numeric(eu1_mac_tre$edge.length)
colnames(eu1_gdnet) <- c("From", "To", "value")

# genetic distance from-to network
na1_gdnet <- as.data.frame(na1_mac_tre$edge)
na1_gdnet$length <- as.numeric(na1_mac_tre$edge.length)
colnames(na1_gdnet) <- c("From", "To", "value")
# ploting genetic distances
par(mfrow=c(2,1))
hist(eu1_gdnet$value, breaks= 100, xlim = c(0,0.125))
abline(v=0.01, col="red")
hist(na1_gdnet$value, breaks= 100, xlim = c(0,0.125))
abline(v=0.01, col="red")
dev.off()
#--------------------------------
#----- Network
x = eu1_gdnet
y = na1_gdnet
GRPHx <- graph_from_data_frame(x, directed = FALSE)
GRPHy <- graph_from_data_frame(y, directed = FALSE)
# adding colors and attributes
rbPal <- colorRampPalette(c("grey", "black"))
counPal <- colorRampPalette(c("red", "yellow", "blue", "white", "brown"), bias = 1)
V(GRPHx)$size=5
V(GRPHy)$size=5

yearx <- dat[which(dat$ID %in% eu1_mac_tre$tip.label),]|>
                  select(Year) |>
                  unlist() |>
                  as.numeric() # make the categories of x into numeric values for color ramp
V(GRPHx)$year <- yearx
V(igraph)$color <- counPal(10)[cut(as.numeric(V(GRPHx)$year),breaks = 10)]
# CounColorx <- unique(cbind(V(GRPHx)$Country, V(GRPHx)$color))
# CounColory <- unique(cbind(V(GRPHy)$Country, V(GRPHy)$color))

E(GRPHx)$xx <- x$value # make the categories of x into numeric values for color ramp
E(GRPHy)$xx <- y$value # make the categories of x into numeric values for color ramp
# E(igraph)$color <- rbPal(10)[cut(as.numeric(E(igraph)$xx),breaks = 3)]

# ploting branch lenght
par(mfrow=c(2,1))
hist(E(GRPHx)$xx, breaks= 100, xlim = c(0,0.02))
abline(v=0.005, col="red")
hist(E(GRPHy)$xx, breaks= 100, xlim = c(0,0.02))
abline(v=0.005, col="red")
dev.off()
# Select conditionals
condx  <- E(GRPHx)[E(GRPHx)$xx >=  0.005]
condy  <- E(GRPHy)[E(GRPHy)$xx >=  0.005]
# Remove edges nodes
GRPHx <- delete.edges(GRPHx, condx)
GRPHy <- delete.edges(GRPHy, condy)

E(GRPHx)$color <- rbPal(10)[cut(as.numeric(E(GRPHx)$xx),breaks = 10)]
E(GRPHy)$color <- rbPal(10)[cut(as.numeric(E(GRPHy)$xx),breaks = 10)]

plot(GRPHx,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',
     edge.curved = F, edge.width=1, layout=layout_with_kk)
# legend(x= 0.9, y= -0.9, CounColor[,1], pch=21,  col="black", pt.bg=CounColor[,2], pt.cex=2,cex=.8, bty="n", ncol=1)
plot(GRPHy,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',
     edge.curved = F, edge.width=1, layout=layout_with_kk)
dev.off()

## put it on a map
library(rworldmap)
library(rworldxtra)

worldmap <- getMap(resolution = "high") #grab the world map
NorthAmerica <- worldmap[which(worldmap$REGION == "North America"),] # grab north america


plot(GRPHx,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_kk)

plot(GRPHy,  edge.arrow.size=.05, vertex.label.cex=.3, vertex.label.color='black',edge.curved=T, edge.width=0.2, layout=layout_with_kk)



lox <- dat[which(dat$ID %in% eu1_mac_tre$tip.label),] |>
  select(Lon, Lat)
loy <- dat[which(dat$ID %in% na1_mac_tre$tip.label),] |>
      select(Lon, Lat)
# Network map plot
plot(NorthAmerica, xlim = c(-124.5, -124.25), ylim = c(42.275, 42.35))
plot(GRPHx, vertex.size = 0.5, edge.arrow.size =.05, vertex.label.cex= 0.3,
     vertex.label.color = "", edge.width = 0.2, edge.curved = TRUE, layout = lox,
     xlim = c(-124.45, -124.2), ylim = c(42.1, 42.325), rescale = FALSE, add = TRUE)

plot(NorthAmerica, xlim = c(-124.5, -124.1), ylim = c(42.125, 42.129))
plot(GRPHy, vertex.size = 0.5, edge.arrow.size =.05, vertex.label.cex= 0.3,
     vertex.label.color = "", edge.width = 0.2, edge.curved = TRUE, layout = loy,
     xlim = c(-124.45, -124.2), ylim = c(42.1, 42.2), rescale = FALSE, add = TRUE)

## very important to set layout to "lo" (or whatever you named it above), rescale = FALSE, add = TRUE
GRPHx
