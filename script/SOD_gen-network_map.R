library(sf)
library(rgeos)
library(ggplot2)
#
library("rnaturalearth")
library(rnaturalearthdata)
library("maps")
library("tools")


setwd("/Users/ricardoi/git_local/genedis_networks/")


# loading SOD data
data <- read.csv("data/hazel_population_data.csv")
data.wgs <- read.csv("data/hazel_population_data.wgs.csv")

dim(data)
dim(data.wgs)

dat <- data[data$ID %in% data.wgs$ID,]
dat

#----- map
# Select world data
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf() +
  geom_point(dat = dat, aes(x = Lon, y = Lat, col = Lin), size = 2,
             shape = 1) +
  # geom_sf(data = counties, fill = NA ) +
  coord_sf(xlim = c(-124.45, -124.2), ylim = c(42.05, 42.4), expand = T)
