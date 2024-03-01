lib <- c("ape", "bite", "rgeos", "sp", "rgdal")
sapply(lib, library, character.only = T)
source("Scripts/mcmc_ADIM.R")

##### Darwin's finches distribution
shp <- readOGR("Application/Data/", "data_0")

tree <- read.tree("Application/Data/Coerebinae.tre")

for(sp in gsub("_", " ", tree$tip.label)){
  png(sprintf("Application/Figures/map_%s.png", sp), width = 40, height = 25, unit = "cm", res = 400)
  plot(shp[which(shp@data$BINOMIAL %in% gsub("_", " ", tree$tip.label)),], main = sp)
  plot(shp[which(shp@data$BINOMIAL %in% sp),], col = "red", add = T)
  dev.off()
}