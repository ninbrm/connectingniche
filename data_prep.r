###XXXXXXXXXXXXXXXXXXXXXX###
### Prepare the data sets ###
###XXXXXXXXXXXXXXXXXXXXXX###

options("rgdal_show_exportToProj4_warnings"="none")

library(sf)
library(rgdal)
library(raster)

areas <- as(st_read("data/reindeer_areas.shp"), 'Spatial')
areas@proj4string<- CRS("+init=epsg:25832")
areas <- spTransform(areas, CRS("+init=epsg:25833"))
plot(areas)

rastt <- brick(raster("data/rsf_sum.tif"), 
          raster("data/ssf_sum.tif"))

rastt <- crop(rastt, areas)
plot(rastt[[1]])
plot(areas, add=T)

rastt <- aggregate(rastt, fact=2)

plot(rastt[[1]])
plot(areas, add=T)

# * Mask to the areas
rastt <- mask(rastt, areas)
plot(rastt)

# * Standardize
rastt <- mask(rastt, areas)

hist(rastt[[1]])
quantile(rastt[[1]], probs=c(0.95, 0.99, 0.999))
hist(rastt[[2]])
quantile(rastt[[2]], probs=c(0.95, 0.99, 0.999))

values(rastt[[1]]) <- ifelse(values(rastt[[1]])>quantile(rastt[[1]], probs=c(0.999)), 1, values(rastt[[1]])/quantile(rastt[[1]], probs=c(0.999)))
values(rastt[[2]]) <- ifelse(values(rastt[[2]])>quantile(rastt[[2]], probs=c(0.999)), 1, values(rastt[[2]])/quantile(rastt[[2]], probs=c(0.999)))

plot(rastt)


# * Clip each area
area_names <- areas$name_area

i=1
for (i in c(1:length(area_names))){
  areaT <- areas[areas$name_area==area_names[i],]
  rastT <- crop(rastt, areaT)
  rastT <- mask(rastT, areaT)
  
  plot(rastT[[1]])
  
  writeRaster(rastT[[1]], filename=paste0("output/rsf_sum_area", i, ".asc"), NAflag=-999, overwrite=T)
  writeRaster(rastT[[2]], filename=paste0("output/ssf_sum_area", i, ".asc"), NAflag=-999, overwrite=T)
}