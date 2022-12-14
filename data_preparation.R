#Script to read and prepare the occurence and environmental data
library(rgbif)
library(CoordinateCleaner)
library(maps)
library(dplyr)
library(sp)
library(terra)
library(mapview)
library(ggplot2)
library(DT)
library(readr)

##Occurrence data
-------------------------------------------------------------------

#Species_list contains the scientific names of all chosen species
species_list <- read.csv2("Species_list.csv")
myspecies <- c(species_list$Species)

#Counting available Occurrence points for species
for (p in 1:length(myspecies)) {
  count <- occ_search(scientificName = myspecies[p], 
                      hasCoordinate = TRUE, 
                      limit = 0)
  species_list$Occurences[p] <- count$meta$count
}
rm(p, count)

#Loading Occurrences (data = raw)
data <- read_delim("C:/Users/luca4/OneDrive/UniUnterlagen/TUM/7.Semester/BA/Daten/gbif_data/data.csv", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)

#Cleaning data points
##Removing records of absence or zero abundance
absence_rows <- which(data$individualCount == 0 |
                      data$occurrenceStatus %in% 
                      c("absent", "Absent", "ABSENT", "ausente", "Ausente", "AUSENTE"))
length(absence_rows)
if (length(absence_rows) > 0) {
  species_data <- data[-absence_rows, ]
}

##Removing invalid/impossible
species_data <- cc_val(species_data,
                  lon = "decimalLongitude",
                  lat = "decimalLatitude",
                  value = "clean",
                  verbose = TRUE)

##Removing falsely georeferenced records
species_data <- clean_coordinates(species_data,
                                  lon = "decimalLongitude",
                                  lat = "decimalLatitude",
                                  species = NULL,
                                  tests = c("capitals", "centroids",
                                            "gbif", "institutions", #"outliers", 
                                            "seas", "zeros"),
                                  value = "clean")

##Removing points with uncertainty > grid cell size of predictor variable 
uncertain <- which(species_data[, "coordinateUncertaintyInMeters"] > 5000 | is.na(species_data[, "coordinateUncertaintyInMeters"]))
length(uncertain)
if (length(uncertain) > 0) {
  species_data <- species_data[-uncertain, ]
}

##Removing unnecessary variables
species_data <- species_data[,c(1,2,3,10,13,14,16,17,18,19,20,
                                21,22,23,24,26,28,30,35,36)]

##Counting cleaned data
for (a in 1:length(myspecies)) {
  species_list$Cleaned[a] <- length(which(species_data[, "species"] == myspecies[a]))
}

##Creating seperate datasets
for (b in 1:length(myspecies)) {
  species_rows <- which(species_data[, "species"] == myspecies[b])
  species_temp <- species_data[species_rows,]
  assign(myspecies[b], species_temp)
  print(myspecies[b])
}
rm(b, species_rows, species_temp)

##Adding climate, disturbance and forestcover parameters
for (c in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[c])
  species_data_temp$ID <- c(1:nrow(species_data_temp))
  ###Vectorizing Data
  species_data_temp_vec <- vect(species_data_temp, 
                           geom=c("decimalLongitude", "decimalLatitude"), 
                           crs="+init=epsg:4326")
  ###Extracting Predictor-data for Occurence point
  species_dib_clim_fc <- terra::extract(dib_clim_fc, species_data_temp_vec)
  ###Joining Predictor data and occurrences
  species_data_temp <- full_join(species_dib_clim_fc, species_data_temp, by="ID")
  assign(myspecies[c], species_data_temp)
  print(myspecies[c])
}
rm(c, species_data_temp_vec, species_dib_clim_fc, species_data_temp)

##Finding rows with all NA
for (d in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[d])
  dib_NA_rows_all <- which(is.na(species_data_temp$dib1) & 
                             is.na(species_data_temp$dib2) & 
                             is.na(species_data_temp$dib3) & 
                             is.na(species_data_temp$dib4))
  length(dib_NA_rows_all)
  if (length(dib_NA_rows_all) > 0) {
    species_data_temp <- species_data_temp[-dib_NA_rows_all, ]
  }
  assign(myspecies[d] ,species_data_temp)
  print(myspecies[d])
}
rm(d, species_data_temp, dib_NA_rows_all)

##Counting cleaned data
for (e in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[e])
  species_list$Cleaned_Dib[e] <- nrow(species_data_temp)
}
rm(e, species_data_temp)

##Cleaning for distance: Removing all occurrence points with <5km distance to each other
dsg <- function(data.coor,mindist){
  z <- cbind(data.coor, 1:nrow(data.coor))
  sample.out <- sample(1:nrow(data.coor), nrow(data.coor))
  xy <- data.coor[sample.out,]
  i <- 1
  repeat{
    keep <- z[i,]
    zx <- which(z[,1]>(z[i,1]-mindist) & z[,1] < (z[i,1]+mindist))
    zy <- which(z[,2]>(z[i,2]-mindist) & z[,2] < (z[i,2]+mindist))
    out <- zy[na.exclude(match(zx,zy))]
    if (length(out) > 1){
      z <- rbind(keep,z[-out,])}
    #print(i)
    i<-i + 1
    if (i >= nrow(z)){break}
  }
  return (z[,3])
}
for (f in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[f])
  subsample <- dsg(species_data_temp[,c("decimalLatitude", "decimalLongitude")], 0.04166667)
  species_data_temp <- species_data_temp[subsample,]
  assign(myspecies[f] ,species_data_temp)
  print(myspecies[f])
}
rm(f, species_data_temp)

##Counting distance-cleaned data
for (g in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[g])
  species_list$D_Cleaned[g] <- nrow(species_data_temp)
}
rm(g, species_data_temp)

##Environmental data
-------------------------------------------------------------------
library(rgdal)
library(terra)

#Stacking bioclim-data
bioclim <- rast("bioclim_europe.tif", lyrs=1)
for (h in 2:19) {
  bioclim <- c(bioclim, rast("bioclim_europe.tif", lyrs=h))
}
names(bioclim) <- paste0("bio", seq(1, 19, 1))

#Stacking disturbance-data
disturbance <- rast("disturbance_metrics.tif", lyrs=1)
for (i in 2:4) {
  disturbance <- c(disturbance, rast("disturbance_metrics.tif", lyrs=i))
}
names(disturbance) <- paste0("dib", seq(1, 4, 1))

#Loading forestcover data
forestcover <- rast("forestcover.tif")

#Combining all predictor-data
dib_clim_fc <- c(disturbance, bioclim, forestcover)





