
# Libraries -----
library(tidyr)
library(dplyr)
library(sf)                 # spatial data wrangling
library(rnaturalearth)     # optional; for base map
library(rnaturalearthdata) # optional; for base map
library(ggplot2)            # optional for mapping
library(adegenet)           # compute diversity metrics for genetic data
library(hierfstat)          # compute diversity metrics for genetic data

# Some extra custom functions:
source('other_functions.R')

microsats <- read.csv('D:/科研软件/data/Global_changing_data/macrocourse_EU_data_2024.csv', header = TRUE)
coordinates <- read.csv('D:/科研软件/data/Global_changing_data/EU_americanus_coordinates.csv', header = TRUE)
taxa <- read.csv2("D:/科研软件/data/Global_changing_data/taxonomy.csv")
microsats <- merge(microsats, taxa)

microsats$bats <- ifelse(microsats$order == "chiroptera", 1, 0)
microsats$nonbats <- ifelse(microsats$order == "chiroptera", 0, 1)

head(microsats)
str(microsats)
# Overview of the data:
head(microsats)
str(microsats)
head(coordinates)
str(coordinates)

microsats_bats <- microsats[microsats$bats == 1,]
microsats_nonbats <- microsats[microsats$nonbats == 1,]

microsats_rhodentia <- microsats[microsats$rhodentia == 1,]
microsats_nonrhodentia <- microsats[microsats$nonrhodentia == 1,]

# 1. Mapping samples -----
plot(coordinates$lon, coordinates$lat)
text(coordinates$lon, coordinates$lat, labels = coordinates$pop)
verviewCoordinateReferenceSystems.pdf

coords_spatial <- st_as_sf(coordinates, coords = c('lon', 'lat'), crs = 4326)
coords_spatial <- st_as_sf(microsats, coords = c('lon', 'lat'), crs = 4326)
bats_spatial <- st_as_sf(microsats_bats, coords = c('lon', 'lat'), crs = 4326)
nonbats_spatial <- st_as_sf(microsats_nonbats, coords = c('lon', 'lat'), crs = 4326)

head(coords_spatial)


library(terra)
library(sf)
library(ggplot2)
library(tidyterra)

gdat <- read.csv('D:/科研软件/data/Global_changing_data/macrocourse_EU_data_2024.csv',h=T) %>% 
  st_as_sf(coords=c('lon', 'lat'), crs=4326)     # convert to sf (crs=4326 is the geographic coordinate system WGS84)

# We need an outline of Europe to mask the altitude raster file:
europe <- rnaturalearth::ne_countries(continent = 'Europe', returnclass = 'sf') %>%
  st_make_valid() %>%               # fix any invalid geometries that give errors
  st_crop(st_bbox(gdat)) %>%        # Crop the Europe map to the same extent as the genetic data
  st_union()                        # Dissolve country borders to get 1 outline

# Plot to check:
ggplot(data=europe) + geom_sf()

# Convert to terra format (faster):
eu <- vect(europe)

# Read in altitude file with terra
altitude <- rast('WC_alt.tif') # equal area projection

euproj <- terra::project(eu, crs(altitude))  # Reproject the Europe outline so it correctly aligns with the raster file

altcrop <- terra::crop(altitude, euproj, mask = TRUE)  # Cut out the Europe portion of the altitude map

# Plot:
ggplot() +
  geom_spatraster(data=altcrop) +
  coord_sf(crs = 'ESRI:102014') +
  #scale_fill_continuous(na.value= 'transparent', name = 'altitude (m)') +
  viridis::scale_fill_viridis(option = 'H', na.value= 'transparent', name = 'altitude (m)') +
  theme_void()



#all species
world_map <- ne_countries(scale = "medium", returnclass = 'sf')
ca_eu <- subset(world_map,continent = "Europe"|name == "Turkey")
xlim <- c(-15, 50)
ylim <- c(35, 70)
ggplot() + 
  geom_sf(data = ca_eu) +
  geom_sf(data = coords_spatial, aes(color=allelic_richness)) +
  coord_sf(xlim = xlim, ylim = ylim) +
  scale_colour_gradient(low ="blue2", high = "orange") 
  theme_minimal()

  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = coords_spatial, aes(color=global_fst)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    scale_colour_gradient(low ="blue2", high = "orange")
  theme_minimal()
  
  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = coords_spatial, aes(color=gene_diversity)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    scale_colour_gradient(low ="blue2", high = "orange")
  theme_minimal()
  
# bats
  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = bats_spatial, aes(color=allelic_richness)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    scale_colour_gradient(low ="blue2", high = "orange")
  theme_minimal()
  
  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = bats_spatial, aes(color=global_fst)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    scale_colour_gradient(low ="blue2", high = "orange")
  theme_minimal()
  
  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = bats_spatial, aes(color=gene_diversity)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    scale_colour_gradient(low ="blue2", high = "orange")
  theme_minimal()
  
 
  microsats$rhodentia <- ifelse(microsats$order == "rhodentia", 1, 0)
  microsats$nonrhodentia <- ifelse(microsats$order == "rhodentia", 0, 1)
  
  microsats_rhodentia <- microsats[microsats$rhodentia == 1,]
  microsats_nonrhodentia <- microsats[microsats$nonrhodentia == 1,]
  
  rhodentia_spatial <- st_as_sf(microsats_rhodentia, coords = c('lon', 'lat'), crs = 4326)
  
#rhodentia
  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = rhodentia_spatial, aes(color=allelic_richness)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    scale_colour_gradient(low ="blue2", high = "orange")
  theme_minimal()
  
  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = rhodentia_spatial, aes(color=global_fst)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    scale_colour_gradient(low ="blue2", high = "orange")
  theme_minimal()
  
  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = rhodentia_spatial, aes(color=gene_diversity)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    scale_colour_gradient(low ="blue2", high = "orange")
  

  
#species
#bats
  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = bats_spatial, aes(color=species)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    theme_minimal()
  
#rhodentia
  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = rhodentia_spatial, aes(color=species)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    theme_minimal()
  
#carnivora
  microsats$carnivora <- ifelse(microsats$order == "carnivora", 1, 0)
  microsats_carnivora <- microsats[microsats$carnivora == 1,]
  carnivora_spatial <- st_as_sf(microsats_carnivora, coords = c('lon', 'lat'), crs = 4326)
  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = carnivora_spatial, aes(color=species)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    theme_minimal()
  
#artiodactyla
  microsats$artiodactyla <- ifelse(microsats$order == "artiodactyla", 1, 0)
  microsats_artiodactyla <- microsats[microsats$artiodactyla == 1,]
  artiodactyla_spatial <- st_as_sf(microsats_artiodactyla, coords = c('lon', 'lat'), crs = 4326)
  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = artiodactyla_spatial, aes(color=species)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    theme_minimal()
  
#eulipotyphla
  microsats$eulipotyphla <- ifelse(microsats$order == "eulipotyphla", 1, 0)
  microsats_eulipotyphla <- microsats[microsats$eulipotyphla == 1,]
  eulipotyphla_spatial <- st_as_sf(microsats_eulipotyphla, coords = c('lon', 'lat'), crs = 4326)
  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = eulipotyphla_spatial, aes(color=species)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    theme_minimal()
  
#lagomorpha
  microsats$lagomorpha <- ifelse(microsats$order == "lagomorpha", 1, 0)
  microsats_lagomorpha <- microsats[microsats$lagomorpha == 1,]
  lagomorpha_spatial <- st_as_sf(microsats_lagomorpha, coords = c('lon', 'lat'), crs = 4326)
  ggplot() + 
    geom_sf(data = ca_eu) +
    geom_sf(data = lagomorpha_spatial, aes(color=species)) +
    coord_sf(xlim = xlim, ylim = ylim) +
    theme_minimal()
  
  