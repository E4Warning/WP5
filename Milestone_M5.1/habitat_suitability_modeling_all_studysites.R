## Author: Murat Biricik, Max Planck Institute for Animal Behaviour, Germany
## Project: E4Warning, WP5, Milestone M5.1
## Date: 2024

# Modeling for the Study Sites
#-------------------------------------------------------------------------------

# set study site

## (introductory way; assign site to bypass selection steps below!..)
site_selected <- "ES"

# skip BLOCK 1 after first use in the session as long as study site not changed!

# BLOCK 1

if (FALSE) {

# Elementary data

# load packages
{
library(corrplot)
library(dataCompareR)
library(dplyr)
library(ENMTools)
library(factoextra)
library(GGally)
library(raster)
library(rgdal)
library(sdm)
library(SDMTools)
library(sf)
library(spThin)
library(terra)
library(usdm)
library(USE)
}

# working directory
setwd("C:/Users/mbiricik/Documents/E4W")

# source
source("funcs.RDS")

# map of Europe
Eu <- readOGR("./Data/Europe/Europe.shp")

#---------------------------------

# Initial settings

# set study site

## (introductory way; assign site to bypass selection steps below!..)
# site_selected <- "ES"

if (FALSE) {
# choose one of the study sites
sites <- c("CH", 
           "DE",
           "ES",
           "GR")
site_selected <- select.list(sites, title = "Choose study site")

# get rid of unnecessaries
rm(sites)

}

# edit SS map file name
SS_file <- paste0("./Data/StudySites/StudySiteBuffer100km_", site_selected, ".gpkg")

# get map
SS <- readOGR(SS_file)

# visualization
# plot(SS)
# plot(Eu, add = T)

# set optimal resolution of the grid for Uniform Sampling (if already calculated)
## OptRes <- 3  # CH
## OptRes <- 3  # DE
## OptRes <- 4  # ES
## OptRes <- 6  # GR
## OptRes <- NA # unknown

ifelse(site_selected == "CH" | site_selected == "DE", OptRes <- 3,
       ifelse(site_selected == "ES", OptRes <- 4,
              ifelse(site_selected == "GR", OptRes <- 6, OptRes <- NA)))

# get rid of unnecessaries
rm(SS_file)

#-------------------------------------------------------------------------------

# Environmental data - Eu

# load rasters of the candidate variables clipped by the study area
# (collinear variables excluded by #)
# (possibly collinear variables in SS excluded by ##)


{
  # bioclimate (WorldClim data) 
#  Bio01         <- raster("./Data/rs_Eu/Bio01.tif")
  Bio02         <- raster("./Data/rs_Eu/Bio02.tif")
  Bio03         <- raster("./Data/rs_Eu/Bio03.tif")
#  Bio04         <- raster("./Data/rs_Eu/Bio04.tif")
#  Bio05         <- raster("./Data/rs_Eu/Bio05.tif")
#  Bio06         <- raster("./Data/rs_Eu/Bio06.tif")
#  Bio07         <- raster("./Data/rs_Eu/Bio07.tif")
  Bio08         <- raster("./Data/rs_Eu/Bio08.tif")
  Bio09         <- raster("./Data/rs_Eu/Bio09.tif")
#  Bio10         <- raster("./Data/rs_Eu/Bio10.tif")
#  Bio11         <- raster("./Data/rs_Eu/Bio11.tif")
#  Bio12         <- raster("./Data/rs_Eu/Bio12.tif")
  Bio13         <- raster("./Data/rs_Eu/Bio13.tif")
#  Bio14         <- raster("./Data/rs_Eu/Bio14.tif")
  Bio15         <- raster("./Data/rs_Eu/Bio15.tif")
#  Bio16         <- raster("./Data/rs_Eu/Bio16.tif")
#  Bio17         <- raster("./Data/rs_Eu/Bio17.tif")
#  Bio18         <- raster("./Data/rs_Eu/Bio18.tif")
#  Bio19         <- raster("./Data/rs_Eu/Bio19.tif")

  # day length on 1st Jul 2022 (calendar day 181)
  DayLength_181 <- raster("./Data/rs_Eu/DayLength_181.tif")
  
  # distance to the nearest water body as Jul 2022
  DistWat_Jul22 <- raster("./Data/rs_Eu/DistWat_Jul22.tif")
  
  # elevation above sea level  (WorldClim data)
  Elev          <- raster("./Data/rs_Eu/Elev.tif")

  # night light (NASA data)
  NLight        <- raster("./Data/rs_Eu/NLight.tif")
  
  # average monthly precipitation (historical WorldClim data)
#  PrecAvrg01    <- raster("./Data/rs_Eu/PrecAvrg01.tif")
#  PrecAvrg02    <- raster("./Data/rs_Eu/PrecAvrg02.tif")
#  PrecAvrg03    <- raster("./Data/rs_Eu/PrecAvrg03.tif")
#  PrecAvrg04    <- raster("./Data/rs_Eu/PrecAvrg04.tif")
  PrecAvrg05    <- raster("./Data/rs_Eu/PrecAvrg05.tif")
#  PrecAvrg06    <- raster("./Data/rs_Eu/PrecAvrg06.tif")
  PrecAvrg07    <- raster("./Data/rs_Eu/PrecAvrg07.tif")
#  PrecAvrg08    <- raster("./Data/rs_Eu/PrecAvrg08.tif")
  PrecAvrg09    <- raster("./Data/rs_Eu/PrecAvrg09.tif")
#  PrecAvrg10    <- raster("./Data/rs_Eu/PrecAvrg10.tif")
#  PrecAvrg11    <- raster("./Data/rs_Eu/PrecAvrg11.tif")
#  PrecAvrg12    <- raster("./Data/rs_Eu/PrecAvrg12.tif")
  
  # soil categories
  Soil          <- raster("./Data/rs_Eu/Soil.tif")
  
  # average monthly temperature (historical WorldClim data)
#  TempAvrg01    <- raster("./Data/rs_Eu/TempAvrg01.tif")
#  TempAvrg02    <- raster("./Data/rs_Eu/TempAvrg02.tif")
#  TempAvrg03    <- raster("./Data/rs_Eu/TempAvrg03.tif")
#  TempAvrg04    <- raster("./Data/rs_Eu/TempAvrg04.tif")
#  TempAvrg05    <- raster("./Data/rs_Eu/TempAvrg05.tif")
#  TempAvrg06    <- raster("./Data/rs_Eu/TempAvrg06.tif")
#  TempAvrg07    <- raster("./Data/rs_Eu/TempAvrg07.tif")
#  TempAvrg08    <- raster("./Data/rs_Eu/TempAvrg08.tif")
#  TempAvrg09    <- raster("./Data/rs_Eu/TempAvrg09.tif")
#  TempAvrg10    <- raster("./Data/rs_Eu/TempAvrg10.tif")
#  TempAvrg11    <- raster("./Data/rs_Eu/TempAvrg11.tif")
#  TempAvrg12    <- raster("./Data/rs_Eu/TempAvrg12.tif")
  
  # land cover (WorldCover data)
    # built-up
  WCovBuilt     <- raster("./Data/rs_Eu/WCovBuilt.tif")
    # crop land
  WCovCrop      <- raster("./Data/rs_Eu/WCovCrop.tif")
    # grass land
  WCovGrass     <- raster("./Data/rs_Eu/WCovGrass.tif")
    # shrub land
  WCovShrub     <- raster("./Data/rs_Eu/WCovShrub.tif")
    # snow and ice
##  WCovSnow      <- raster("./Data/rs_Eu/WCovSnow.tif")
    # bare / sparse vegetation
##  WCovSparBare  <- raster("./Data/rs_Eu/WCovSparBare.tif")
    # tree cover
  WCovTree      <- raster("./Data/rs_Eu/WCovTree.tif")
  
  # human population (WorldPop data)
  WorldPop      <- raster("./Data/rs_Eu/WorldPop.tif")

} 

# create a raster stack

rs_Eu <- stack(
  # Bio01, 
  Bio02, Bio03, 
  # Bio04, Bio05, Bio06, Bio07, 
  Bio08, Bio09,
  # Bio10,  Bio11, Bio12, 
  Bio13, 
  # Bio14, 
  Bio15, 
  # Bio16, Bio17, Bio18, Bio19,
  DayLength_181,
  DistWat_Jul22,
  Elev,
  NLight,
  # PrecAvrg01, PrecAvrg02, PrecAvrg03, PrecAvrg04, 
  PrecAvrg05, 
  # PrecAvrg06,
  PrecAvrg07, 
  # PrecAvrg08, 
  PrecAvrg09,
  # PrecAvrg10, PrecAvrg11, PrecAvrg12,
  Soil,
  # TempAvrg01, TempAvrg02, TempAvrg03, TempAvrg04, TempAvrg05, TempAvrg06,
  # TempAvrg07, TempAvrg08, TempAvrg09, TempAvrg10, TempAvrg11, TempAvrg12,
  WCovBuilt, WCovCrop, WCovGrass, WCovShrub, 
  ## WCovSnow, WCovSparBare,
  WCovTree,
  WorldPop
)

# set coordinate reference system
crs(rs_Eu) <- "+init=epsg:4326 +WGS84"

# get rid of single rasters
rm(
  Bio01, Bio02, Bio03, Bio04, Bio05, Bio06, Bio07, Bio08, Bio09, Bio10,
  Bio11, Bio12, Bio13, Bio14, Bio15, Bio16, Bio17, Bio18, Bio19,
  DayLength_181,
  DistWat_Jul22,
  Elev,
  NLight,
  PrecAvrg01, PrecAvrg02, PrecAvrg03, PrecAvrg04, PrecAvrg05, PrecAvrg06,
  PrecAvrg07, PrecAvrg08, PrecAvrg09, PrecAvrg10, PrecAvrg11, PrecAvrg12,
  Soil,
  TempAvrg01, TempAvrg02, TempAvrg03, TempAvrg04, TempAvrg05, TempAvrg06,
  TempAvrg07, TempAvrg08, TempAvrg09, TempAvrg10, TempAvrg11, TempAvrg12,
  WCovBuilt, WCovCrop, WCovGrass, WCovShrub, WCovSnow, WCovSparBare, WCovTree,
  WorldPop
)

#----------------------------------

# Environmental data - Study Site

# crop the raster stack
rs_SS <- mask(rs_Eu, SS)
rs_SS <- crop(rs_SS, SS)

# visualization
# plot(rs_SS)

# get rid of unnecessaries
rm(rs_Eu)

#-------------------------------------------------------------------------------

# NDVI data source

## creating the full list of available dates of the NDVI files
## (skip this block after single use!..)
# NDVIfiles <- as.data.frame(list.files("F:/GISData/E4W/NDVI_Eu_resampled"))
# names(NDVIfiles) <- "NDVIfiles"
# NDVIfiles$NDVIfiles <- as.numeric(substring(NDVIfiles$NDVIfiles, 4, 11))
# saveRDS(NDVIfiles, "NDVIfiles.RDS")

# get the list of dates of available NDVI files
NDVIfiles <- readRDS("NDVIfiles.RDS")
NDVIfiles <- as.numeric(unlist(NDVIfiles))

# format dates of NDVI files in YYYY-MM-DD
NDVIfiles_formatted <- as.Date(paste0(substring(NDVIfiles, 1, 4), "-",
                                      substring(NDVIfiles, 5, 6), "-",
                                      substring(NDVIfiles, 7, 8)))

# convert in Julian dates
NDVIfiles_julian <- julian(NDVIfiles_formatted)

# get rid of unnecessaries
rm(NDVIfiles_formatted)

#-----------------

# Adding NDVI raster to the raster stack

# get NDVI raster (which were created at later steps)
# ("NDVI.tif" = "NDVI_Eu_yearly_average.tif"!..)
NDVI <- raster("./Data/rs_Eu/NDVI.tif")

# clipping
NDVI <- mask(NDVI, SS)
NDVI <- crop(NDVI, SS)

# add NDVI data to the stack
rs_SS_w_NDVI <- addLayer(rs_SS, NDVI)

# get rid of unnecessaries
rm(NDVI)

#-------------------------------------------------------------------------------

# Bird data

# edit SS-observations file name
SS_file <- paste0("./Data/Obs/", site_selected, "_allspecies_light.CSV")

# get all observations in SS
SS_allObs <- read.csv(SS_file)

# get list of bird species which observed in any of SS at least once!
# sp_SS <- read.csv("./Data/Obs/SSBirds.csv")

# get list of bird species which observed in SS at least once!

# edit SS-observations file name
SS_list_file <- paste0("./Data/Obs/", site_selected, "_Birds.CSV")
splist_SS <- read.csv(SS_list_file)
splist_SS <- as.vector(splist_SS$species)

# get rid of unnecessaries
rm(SS_file, SS_list_file)

# end of BLOCK 1
}


# select working species
# (list sorted from the most frequently recorded species to the least recorded species!..)
sp <- select.list(splist_SS, title = "Choose species")


## (this is only for tracking of the selected species on Global Environment!..)
sp_list_no <- as.numeric(match(sp, splist_SS))


for (i in 11:15) {

sp <- splist_SS[i]

# BLOCK 2
{

################################################################################

# NDVI Data
#-----------

# Data arrangement

# select species observation data from all observations in SS
sp_Occ_SS <- subset(SS_allObs, species == sp)

# select columns to be used
sp_Occ_SS <- dplyr::select(sp_Occ_SS, c(decimalLongitude, decimalLatitude, day, month, year))

# remove NAs
sp_Occ_SS <- na.omit(sp_Occ_SS)

# remove duplicates
sp_Occ_SS <- sp_Occ_SS %>% distinct()

# format dates of observation
sp_Occ_SS$eventDate <- as.Date(paste(sp_Occ_SS$year, sp_Occ_SS$month, sp_Occ_SS$day, sep = "-"))

# add a column of Julian date
sp_Occ_SS$eventDate_julian <- julian(sp_Occ_SS$eventDate)

# select columns to be used
sp_Occ_SS <- dplyr::select(sp_Occ_SS, c(decimalLongitude, decimalLatitude, eventDate, eventDate_julian))

#---------------------------------

# Getting NDVI data for the most close date of observation

# detect the most close date of the NDVI file to the observation date
fileindex <- 1:nrow(sp_Occ_SS)
for (i in fileindex) {
  fileindex[i] <- which.min(abs(NDVIfiles_julian - sp_Occ_SS$eventDate_julian[i]))
}

# add a column of most close NDVI-dates
sp_Occ_SS$NDVIfile <- NDVIfiles[fileindex]

# define link to corresponding NDVI file
sp_Occ_SS$NDVIfile <- paste0("F:/GISData/E4W/NDVI_Eu_resampled/Eu_", sp_Occ_SS$NDVIfile, ".tif")

# sort NDVI files (to speed up)
sp_Occ_SS <- sp_Occ_SS[order(sp_Occ_SS$NDVIfile),]

# remove the column unneeded
sp_Occ_SS <- subset(sp_Occ_SS, select = -c(eventDate_julian))

# add a species column to define occurrence
sp_Occ_SS$species = 1

# relocate species column
sp_Occ_SS <- relocate(sp_Occ_SS, species)

# get rid of unnecessaries
rm(fileindex, i)

#-------------------------------------------------------------------------------

# Uniform Sampling of the Pseudo-absences (PAs)

#----------------------------------

# Getting optimal resolution of the grid for sampling the environmental space
## (skip after single use!)

if (is.na(OptRes)){

# Generating the PCA environmental space

df_PCA_SS_st <- UniformSampling_PCAenv(rs_SS)

# Calculation of optimal resolution
## (be patient!..)

optRes_SS <- optimRes(sdf      = df_PCA_SS_st,
                      grid.res = c(1:10),
                      perc.thr = 20,
                      showOpt  = T, 
                      cr       = 5)

OptRes <- as.numeric(optRes_SS$Opt_res)

# get rid of unnecessaries
rm(df_PCA_SS_st, optRes_SS)

}

#----------------------------------

# Uniform sampling of the pseudo-absences (PA) within the environmental space

# get occurrence coordinates
sp_Occ_Coord_SS <- dplyr::select(sp_Occ_SS, c(decimalLongitude, decimalLatitude))

# make occurrences SS spatial
sp_Occ_SS_st <- st_as_sf(sp_Occ_Coord_SS, 
                         coords = c("decimalLongitude", "decimalLatitude"))

#----------------------------------

# Sampling

# set threshold
thres <- 0.75

# reduce kernel density boundary threshold and retry if error occurs
## (be patient!..)

set.seed(12345)

while (!exists("PA_SS") & thres > 0) {
  try(
    PA_SS <- paSampling(env.rast  = rs_SS,        # raster stack of environmental variables
                        pres      = sp_Occ_SS_st, # presences
                        thres     = thres,        # boundary of the kernel density estimate
                        H         = NULL,         # kernel bandwidth (NULL for estimating automatically)
                        grid.res  = OptRes,       # as.numeric(optRes_SS$Opt_res),
                        # resolution of the sampling grid
                        # n.tr      = 5,            # no. PAs for the training dataset to sample in each cell of the sampling grid
                        prev      = 0.07,         # prevalence value to be specified instead of n.tr and n.ts
                        # sub.ts    = T,            # sample the validation PAs 
                        # n.ts      = 5,            # no. PAs for the testing dataset to sample in each cell of the sampling grid
                        plot_proc = F,            # plot sampling progress
                        verbose   = F),
    thres <- thres - 0.01
    )
}

# get rid of unnecessaries
rm(sp_Occ_SS_st, sp_Occ_Coord_SS, thres)

#-------------------------------------------------------------------------------

# Checking modeling achievement

# check whether PAs could be created
ifelse(exists("PA_SS"), (modeling = T), (modeling = F))

try(
    if (!modeling) {
       # create line
       line <- paste(site_selected, sp, modeling, NA, sep = "\t")

       # add line to the log file 
       write(line, "logfile.txt", append = T, sep = "\n")
  
       # halt execution
       stop("MODELING FAILED!..")
     }
    )

# get rid of unnecessaries
rm(modeling)

#-------------------------------------------------------------------------------

# Data arrangement

# get coordinates of PAs for SS
sp_PA_SS <- data.frame("decimalLongitude" = PA_SS$x, "decimalLatitude" = PA_SS$y)

# define 8-times more Date & NDVI files for PA than that for Occ
eventDate <- rep(sp_Occ_SS$eventDate, 8)
NDVIfile  <- rep(sp_Occ_SS$NDVIfile, 8)

# get random coordinates
set.seed(12345)
sp_PA_SS <- sp_PA_SS[(sample(nrow(sp_PA_SS), length(eventDate), replace = F)), ]

# add columns
sp_PA_SS <- data.frame(sp_PA_SS, eventDate, NDVIfile)

# add a species column to define pseudo-absences (PAs)
sp_PA_SS$species = 0

# relocate species column
sp_PA_SS <- relocate(sp_PA_SS, species)

# check PAs
# head(sp_PA_SS)

# check occurrences (Occs)
# head(sp_Occ_SS)

# combine Occs & PAs
sp_SS <- rbind(sp_Occ_SS, sp_PA_SS)

# sort NDVI files (to speed up)
sp_SS <- sp_SS[order(sp_SS$NDVIfile),]

# get rid of unnecessaries
rm(sp_Occ_SS, sp_PA_SS, PA_SS, eventDate, NDVIfile)

#----------------------------------

# getting NDVI data
## (be patient!..)
  
for (i in 1:nrow(sp_SS)) {
  r <- raster(sp_SS$NDVIfile[i])
  rs <- stack(r)
  crs(rs) <- "+init=epsg:4326 +WGS84"
  sp_SS$NDVI[i] <- extract(rs, sp_SS[i, c("decimalLongitude", "decimalLatitude")])
}

# remove NAs
sp_SS <- na.omit(sp_SS)

# remove the columns unneeded
sp_SS_NDVI <- subset(sp_SS, select = -c(eventDate, NDVIfile))

# convert NDVI data
sp_SS_NDVI$NDVI <- as.numeric(sp_SS_NDVI$NDVI)

# get rid of unnecessaries
rm(i, r, rs, sp_SS)

#----------------------------------

# Spatial thinning - NDVI data

# data arrangement - occurrences
sp_Occ_SS_NDVI <- subset(sp_SS_NDVI, species == 1)
sp_SS_thin <- subset(sp_Occ_SS_NDVI, select = c(decimalLongitude, decimalLatitude))
sp_SS_thin$species <- sp

# thinning occurrences
thin(sp_SS_thin, 
     long.col       = "decimalLongitude", 
     lat.col        = "decimalLatitude", 
     spec.col       = "species",
     thin.par       = 3,                             # set minimum distance (km)
     reps           = 1,
     out.dir        = getwd(),
     write.log.file = F,
     max.files      = 1)

# get thinned data and delete the file created by 'spThin'
sp_Occ_SS_NDVI_thinned <- read.csv("./thinned_data_thin1.csv")
file.remove("./thinned_data_thin1.csv")

# data arrangement - pseudo absences (PAs)
sp_PA_SS_NDVI <- subset(sp_SS_NDVI, species == 0)
sp_SS_thin <- subset(sp_PA_SS_NDVI, select = c(decimalLongitude, decimalLatitude))
sp_SS_thin$species <- sp

# thinning PAs
thin(sp_SS_thin, 
     long.col       = "decimalLongitude", 
     lat.col        = "decimalLatitude", 
     spec.col       = "species",
     thin.par       = 3,                             # set minimum distance (km)
     reps           = 1,
     out.dir        = getwd(),
     write.log.file = F,
     max.files      = 1)

# get thinned data and delete the file created by 'spThin'
sp_PA_SS_NDVI_thinned <- read.csv("./thinned_data_thin1.csv")
file.remove("./thinned_data_thin1.csv")

# re-assign species
sp_Occ_SS_NDVI_thinned$species <- 1
sp_PA_SS_NDVI_thinned$species  <- 0

# combine data
sp_SS_NDVI_thinned <- rbind(sp_Occ_SS_NDVI_thinned, sp_PA_SS_NDVI_thinned)

# get rid of unnecessaries
# rm(sp_SS_thin, sp_Occ_SS_NDVI, sp_PA_SS_NDVI)

#-------------------------------------------------------------------------------

# Checking modeling achievement

# check no. occurrences
ifelse(nrow(sp_Occ_SS_NDVI_thinned) > 11, (modeling = T), (modeling = F))

try(
    if (!modeling) {
       # create line
       line <- paste(site_selected, sp, modeling, NA, sep = "\t")

       # add line to the log file 
       write(line, "logfile.txt", append = T, sep = "\n")
  
       # halt execution
       stop("MODELING FAILED!..")
     }
    )

# get rid of unnecessaries
rm(modeling)

#-------------------------------------------------------------------------------




#----------------------------------

# Removing PAs that are too close (<0.04 degrees) to the occurrences

# data arrangement
sp_Occ_SS_NDVI <- subset(sp_SS_NDVI_thinned, species == 1)
sp_Occ_SS_NDVI <- subset(sp_Occ_SS_NDVI_thinned, 
                         select = c(decimalLongitude, decimalLatitude))
names(sp_Occ_SS_NDVI) <- c("X", "Y")

sp_PA_SS_NDVI <- subset(sp_SS_NDVI_thinned, species == 0)
sp_PA_SS_NDVI <- subset(sp_PA_SS_NDVI_thinned, 
                        select = c(decimalLongitude, decimalLatitude))
names(sp_PA_SS_NDVI) <- c("X", "Y")

# thinning PAs
sp_PA_SS_NDVI_thinned <- RemoveClosePAs(sp_Occ_SS_NDVI, sp_PA_SS_NDVI, 0.04)

# add species column
sp_PA_SS_NDVI_thinned$species <- 0

# relocate species column
sp_PA_SS_NDVI_thinned <- relocate(sp_PA_SS_NDVI_thinned, species)

# rename columns
names(sp_PA_SS_NDVI_thinned) <- c("species", "decimalLongitude", "decimalLatitude")

# get rid of unnecessaries
rm(sp_SS_NDVI_thinned, sp_Occ_SS_NDVI, sp_PA_SS_NDVI)

#----------------------------------

# Selection of NDVI data for thinned points (both Occs & PAs)

# check Occs
# head(sp_Occ_SS_NDVI_thinned)

# check PAs
# head(sp_PA_SS_NDVI_thinned)

# combine Occs & PAs
sp_SS_NDVI_thinned <- rbind(sp_Occ_SS_NDVI_thinned, sp_PA_SS_NDVI_thinned)

# remove species column
sp_SS_NDVI_thinned <- subset(sp_SS_NDVI_thinned, select = -species)

# check NDVI data
# head(sp_SS_NDVI)

# select NDVI data
sp_SS_NDVI <- merge(sp_SS_NDVI_thinned, sp_SS_NDVI, 
                    by = c("decimalLongitude", "decimalLatitude"))

# select columns to be used
sp_SS_NDVI <- subset(sp_SS_NDVI, 
                     select = c(species, decimalLongitude, decimalLatitude, NDVI))

# get rid of unnecessaries
rm(sp_Occ_SS_NDVI_thinned, sp_PA_SS_NDVI_thinned, sp_SS_NDVI_thinned)

################################################################################

# Modeling
#----------

# Data arrangement

# get coordinates
coords_SS <- subset(sp_SS_NDVI, select = c(decimalLongitude, decimalLatitude))

# sample environmental variables (other than NDVI)
SS_data <- data.frame(extract(rs_SS, coords_SS))

# add NDVI data
SS_data <- na.omit(data.frame(sp_SS_NDVI, SS_data))

# check no. Occs & PAs
table(SS_data$species)

# get rid of unnecessaries
rm(coords_SS, sp_SS_NDVI)

#----------------------------------

# Removing collinear variables

SS_data <- CheckCollinearity(SS_data)

#-------------------------------------------------------------------------------

# Creating PCA rasters

# select environmental variables
SS_data_pca <- subset(SS_data, select = -c(species, decimalLongitude, decimalLatitude))

# check NAs
# apply(is.na(SS_data_pca), 2, which)

# pca
SS_pca <- prcomp(SS_data_pca, scale = TRUE)

# summary(SS_pca)

eigval <- get_eigenvalue(SS_pca)

# scree plot
# fviz_eig(SS_pca)

# define the number of PCA rasters to be inputted in the modeling
pca_no <- 0

for (i in 1:nrow(eigval)) {
  if (eigval[i, "eigenvalue"] > 1) {pca_no <- pca_no + 1}
}

#-----------------

# create pca rasters in defined numbers
rs_SS_pca <- raster.pca(rs_SS_w_NDVI, pca_no)

# create a stack of pca rasters
rs_SS_pca <- stack(rs_SS_pca["rasters"])

# visualization (an example)
# plot(rs_SS_pca@layers[[3]])

# get rid of unnecessaries
rm(SS_pca, SS_data_pca, eigval, pca_no)

#----------------------------------

# Sampling

# sample predictor data
SS_alldata <- extract(rs_SS_pca, SS_data[, c("decimalLongitude", "decimalLatitude")])

# merge columns that indicate Occ and PA along with their coordinates
SS_alldata <- as.data.frame(cbind(species          = SS_data$species,
                                  decimalLongitude = SS_data$decimalLongitude,
                                  decimalLatitude  = SS_data$decimalLatitude,
                                  SS_alldata))

# check NAs
# apply(is.na(SS_alldata), 2, which)

# get rid of unnecessaries
rm(SS_data)

#----------------------------------

# Modeling (using package SDM)

# create "sdmData" required by the sdm package
d_SS <- sdmData(species ~ . + coords(decimalLongitude + decimalLatitude), train = SS_alldata)

# model building

modelSDM_SS <- sdm(species ~ . + coords(decimalLongitude + decimalLatitude), d_SS,
                   methods           = c("GAM",
                                         "GLM", 
                                         "GBM", 
                                         "Maxent",
                                         "RF"),
                   interaction.depth = 3,
                   replication       = c("boot", "subs"),
                   test.p            = 30,
                   n                 = 3,
                   parallelSetting   = list(ncore = 5, method = "parallel"))
modelSDM_SS

# get rid of unnecessaries
rm(d_SS)

#---------------------------------

# Evaluation

# (skip for speeding up!..)

if (F) {

  # average variable importance for all models
(vi_SS <- getVarImp(modelSDM_SS))

# plot variable importance (using R standard plot instead of ggplot)
par(mfrow = c(1, 1))
plot(vi_SS, main = "Model-SS", gg = F)

# response curve
rcurve(modelSDM_SS, main = "Model-SS")

}

#---------------------------------

# Ensemble model SS

# arrange species name
# replace " " (blank) with "_" (dash)
i <- regexpr(" ", sp)
regmatches(sp, i) <- "_"

# edit species file name & its path
sp_file <- paste0("./Output/", site_selected, "/", sp, ".tif")

ens_SS <- ensemble(modelSDM_SS,
                   id = c,
                   newdata  = rs_SS_pca,
                   filename = sp_file,
                   setting  = list(method = "weighted", stat = "TSS"))

plot(ens_SS, main = sp)
# plot(Eu, add = T)
# occs <- subset(SS_alldata, species == 1)
# points(occs[, c("decimalLongitude", "decimalLatitude")])

# get rid of unnecessaries
# rm(occs)
rm(modelSDM_SS, rs_SS_pca)

#---------------------------------

# Ensemble model evaluation

# get coordinates
coord <- SS_alldata[, c("decimalLongitude", "decimalLatitude")]

# sampling
ens_SS_eval <- extract(ens_SS, coord)

# evaluation
ev_ens_SS <- evaluates(SS_alldata[, "species"], ens_SS_eval)

# get AUC value
auc_ens <- ev_ens_SS@statistics$AUC

# visualization
# plot(ev_ens_SS, main = sp)

#-------------------------------------------------------------------------------

# Checking modeling achievement

# (set threshold for success!..)
ifelse(auc_ens > 0.8, (modeling = T), (modeling = F))

try(
    if (!modeling) {
       # create line
       line <- paste(site_selected, sp, modeling, NA, sep = "\t")

       # add line to the log file 
       write(line, "logfile.txt", append = T, sep = "\n")
  
       # halt execution
       stop("MODELING FAILED!..")
     }
    )

# recording final state of achievement when modeling succeeded

# create line
line <- paste(site_selected, sp, modeling, auc_ens, sep = "\t")

# add line to the log file 
write(line, "logfile.txt", append = T, sep = "\n")

# get rid of unnecessaries
rm(ens_SS_eval, modeling, auc_ens, line)

#-------------------------------------------------------------------------------

# Suitability map

# Selecting occurrence threshold 

# get observations
obs <- SS_alldata$species

# predictions
obs_coord <- SS_alldata[, c("decimalLongitude", "decimalLatitude")]
p <- extract(ens_SS, obs_coord)

# get predictions & remove NAs
df <- na.omit(data.frame(obs, p))

# sort data
df <- df[order(df$obs, df$p), ]

# calculate optimum threshold
OptTh <- optim.thresh(df$obs, df$p)

# select mean occurrence prediction as threshold metric
OptTh <- OptTh$mean.occurence.prediction

# get rid of unnecessaries
rm(obs, obs_coord, df, p)

#---------------------------------

# Creating suitability raster

df <- data.frame(values(ens_SS))

# duplicate column
df$replace <- as.numeric(df[, 1])

# rename columns
names(df) <- c("original", "replace")

# drop out values below threshold
df$replace <- ifelse(df$replace < OptTh, df$replace * 0, df$replace)

# keep unique values
df <- df[!duplicated(df$original, ), ]

# substitute raster values & save
# (replacing existing sp_file!..)

subs(ens_SS, df,
     by         = "original",
     which      = "replace", 
     subsWithNA = F,
     filename   = sp_file,
     overwrite  = T)

# visualization
plot(raster(sp_file), main = sp)

# get rid of unnecessaries
rm(df, SS_alldata, OptTh, sp_file)

# end of BLOCK 2
}

}

graphics.off()

#-------------------------------

# rm(list = ls())

################################################################################
################################################################################
