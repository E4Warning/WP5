## Author: Murat Biricik, Max Planck Institute for Animal Behaviour, Germany
## Project: E4Warning, WP5, Milestone M5.1
## Date: 2024

# Habitat Suitability Modeling for Birds in a Certain Study Site & Period

#-------------------------------------------------------------------------------

# Study site selection

# set study site ("CH" / "DE" / "ES" / "GR")
site_selected <- "..."

# (or...)
if (FALSE) {
# choose one of the study sites
sites <- c("CH", "DE", "ES", "GR")
site_selected <- select.list(sites, title = "Choose study site")
rm(sites)
}

#-------------------------------------------------------------------------------
# (use this part only once per session!)

# START skip 1
{

# Program data

# load packages
{
library(biomod2)
library(dplyr)
library(ecospat)
library(ENMTools)
library(factoextra)
library(modEvA)
library(pROC)
library(raster)
library(SDMTools)
library(sf)
library(spThin)
library(terra)
library(USE) 
}

# set the working directory
setwd("C:/...")

# source
source("funcs.RDS")

#---------------------------------

# Initial settings

# edit path & name of the study site (SS) shapefile
SS_file <- paste0("./Data/StudySites/...", site_selected, ".gpkg")

# get the shape file
SS <- read_sf(SS_file)

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

}
# END skip 1

#-------------------------------------------------------------------------------
# (skip this part after first use as long as the selected study site remains constant!..)

# START skip 2
{

# Environmental data - Eu

# load rasters of the candidate variables clipped by the study area

{
  # bioclimate (WorldClim data) 
  Bio02         <- raster("./Data/rs_Eu/Bio02.tif")
  Bio03         <- raster("./Data/rs_Eu/Bio03.tif")
  Bio08         <- raster("./Data/rs_Eu/Bio08.tif")
  Bio09         <- raster("./Data/rs_Eu/Bio09.tif")
  Bio13         <- raster("./Data/rs_Eu/Bio13.tif")
  Bio15         <- raster("./Data/rs_Eu/Bio15.tif")

  # average day length
  DayLength     <- raster("./Data/rs_Eu/DayLength.tif")

  # distance to the nearest water body (as Jul 2022)
  DistWat       <- raster("./Data/rs_Eu/DistWat.tif")
  
  # elevation above sea level (WorldClim DEM data)
  Elev          <- raster("./Data/rs_Eu/Elev.tif")

  # average NDVI (2013-2022)
  NDVI          <- raster("./Data/rs_Eu/NDVI.tif")
  
  # night light (NASA data)
  NLight        <- raster("./Data/rs_Eu/NLight.tif")
  
  # average monthly precipitation (historical WorldClim data)
  PrecAvrg      <- raster("./Data/rs_Eu/PrecAvrg.tif")

  # soil categories
  Soil          <- raster("./Data/rs_Eu/Soil.tif")
  
  # average monthly temperature (historical WorldClim data)
  TempAvrg      <- raster("./Data/rs_Eu/TempAvrg.tif")
  
  # land cover (WorldCover data)
    # built-up
  WCovBuilt     <- raster("./Data/rs_Eu/WCovBuilt.tif")
    # crop land
  WCovCrop      <- raster("./Data/rs_Eu/WCovCrop.tif")
    # grass land
  WCovGrass     <- raster("./Data/rs_Eu/WCovGrass.tif")
    # shrub land
  WCovShrub     <- raster("./Data/rs_Eu/WCovShrub.tif")
    # tree cover
  WCovTree      <- raster("./Data/rs_Eu/WCovTree.tif")
  
  # human population (WorldPop data)
  WorldPop      <- raster("./Data/rs_Eu/WorldPop.tif")
}

# create a raster stack
rs_Eu <- stack(Bio02, Bio03, Bio08, Bio09, Bio13, Bio15, DayLength, DistWat,
               Elev, NDVI, NLight, PrecAvrg, Soil, TempAvrg,
               WCovBuilt, WCovCrop, WCovGrass, WCovShrub, WCovTree, WorldPop)

# set coordinate reference system
crs(rs_Eu) <- "+init=epsg:4326 +WGS84"

#----------------------------------

# Environmental data - Study Site (SS)

# crop the raster stack
rs_SS <- mask(rs_Eu, SS)
rs_SS <- crop(rs_SS, SS)

# visualization
# plot(rs_SS)

# get rid of unnecessaries
rm(rs_Eu)

#----------------------------------

# Bird data

# edit path & file name of the SS-observations
SS_file <- paste0("./Data/Obs/", site_selected, "...RDS")

# get all observations in the study site
SS_allObs <- readRDS(SS_file)

# list of species for SS
splist_SS <- unique(SS_allObs$species)

# get rid of unnecessaries
rm(SS_file)

#----------------------------------

}
# END skip 2

#-------------------------------------------------------------------------------

# select the working species
# sp <- select.list(splist_SS, title = "Choose species")

for (i in 1:length(splist_SS)) {
  sp <- splist_SS[i]

# START main block
{

# select species observation data from all observations in the SS
sp_Occ_SS <- subset(SS_allObs, species == sp)

# select columns to be used
sp_Occ_SS <- dplyr::select(sp_Occ_SS, c(decimalLongitude, decimalLatitude))

# remove NAs, if any
sp_Occ_SS <- na.omit(sp_Occ_SS)

# remove duplicates
sp_Occ_SS <- sp_Occ_SS %>% distinct()

# first check of the number of observations (should be 12 at least!..)
if(nrow(sp_Occ_SS) < 12) { print(sp); print("low number of occurrence!.."); next }

# add a species column to define occurrences
sp_Occ_SS <- data.frame(species = 1, sp_Occ_SS)

#----------------------------------

# Spatial thinning

# thinning
thin(sp_Occ_SS,
     long.col       = "decimalLongitude",
     lat.col        = "decimalLatitude",
     spec.col       = "species",
     thin.par       = 2,                             # set minimum distance (km)
     reps           = 1,
     out.dir        = getwd(),
     write.log.file = F,
     max.files      = 1)

# get thinned data & delete the file created by 'spThin'
sp_Occ_SS <- read.csv("./thinned_data_thin1.csv")
file.remove("./thinned_data_thin1.csv")

#----------------------------------

# re-check the number of observations (should be 10 at least!..)
if(nrow(sp_Occ_SS) < 10) { print(sp); print("low number of occurrence!.."); next }

# rename columns
colnames(sp_Occ_SS) <- c("species", "decimalLongitude", "decimalLatitude")

#-------------------------------------------------------------------------------

# Uniform Sampling of the Pseudo-absences (PAs)

#----------------------------------

# make occurrences spatial
sp_Occ_SS_st <- st_as_sf(sp_Occ_SS, 
                         coords = c("decimalLongitude", "decimalLatitude"))

# Sampling

# set the kernel density boundary threshold
thres <- 0.75

# reduce the threshold and retry if error occurs
# (be patient!..)

set.seed(12345)

while (!exists("PA_SS") & thres > 0.5) {
  try(
    PA_SS <- paSampling(env.rast  = rs_SS,
                        pres      = sp_Occ_SS_st,
                        thres     = thres,
                        grid.res  = OptRes,
                        prev      = 0.1,
                        plot_proc = F),
  )
  thres <- thres - 0.02
}

# check task fulfillment
if(!exists("PA_SS")) { print(sp); print("low number of occurrence!.."); next }

# get rid of unnecessaries
rm(sp_Occ_SS_st, thres)

#----------------------------------

# get coordinates of PAs for SS
sp_PA_Coord_SS <- data.frame("decimalLongitude" = PA_SS$x,
                             "decimalLatitude"  = PA_SS$y)

# add a species column to define pseudo-absences (PAs)
sp_PA_SS <- data.frame(species = 0, sp_PA_Coord_SS)

# get rid of unnecessaries
rm(PA_SS)

#----------------------------------

# Removing PAs that are too close (<0.03 degrees) to the occurrences

# rename columns of coordinates
names(sp_Occ_SS)       <- c("X", "Y")
names(sp_PA_Coord_SS)  <- c("X", "Y")

# thinning PAs
sp_PA_Coord_SS <- RemoveClosePAs(sp_Occ_SS, sp_PA_Coord_SS, 0.03)

# add a species column
sp_PA_SS <- data.frame(species = 0, sp_PA_Coord_SS)

# rename columns back
names(sp_Occ_SS) <- c("species", "decimalLongitude", "decimalLatitude")
names(sp_PA_SS)  <- c("species", "decimalLongitude", "decimalLatitude")

#-------------------------------------------------------------------------------

# Modeling

# combine Occs & PAs
sp_SS <- rbind(sp_Occ_SS, sp_PA_SS)

# coordinates
sp_SS_Coord <- subset(sp_SS, select = c(decimalLongitude, decimalLatitude))

# sample environmental variables
SS_data <- data.frame(extract(rs_SS, sp_SS_Coord))

# add species & coordinate columns
SS_data <- data.frame(sp_SS, SS_data)

# remove NAs, if any
SS_data <- na.omit(SS_data)

# check the number of observations (should be 10 at least!..)
if(nrow(SS_data[SS_data$species == 1, ]) < 10) { print(sp); print("low number of occurrence!.."); next }

# get rid of unnecessaries
rm(sp_PA_Coord_SS, sp_PA_SS, sp_SS, sp_SS_Coord)

#----------------------------------

# Removing collinear variables

SS_data <- CheckCollinearity(SS_data)

#----------------------------------

# Creating PCA rasters

# select environmental variables
SS_data_pca <- subset(SS_data, select = -c(species, decimalLongitude, decimalLatitude))

# check NAs
# apply(is.na(SS_data_pca), 2, which)

# check & remove variables with no variance (relevant for small areas such as SS!)
j <- 1
while (j <= ncol(SS_data_pca)) {
  if(var(SS_data_pca[, j]) == 0) { SS_data_pca[, j] <- NULL }
  j = j + 1
}

# list of final variables
EnvVars <- colnames(SS_data_pca)

# create a final raster stack of PCA components
rs_SS_pca <- subset(rs_SS, EnvVars)

# PCA
SS_pca <- prcomp(SS_data_pca, scale = T)

# PCA results
# summary(SS_pca)

# eigenvalues
eigval <- get_eigenvalue(SS_pca)

# define the number of PCA rasters to be inputted in the modeling
pca_no <- 0

for (j in 1:nrow(eigval)) {
  if (eigval[j, "eigenvalue"] > 1) {pca_no <- pca_no + 1}
}

# create PCA rasters in defined numbers
rs_SS_pca <- raster.pca(rs_SS_pca, pca_no)

# create a stack of pca rasters
rs_SS_pca <- stack(rs_SS_pca["rasters"])

# visualization (an example)
# plot(rs_SS_pca@layers[[3]])

# get rid of unnecessaries
rm(SS_pca, SS_data_pca, eigval)

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

# arrange species name
# replace " " (blank) with "_" (dash)
j <- regexpr(" ", sp)
regmatches(sp, j) <- "_"

# sample training & testing data
k <- 1
while (k < 5) {
  # randomly select 70% of the rows as training data
  j <- sample(seq_len(nrow(SS_alldata)), size = round(0.7 * nrow(SS_alldata)))

  # subset the training and testing data
  data_train <- SS_alldata[j, ]
  data_test  <- SS_alldata[-j, ]

  # check the number of observations in the testing data (should be 2 at least!..), if not, re-sample data
  ifelse((nrow(data_test[data_test$species == 1, ]) < 2), k <- k + 1, k <- 5)
}

# check the number of observations in the testing data (should be 2 at least!..)
if(nrow(data_test[data_test$species == 1, ]) < 2) { print(sp); print("low number of occurrence in the testing data!.."); next }

# get corresponding Occ & PA data
spec_train <- as.numeric(data_train[, "species"])
spec_test  <- as.numeric(data_test[, "species"])

# get corresponding coordinates
coord_train <- data_train[, c("decimalLongitude", "decimalLatitude")]
coord_test  <- data_test[, c("decimalLongitude", "decimalLatitude")]

# get rid if unnecessaries
rm(data_train, data_test, j, k)

#----------------------------------

# format data for 'biomod2'
(data_bm <- BIOMOD_FormatingData(resp.var      = spec_train,
                                 resp.xy       = coord_train,
                                 expl.var      = rs_SS_pca,
                                 resp.name     = sp,
                                 eval.resp.var = spec_test,
                                 eval.resp.xy  = coord_test,
                                 eval.expl.var = rs_SS_pca,
                                 na.rm         = T))

# set modeling options
options_bm <- BIOMOD_ModelingOptions(GAM = list(k                 =-1,
                                                interaction.level = 1),
                                     GBM = list(distribution      = "bernoulli",
                                                n.trees           = 100,
                                                interaction.depth = 2,
                                                shrinkage         = 0.1,
                                                bag.fraction      = 0.5),
                                     GLM = list(interaction.level = 1),
                                     MAXENT = NULL,
                                     RF  = list(mtry              = 3,
                                                ntree             = 100,
                                                nodesize          = 10))

# single models
model_bm <- BIOMOD_Modeling(bm.format        = data_bm,
                            modeling.id      = "AllModels",
                            models           = c("GAM",
                                                 "GBM",
                                                 "GLM", 
                                                 "MAXENT",
                                                 "RF"),
                            bm.options       = options_bm,
                            metric.eval      = c("TSS", "ROC", "KAPPA"),
                            nb.rep           = 10,
                            data.split.table = NULL,
                            data.split.perc  = 70,
                            var.import       = 5,
                            do.full.models   = T,
                            seed.val         = 12345)

# get rid if unnecessaries
rm(coord_test, coord_train, data_bm, spec_test)

#-------------------------------------------------------------------------------
# Evaluation - single models skipped!
#-------------------------------------------------------------------------------

# Ensemble model

thres <- 0.7
try(
  while (thres > 0.6) {
    model_ens_bm <- BIOMOD_EnsembleModeling(bm.mod               = model_bm,
                                            models.chosen        = "all",
                                            em.by                = "all",
                                            metric.select        = "TSS",
                                            metric.select.thresh = thres,
                                            metric.eval          = c("TSS", "ROC", "KAPPA"),
                                            var.import           = 5,
                                            em.algo              = "EMwmean",
                                            seed.val             = 12345)
    ifelse((!exists("model_ens_bm") | model_ens_bm@em.computed == "none"),
            thres <- thres - 0.02,
            thres <- 0)
  }
)

# check task fulfillment
if(!exists("model_ens_bm") | model_ens_bm@em.computed == "none") { print(sp); print("ensemble modeling failed!.."); next }

#-------------------------------------------------------------------------------

# Evaluation - ensemble model

# START skip Evaluation - ensemble model
if (FALSE) {

# Calculating ROC

# get predictions
df <- get_predictions(model_ens_bm)

# methods used
methods <- unique(df$algo)

# Calculating & saving

for (i in 1:length(methods)) {
  # select method
  method <- methods[i]
  
  # get & arrange prediction data
  df1 <- subset(df, algo == method)
  df1 <- data.frame(points = df1$points, pred = df1$pred)
  df1 <- as.vector(df1)
  prd <- split(df1$pred, df1$points)
  preds <- 1:length(prd)
  
  # calculate average predictions
  for (j in preds) {
    preds[j] <- mean(unlist(prd[j]))
  }

  # get observations & calculate ROC
  r <- roc(spec_train ~ preds)
  
  # create random predictions for comparison
  set.seed(12345)
  preds_random <- runif(length(prd), min = 0, max = 1)

  # calculate ROC for random predictions
  r_random <- roc(spec_train ~ preds_random)

  # compare differences in area under ROC
  (r_test <- roc.test(r, r_random, method = "bootstrap", paired = T))
  
  # create line
  line <- paste(method, as.numeric(r_test$estimate[1]), r_test$p.value, sep = "\t")

  # add line to the results file 
  write(line, "./E4W_AlgoEval.txt", append = T)
}

# visualization
# plot(r)
# lines(r_random, col = "red")

# get rid if unnecessaries
rm(df, df1, i, j, line, method, methods, preds, preds_random, prd, r, r_random, r_test)

}
# END skip Evaluation - ensemble model

#-------------------------------------------------------------------------------
}
# END main block

#-------------------------------------------------------------------------------

# Prediction 2022

# start with the first month
mon <- 1
while(mon <= 12) {

# Creating new data

# creating new environmental raster stack for prediction

{
# edit new environmental data file names...

# day length
DayLengthFile <- paste0("F:/GISData/E4W/DayLength2022_Eu_month/DayLength2022_month_", mon, ".tif")

# NDVI
NDVIFile <- paste0("F:/GISData/E4W/NDVI2022_Eu_month/NDVI2022_month_", mon, ".tif")

# average precipitation
PrecipitationFile <- paste0("F:/GISData/E4W/Precipitation_average_Eu_month/Precipitation_average_", mon, ".tif")

# average temperature
TemperatureFile <- paste0("F:/GISData/E4W/Temperature_average_Eu_month/Temperature_average_", mon, ".tif")

# get new (dynamic) environmental data
DayLength     <- raster(DayLengthFile)
NDVI          <- raster(NDVIFile)
PrecAvrg      <- raster(PrecipitationFile)
TempAvrg      <- raster(TemperatureFile)
} 

# create a new raster stack
rs_Eu_new <- stack(mget(EnvVars))

# set coordinate reference system
crs(rs_Eu_new) <- "+init=epsg:4326 +WGS84"

# get rid of unnecessaries
rm(DayLengthFile, NDVIFile, PrecipitationFile, TemperatureFile)

#----------------------------------

# new environmental data - Study Site

# crop the raster stack
rs_SS_new <- mask(rs_Eu_new, SS)
rs_SS_new <- crop(rs_SS_new, SS)

#----------------------------------

# create new pca rasters in the same numbers already defined
rs_SS_pca_new <- raster.pca(rs_SS_new, pca_no)

# create a new stack of pca components
rs_SS_pca_new <- stack(rs_SS_pca_new["rasters"])

# get rid of unnecessaries
rm(rs_Eu_new, rs_SS_new)

#-------------------------------------------------------------------------------

# Prediction

# ensemble model projection
(proj_ens_bm <- BIOMOD_EnsembleForecasting(bm.em         = model_ens_bm,
                                           new.env       = rs_SS_pca_new,
                                           proj.name     = "Ensemble_Model_Prediction",
                                           models.chosen = "all",
                                           output.format = ".tif"))

#-------------------------------------------------------------------------------

# Prediction map

# name & location of ensemble map
em <- raster(proj_ens_bm@proj.out@link)

# visualization
# plot(em)

#-------------------------------------------------------------------------------

# Suitability map

# Selecting occurrence threshold 

# get observations
obs <- SS_alldata$species

# predictions
obs_coord <- SS_alldata[, c("decimalLongitude", "decimalLatitude")]
p <- extract(em, obs_coord)

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

# get raster values
df <- data.frame(values(em))

# duplicate column
df$replace <- as.numeric(df[, 1])

# rename columns
names(df) <- c("original", "replace")

# drop out values below threshold
df$replace <- ifelse(df$replace < OptTh, df$replace * 0, df$replace)

# keep unique values
df <- df[!duplicated(df$original, ), ]

# correct probability values
df$replace <- df$replace / 1000

# substitute raster values & save
# (replacing existing sp_file!..)

# edit species file name & its path
sp_file <- paste0("C:/.../Output/", site_selected, "_2022/", sp, "_", mon, ".tif")

subs(em, df,
     by         = "original",
     which      = "replace",
     subsWithNA = F,
     filename   = sp_file,
     overwrite  = T)

# visualization
plot(raster(sp_file), main = c(sp, mon))

# get rid of unnecessaries
rm(df, OptTh, sp_file)

mon <- mon + 1

}

rm(model_ens_bm, pca_no)

}

################################################################################
