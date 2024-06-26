## Author: Murat Biricik, Max Planck Institute for Animal Behaviour, Germany
## Project: E4Warning, WP5, Milestone M5.1
## Date: 2024

# Species composition for specified months of 2022 in the Study Site

#-------------------------------------------------------------------------------

# Main data

# load packages
{
library(raster)
}

# set the working directory
setwd("C:/...")

#---------------------------------

# Initial settings

# select study site ("CH" / "DE" / "ES" / "GR")
site_selected <- "..."

# set link to the folder that contains all available raster files
# (format of files in the folder: "[Genus]_[species]_[MonthNr].tif")
# path <- "..."
path <- paste0("C:/.../", site_selected, "_2022")

# select month (1 - 12)
mon_selected <- ...

#-------------------------------------------------------------------------------

{

# Creating information summaries

# get list of available rasters
list_rast <- list.files(path)

# create a data frame with rasters of species & months
{
list_sp <- unlist(strsplit(list_rast, "_"))
df_rast <- data.frame(split(list_sp, c(1:3)))
df_rast$sp <- paste0(df_rast$X1, "_", df_rast$X2)
df_rast$mon <- as.numeric(gsub(".tif", "", df_rast$X3))
df_rast  <- subset(df_rast[, 4:5])
}

# table with unique species & binary availability for all months 2022
# table_2022 <- table(df_rast)

# save table
# write.table(table_2022, paste0(site_selected, "_MonthlySpeciesAvailability2022.txt"))

# add raster file names to the data frame
df_rast$rast <- list_rast

# available rasters for the selected month
df_mon <- subset(df_rast, mon == mon_selected)

# get rid of unnecessaries
rm(list_rast, list_sp)

#-------------------------------------------------------------------------------

# Creating raster stacks

# list of rasters
r <- 1:nrow(df_mon)

for (i in 1:nrow(df_mon)) {
  r[i] <- paste0(path, "/", df_mon$rast[i])
}

# stack
rs <- stack(r)

# visualization (an example)
# plot(rs[[1]])

#---------------------------------

# average probability of occurrences for the selected month 
# plot(mean(rs), main = paste(site_selected, ", ", mon_selected, "/ 2022 - Average suitability"))

# save raster
# writeRaster(mean(rs), paste0(site_selected, "_AvrgSuit_", mon_selected, ".tif"))

#---------------------------------

# Binary raster of suitability

# create a reclassification matrix
rcl_mtx <- matrix(c(0, 1, 1), ncol = 3)

# reclassify raster stack
rs_bin <- reclassify(rs, rcl_mtx)

# visualization
# plot(rs_bin)
# plot(rs_bin[[1]])

#-------------------------------

# total species for the selected month 
plot(sum(rs_bin), main = paste(site_selected, ", ", mon_selected, "/ 2022 - Number of bird species"))

# save raster
# writeRaster(sum(rs_bin), paste0(site_selected, "_SpecNmbr_", mon_selected, ".tif"))

# get rid of unnecessaries
rm(i, r, rcl_mtx)

}

#-------------------------------------------------------------------------------

# List of species at a specified location

# (if not already plotted...)
# plot(sum(rs_bin), main = paste(site_selected, ", ", mon_selected, "/ 2022 - Number of bird species"))

# IMPORTANT! Screen display scale should be 100% for a proper locating!
# (Windows Settings -> System -> Display -> Scale)

# select location
# (mouse click & 'Esc')

{
# catch coordinates
coord <- locator()

# list of rasters
list_rast <- extract(rs_bin, data.frame(coord))

# list of species
list_sp <- gsub("_", " ", df_mon$sp)

# list of species at the selected locality
(list_sp <- list_sp[which(list_rast == 1)])
}

################################################################################
