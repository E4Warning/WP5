
## Author: Anne K. Scharf, Max Planck Institute for Animal Behaviour, Germany
## Project: E4Warning, WP5, Milestone M5.1
## Date: 2024

## This script downloads all available avian species per study site. The criteria are:
# * scientific name: Aves
# * occurrence status: present
# * has geospatial issue:F
# * has coordinates:T
# * coordinate uncertainty in meters: 0-100
# * year: between 1.Jan.2013 and 31.Dec.2022

library(sf)
library(rgbif)

# Coordinates of study sites – these points surrounded by 100Km
# # DE: 47.730333, 8.998664
# # CH: 46.155869, 8.860095
# # GR: 38.150904, 24.022343
# # ES: 42.226102, 3.087583
# Projection: "EPSG:4326"

## get polygons of study sites in, will have to be transformed to WKT for gbif with st_as_text()
de <- st_geometry(st_read("./Study_sites/StudySiteBuffer100km_DE.gpkg"))
ch <- st_geometry(st_read("./Study_sites/StudySiteBuffer100km_CH.gpkg"))
es <- st_geometry(st_read("./Study_sites/StudySiteBuffer100km_ES.gpkg"))
gr <- st_geometry(st_read("./Study_sites/StudySiteBuffer100km_GR.gpkg"))

## set credentials for gbif
usethis::edit_r_environ()
# GBIF_USER="****"
# GBIF_PWD="****"
# GBIF_EMAIL="****@ab.mpg.de"

## get data for study site DE
st_as_text(de) # have to paste in gbif--something changes
occ_de <- occ_download(
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"), 
  pred_and(pred_gte("year", 2013),pred_lte("year", 2022)),
  pred_within("MULTIPOLYGON(((10.33211 47.7226,10.27365 48.00131,10.08856 48.25402,9.79363 48.45547,9.41766 48.58524,8.99864 48.63002,8.57962 48.58523,8.20366 48.45545,7.90874 48.254,7.72366 48.00128,7.66522 47.72257,7.73717 47.44533,7.93061 47.19645,8.22554 46.99975,8.59317 46.87382,8.99869 46.8305,9.4042 46.87383,9.77182 46.99976,10.06674 47.19648,10.26017 47.44536,10.33211 47.7226)))"), 
  pred("taxonKey", 212), # dowload all of class Aves
  pred_lte("coordinateUncertaintyInMeters","100"),
  format = "SIMPLE_CSV"
)

occ_download_wait(occ_de)
dwn <- occ_download_get(occ_de,path="/home/ascharf/Documents/Projects/E4Warning")
# dta <- occ_download_import(dwn)
dta_de <- occ_download_import(as.download("/home/ascharf/Documents/Projects/E4Warning/0005972-240314170635999.zip"))
length(unique(dta_de$species))
# 302

## get data for study site CH
st_as_text(ch) # have to paste in gbif--something changes
occ_ch <- occ_download(
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"), 
  pred_and(pred_gte("year", 2013),pred_lte("year", 2022)),
  pred_within("MULTIPOLYGON(((10.15507 46.15012,10.09724 46.4288,9.91661 46.68127,9.62967 46.8823,9.26443 47.01155,8.85776 47.0558,8.45133 47.01054,8.0868 46.88039,7.80093 46.67866,7.62163 46.42576,7.56523 46.14695,7.63544 45.86968,7.82369 45.62086,8.11057 45.42432,8.46809 45.29869,8.86232 45.25579,9.25635 45.29964,9.61328 45.42614,9.89921 45.62338,10.08624 45.87267,10.15507 46.15012)))"), 
  pred("taxonKey", 212), # dowload all of class Aves
  pred_lte("coordinateUncertaintyInMeters","100"),
  format = "SIMPLE_CSV"
)

occ_download_wait(occ_ch)
dwn <- occ_download_get(occ_ch,path="/home/ascharf/Documents/Projects/E4Warning")
# dta_ch <- occ_download_import(dwn)
dta_ch <- occ_download_import(as.download("/home/ascharf/Documents/Projects/E4Warning/0006029-240314170635999.zip"))
length(unique(dta_ch$species))
# 297


## get data for study site ES
st_as_text(es) # have to paste in gbif--something changes
occ_es <- occ_download(
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"), 
  pred_and(pred_gte("year", 2013),pred_lte("year", 2022)),
  pred_within("MULTIPOLYGON(((4.29468 42.28214,4.21443 42.55702,4.02207 42.79978,3.73574 42.9861,3.38376 43.09719,3.00175 43.12179,2.62872 43.0575,2.30253 42.91094,2.05572 42.69705,1.91203 42.43731,1.8843 42.15748,1.97373 41.88491,2.17031 41.64594,2.45428 41.46342,2.79823 41.35471,3.1696 41.33013,3.53345 41.39211,3.85542 41.53488,4.10459 41.74497,4.25632 42.00233,4.29468 42.28214)))"), 
  pred("taxonKey", 212), # dowload all of class Aves
  pred_lte("coordinateUncertaintyInMeters","100"),
  format = "SIMPLE_CSV"
)

occ_download_wait(occ_es)
dwn <- occ_download_get(occ_es,path="/home/ascharf/Documents/Projects/E4Warning")
# dta_es <- occ_download_import(dwn)
dta_es <- occ_download_import(as.download("/home/ascharf/Documents/Projects/E4Warning/0006071-240314170635999.zip"))
length(unique(dta_es$species))
# 331


## get data for study site GR
st_as_text(gr) # have to paste in gbif--something changes
occ_gr <- occ_download(
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"), 
  pred_and(pred_gte("year", 2013),pred_lte("year", 2022)),
  pred_within(" MULTIPOLYGON(((25.12055 38.00168,25.12739 38.27677,25.02602 38.54018,24.82513 38.76605,24.54362 38.93188,24.2091 39.02085,23.85514 39.02373,23.51767 38.94,23.23097 38.77794,23.02378 38.55386,22.91626 38.29024,22.91807 38.01338,23.0278 37.75059,23.23356 37.5275,23.51458 37.36556,23.84346 37.28009,24.18878 37.27898,24.51789 37.36209,24.79972 37.52136,25.00739 37.74156,25.12055 38.00168)))"), 
  pred("taxonKey", 212), # dowload all of class Aves
  pred_lte("coordinateUncertaintyInMeters","100"),
  format = "SIMPLE_CSV"
)

occ_download_wait(occ_gr)
dwn <- occ_download_get(occ_gr,path="/home/ascharf/Documents/Projects/E4Warning")
# dta_gr <- occ_download_import(dwn)
dta_gr <- occ_download_import(as.download("/home/ascharf/Documents/Projects/E4Warning/0006133-240314170635999.zip"))
length(unique(dta_gr$species))
# 191