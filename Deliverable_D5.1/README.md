##### **Deliverable_D5.1.**: Modelled monthly bird habitat suitability and species composition for the study sites

1.  We downloaded all occurrence data available on GBIF for all avian species over for 10 years for our 4 study sites. Script: [`Milestone_M5.1/DownloadAllSpsPerStudySite.R`](https://github.com/E4Warning/WP5/blob/main/Milestone_M5.1/DownloadAllSpsPerStudySite.R)

2.  For each species at each study site we modeled the habitat suitability. Scripts:

    a. `E4Warning_HSMBirds_4StudySites_MonthlyPred2022.R`: The script uses GBIF observational records to calculate habitat suitability and to create the rasters for the four study areas in four European countries (Switzerland, Germany, Spain and Greece). Each raster corresponds to the monthly predicted habitat suitability for the year 2022 for each bird species.

    b. `E4Warning_HSMBirds_4StudySites_MonthlyPred2022_SpecComp.R`: The script evaluates the created rasters to illustrate the predicted suitable bird habitats for specific months of 2022 at each study site. It also allows the user to generate a list of species, for which the rasters are available, for each grid cell freely chosen by the user.

3.  All created rasters files are available under doi_link_edmond
