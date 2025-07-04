#'  #####################################################################################################################
#'  #####################################################################################################################
#'  #####################################################################################################################
#'  #####################################################################################################################
#'  ######## *Perform forecast seasonal dynamics of WNV and EEEV seroconversion in Florida* #############################
#'  #####################################################################################################################
#'  #####################################################################################################################
#'  #####################################################################################################################
#'  #####################################################################################################################
#'  
#'  # Directory of script: 
#'  *1. Load necessary packages and functions*: 
#'  *2. Load viral seroconversion monitoring data*
#'    - `./functions/load_data.R`:
#'  *3. Assemble seasonal climate predictor data*
#'    - `./functions/climate_predictors.R`: 
#'  *4. Assemble landscape predictor data*
#'    - `./functions/lulc_predictors.R`: 
#'  *5. Create spatial mesh*
#'    - `./functions/create_mesh.R`
#'  
#'  
#'  # To do list
#'  1. add options for which climate variables to gather using lots of if else statements
#'  2. figure out how to specify different lags to climate variables
#'  3. figure out how to specify which habitat categories to extract
#'  4. figure out how to incorporate forecasted climate variables into the pipeline!!!!
#'  
#'  
#'    
#'  #####################################################################################################################
#'  ######## *1. Load necessary packages and functions* ###############################################################################
#'  #####################################################################################################################
#'  
#'  Load packages:
pacman::p_load(
  tidyverse
  , sf
  # ,
  # ,
); conflicted::conflict_prefer_all("dplyr", quiet = T)

#'  Functions include: 
#'  - `./functions/load_data.R`: *gather and process seroconversion monitoring data*
#'  - `./functions/climate_predictors.R`: *download, rasterize, and aggregate daymet climate variables* 
#'  - `./functions/lulc_predictors.R`: *download nlcd data, calculate cover proportions, extract site-level data* 
#'  - `./functions/create_mesh.R`: *generate a spatial mesh to fit `fmesher`*

functions <- paste0("./functions/", list.files("./functions"))

sapply(functions, source, .GlobalEnv)

#'  #####################################################################################################################
#'  ######## *2. Load WNV seroconversion monitoring data* ############################################################### 
#'  #####################################################################################################################
#'  
#' Monitoring data should be formatted as a `.csv` and stored in `./data/chickens/` and include the following variables
#'  - `Sample month`: *the month at which the sample was reported* 
#'    - formatted as a integer 1-12
#'  - `Sample year`: *the year at which the sample was reported*
#'    - formatted as a integer 2001 - present
#'  - `Testing`: *indicates if sampling occurred during the week*
#'    - formatted as binary: 1 = sampling occurred, 0 = sampling did not occur
#'  - `Viral Seroconversion result`: *indicates seroconversion status of viral sample*
#'    - possible values: NA = sampling didn't occur; 0 = sample was seronegative; 1 = sample was seropositive 
#'  - `Latitude`: *the latitude of the monitoring site*
#'    - formatted as a numeric 
#'  - `Longitude`: *the latitude of the monitoring site*
#'    - formatted as a numeric
#'  - `County`: *the county of the monitoring site*
#'    - formatted as a character
#'  - `Site`: *the county of the monitoring site*
#'    - formatted as a character
#'  - `Coordinate reference system`: *the crs of the desired spatial data*
#'    - formatted as a character string
#'  
#'  Load and process WNV monitoring data:

wnv_monthly <- load_data(
  filename = "compiled_chicken_data_clean_v2.csv", 
  wnv = n_positive_wnv, 
  active_months = 6:12,
  site = ID
)

wnv_seasonal <- wnv_monthly %>% 
  st_drop_geometry() %>%
  group_by(county, ID, year, lat, lon) %>%
  summarize(
    wnv = sum(wnv, na.rm = T),
    testing = sum(testing, na.rm = T)) %>%
  ungroup() %>%
  st_as_sf(x = .,
           coords = c("lon", "lat"),
           crs = crs) %>%
  bind_cols(
    st_coordinates(.)) %>%
  rename(lat = "Y", lon = "X") %>%
  select(county, ID, date, year, lat, lon, wnv, geometry)

  
#'  #####################################################################################################################
#'  ######## *3. Assemble seasonal climate predictor data* ############################################################## 
#'  #####################################################################################################################
#'  
#'  Calculate seasonal climate variables:
#'  `./functions/climate_predictors.R`:
#'  Arguments:
#'  - `monthly_data`, 
#'  - `seasonal_data`, 
#'  - `state` = "Florida", 
#'  - `start_year`,
#'  - `end_year` = as.numeric(format(Sys.Date(), %Y")), 
#'  - `lags` = c(1, 2), 
#'  - `active_months` = 6:12, 
#'  - `crs` = 32617, 
#'  - `return_rasters` = F, 
#'  - `return_monthly_extract` = F, 
#'  - `return_seasonal_extract` = T)
#'  

wnv_clim <- climate_predictors(
  monthly_data = wnv_monthly
  , seasonal_data = wnv_seasonal
  , start_year = 2001
)




#'  #####################################################################################################################
#'  ######## *4. Assemble landscape predictor data* ##################################################################### 
#'  #####################################################################################################################
#'  
#'  Gather and extract land use land cover from NLCD
#'  `./functions/lulc_predictors.R`:
#'  Arguments:
#'  - `template` = "data/utils/raster.tif"
#'  - `data`
#'  - `start_year`
#'  - `end_year` = as.integer(format(Sys.Date(), "%Y"))
#'  - `buffer` = 2500
#'  - `return_rasters` = F
#'  

lulc <- lulc_predictors(
  data = wnv_clim[["seasonal_extract"]]
  , start_year = 2021
)

data <- lulc %>% 
  mutate(offset = log(testing + 1)),
  select(county, site, year, testing, wnv, offset,
         tmax, tmin, prcp, tmax_lag1, tmin_lag1, prcp_lag1, tmax_lag2, tmin_lag2, prcp_lag2,
         developed, cropland, natural, forest, wetlands)

#'  #####################################################################################################################
#'  ######## *5. Create spatial mesh* ################################################################################### 
#'  #####################################################################################################################
#'  
#'  Arguments for creating mesh: 
#'  - `data`
#'  - `poly`
#'  - `buffer = 2`
#'  - `poly_sample = 8000`
#'  - `resolution = c(56,51)`
#'  - `convex = -0.02`
#'  - `offset_fact1 = 1`
#'  - `offset_fact2 = 20`
#'  - `maxedge_fact = (3*5)`
#'  - `cutoff_fact = 5`
#'  - `n_knots = 100`
#'  
#'  Create spatial mesh:

mesh <- create_mesh(
  data = data, 
  poly = "poly.rds")




















#'  #####################################################################################################################
#'  ##################################################### fit model ##################################################### 
#'  #####################################################################################################################
#'  
#'  