#'  #####################################################################################################################
#'  #####################################################################################################################
#'  #####################################################################################################################
#'  #####################################################################################################################
#'  ######## *Perform forecast seasonal dynamics of West Nile virus seroconversion in Florida* ##########################
#'  #####################################################################################################################
#'  #####################################################################################################################
#'  #####################################################################################################################
#'  #####################################################################################################################
#'  
#'  Directory of script: 
#'  *1. Load necessary packages*
#'  *2. Load necessary functions*
#'  *3. Load WNV seroconversion monitoring data*
#'  *4. Assemble seasonal climate predictor data*
#'  *5. Assemble landscape predictor data*
#'  *6. Create spatial mesh*
#'  
#'  #####################################################################################################################
#'  ######## *1. Load necessary packages* ###############################################################################
#'  #####################################################################################################################
#'  
#'  Load packages:
pacman::p_load(
  tidyverse
  , sf
  # ,
  # ,
); conflicted::conflict_prefer_all("dplyr", quiet = T)

#'  #####################################################################################################################
#'  ######## *2. Load necessary functions* ##############################################################################
#'  #####################################################################################################################
#'  
#'  Functions include: 
#'  - `load_data.R`: *gather and process WNV seroconversion monitoring data by aggregating to the seasonal scale*
#'  - `create_mesh.R`: *generate a spatial mesh to fit the stochastic partial differential equations with `fmesher`*
#'  - `climate_predictors.R`: *download, rasterize, and aggregate daymet climate variables* 

functions <- paste0("./functions/", list.files("./functions"))

sapply(functions, source, .GlobalEnv)

#'  #####################################################################################################################
#'  ######## *3. Load WNV seroconversion monitoring data* ############################################################### 
#'  #####################################################################################################################
#'  
#' Monitoring data should be formatted as a `.csv` and stored in `./data/chickens/` and include the following variables
#'  - `Sample month`: *the month at which the sample was reported* 
#'    - formatted as a integer 1-12
#'  - `Sample year`: *the year at which the sample was reported*
#'    - formatted as a integer 2001 - present
#'  - `Testing`: *indicates if sampling occurred during the week*
#'    - formatted as binary: 1 = sampling occurred, 0 = sampling did not occur
#'  - `WNV Seroconversion result`: *indicates seroconversion status of WNV sample*
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

wnv <- load_data(
  filename = "compiled_chicken_data_clean_v2.csv", 
  month = month, 
  year = year, 
  testing = testing, 
  wnv = n_positive_wnv, 
  active_months = 6:12,
  lat = lat, 
  lon = lon, 
  county = county, 
  site = ID, 
  crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
)

#'  #####################################################################################################################
#'  ######## *4. Assemble seasonal climate predictor data* ############################################################## 
#'  #####################################################################################################################
#'  
#'  Calculate seasonal climate variables:

clim <- climate_predictors(
  state = "Florida"
  , start_year = data
  , end_year = format(Sys.Date(), "%Y")
  , lags = c(1, 2)
  , active_months = 6:12
  , crs = 32617
)


#'  #####################################################################################################################
#'  ######## *5. Create spatial mesh* ################################################################################### 
#'  #####################################################################################################################
#'  
#'  Arguments for creating mesh: 
#'  - data
#'  - poly
#'  - buffer = 2
#'  - poly_sample = 8000
#'  - resolution = c(56,51)
#'  - convex = -0.02
#'  - offset_fact1 = 1
#'  - offset_fact2 = 20
#'  - maxedge_fact = (3*5)
#'  - cutoff_fact = 5
#'  - n_knots = 100
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