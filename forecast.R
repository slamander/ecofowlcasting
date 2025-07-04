
#'  ###########################################################################################
#'  ###########################################################################################
#'  ######## *Perform forecast seasonal dynamics of WNV and EEEV seroconversion in Florida* ###
#'  ###########################################################################################
#'  ###########################################################################################
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
#'  # To do list
#'  1. add options for which climate variables to gather using lots of if else statements
#'  2. figure out how to specify different lags to climate variables
#'  3. figure out how to specify which habitat categories to extract
#'  4. figure out how to incorporate forecasted climate variables into the pipeline!!!!
#'  
#'  
#'  ###########################################################################################
#'  ######## *1. Load necessary packages and functions* #######################################
#'  ###########################################################################################
#'  
#'  Load packages:
pacman::p_load(
  tidyverse
  , sf
  , tinyVAST
  # ,
  # ,
); conflicted::conflict_prefer_all("dplyr", quiet = T)

#'  Functions include: 
#'  - `./functions/load_data.R`: *gather and process seroconversion monitoring data*
#'  - `./functions/climate_predictors.R`: *download, process, and aggregate daymet climate**
#'  - `./functions/lulc_predictors.R`: *download nlcd, cal proportions, extract site data* 
#'  - `./functions/create_mesh.R`: *generate a spatial mesh to fit `fmesher`*

functions <- paste0("./functions/", list.files("./functions"))

sapply(functions, source, .GlobalEnv)

#'  ###########################################################################################
#'  ######## *2. Load WNV seroconversion monitoring data* #####################################
#'  ###########################################################################################
#'  
#' Monitoring data should be formatted as a `.csv` and stored in `./data/chickens/` 
#' Variables should include the following:
#'  - `Sample month`: *the month at which the sample was reported* 
#'    - formatted as a integer 1-12
#'  - `Sample year`: *the year at which the sample was reported*
#'    - formatted as a integer 2001 - present
#'  - `Testing`: *indicates if sampling occurred during the week*
#'    - formatted as binary: 1 = sampling occurred, 0 = sampling did not occur
#'  - `Viral Seroconversion result`: *indicates seroconversion status of viral sample*
#'    - possible values: NA = sampling didn't occur; 
#'    - 0 = sample was seronegative; 
#'    - 1 = sample was seropositive 
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

  
#'  ###########################################################################################
#'  ######## *3. Assemble seasonal climate predictor data* ####################################
#'  ###########################################################################################
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

#'  ###########################################################################################
#'  ######## *4. Assemble landscape predictor data* ###########################################
#'  ###########################################################################################
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
  mutate(offset = log(testing + 1)) %>%
  select(county, site, year, testing, wnv, offset,
         tmax, tmin, prcp, tmax_lag1, tmin_lag1, prcp_lag1, tmax_lag2, tmin_lag2, prcp_lag2,
         developed, cropland, natural, forest, wetlands)

data <- read_rds("data/chickens/seasonal/wnv_eeev_env_covs.rds")  %>%
  mutate_at(vars(starts_with(c("prcp", "tmax", "tmin")), "developed", "natural", "wetlands"),
            .funs = function(x){as.numeric(scale(x))}) %>%
  mutate_at(vars("season", "ID", "county"), as.factor) %>%
  mutate(season = ifelse(season == 1, "inactive", "active"),
         geometry = geometry / 1000,
         offset = log(testing + 1),
         disease = "wnv") %>%
  rename(site = ID) %>%
  bind_cols(as.data.frame(st_coordinates(.))) %>%
  select(county, site, year, testing, wnv, offset, disease,
         tmax, tmin, prcp, tmax_lag1, tmin_lag1, prcp_lag1, tmax_lag2, tmin_lag2, prcp_lag2,
         developed, natural, wetlands)

#'  ###########################################################################################
#'  ######## *5. Create spatial mesh* #########################################################
#'  ###########################################################################################
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
  data
  , resolution = c(55,53))

plot(mesh)

#'  ###########################################################################################
#'  ######## *5. Create spatial mesh* #########################################################
#'  ###########################################################################################
#'  

formula <-   wnv ~ 
  1 + 
  poly(prcp, 2) + poly(prcp_lag1, 2) + poly(tmax_lag2, 2) + poly(tmin_lag1, 2) + 
  s(county, bs="re") + 
  s(site, bs="re")

dsem <- "wnv -> wnv, 1, ar1
    wnv <-> wnv, 0, sqrt_var"

calib_data <- data %>% 
  bind_cols(
    st_coordinates(.)) %>%
  rename(lat = "Y", lon = "X") %>%
  st_drop_geometry() %>%
  droplevels() %>%
  as.data.frame()

sdmtmb_mod <- sdmTMB(
  wnv ~ 1 + poly(prcp, 2) + poly(prcp_lag1, 2) + poly(tmax_lag2, 2) + poly(tmin_lag1, 2) + 
    (1|county) + (1|site), 
  offset = "offset",
  time = "year",
  mesh = mesh,
  spatial = "off",
  spatiotemporal = "ar1",
  family = nbinom1(), 
  data = calib_data
)

tidyvast_mod <- tinyVAST(
  formula, 
  time_column = "year",
  variable_column = "disease",
  spatial_domain = mesh$mesh,
  space_term = NULL,
  spacetime_term = dsem,
  space_columns = c("lon", "lat"),
  family = nbinom1(),
  data = calib_data,
  control = tinyVASTcontrol(trace=1))


