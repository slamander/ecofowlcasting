
climate_predictors <- function(
    monthly_data
    , seasonal_data
    , state = "Florida"
    , start_year
    , end_year = as.numeric(format(Sys.Date(), "%Y"))
    , lags = c(1, 2)
    , active_months = 6:12
    , crs = 32617
    , return_rasters = F
    , return_monthly_extract = F
    , return_seasonal_extract = T
    # , 
    # , 
    ){

  pacman::p_load(
    tidyverse
    , terra
    , sf
    , lubridate
    , sp
    # ,
    # ,
  ); 
  conflicted::conflict_prefer_all("dplyr", quiet = T); 
  conflicted::conflict_prefer_all("terra", quiet = T)  
  
  poly <- states() %>%
    filter(NAME == paste("Florida")) %>%
    st_transform(crs)
  
  dates <- tibble(
    date = seq(as.Date(paste0(start_year, "-01-01")), as.Date(paste0(end_year, "-12-31")), by = "day"),
    year = year(date),
    month = month(date), 
    day = day(date)) %>%
    group_by(year, month) %>%
    summarize(
      start = min(date), 
      end = max(date)) %>%
    mutate(year = year(start),
           month = month(start),
           month_year = paste(month, year, sep = "-"), 
           season = if_else(month %in% active_months, 1, 2),
           season_name = if_else(season == 1, "active", "inactive"))
  
  prcp <- tmax <- tmin <- list()
  
  for(i in 1:nrow(dates)){

    prcp[[paste(dates$month_year[i])]] <- getDaymet(AOI = poly, 
                                                    varname = "prcp",
                                                    startDate = paste(dates$start[i]),
                                                    endDate = paste(dates$end[i]),
                                                    dryrun = F)[[1]] %>%
      sum()

    tmax[[paste(dates$month_year[i])]] <- getDaymet(AOI = poly,
                                                    varname = "tmax",
                                                    startDate = paste(dates$start[i]),
                                                    endDate = paste(dates$end[i]),
                                                    dryrun = F)[[1]] %>%
      mean()

    tmin[[paste(dates$month_year[i])]] <- getDaymet(AOI = poly,
                                                    varname = "tmin",
                                                    startDate = paste(dates$start[i]),
                                                    endDate = paste(dates$end[i]),
                                                    dryrun = F)[[1]] %>%
      mean()
    
    cat("\r processing monthly climate: ", round((i/rnow(dates))*100, 0), "% complete", sep = "")
    
  }
  
  prcp_monthly <- rast(prcp) %>%
    project(poly, res = 1000) %>%
    crop(poly, mask = T)
  
  tmax_monthly <- rast(tmax) %>%
    project(poly, res = 1000) %>%
    crop(poly, mask = T)
  
  tmin_monthly <- rast(tmin) %>%
    project(poly, res = 1000) %>%
    crop(poly, mask = T)
  
  raster_map <- data.frame(
    month = str_split_fixed(names(prcp_monthly), "-", 2)[,1],
    year = str_split_fixed(names(prcp_monthly), "-", 2)[,2] %>% as.integer()) %>%
    mutate(season = if_else(month %in% active_months, 1, 2),
           season_name = if_else(season == 1, "active", "inactive"),
           layer = row_number())
  
  years <- unique(raster_map$year)
  
  raster_annual_map <- raster_map %>%
    group_by(year, season) %>%
    summarize(year = unique(year),
              season = unique(season)) %>%
    ungroup() %>%
    mutate(season_name = if_else(season == 1, "active", "inactive"),
           year = as.character(year),
           season = as.character(season),
           lag0 = row_number(),
           lag1 = lag0 - lags[1], 
           lag2 = lag0 - lags[2])
  
  prcp_seasonal <- tmax_seasonal <- tmin_seasonal <- list()
  
  for(i in 1:length(years)){
    
    active <- raster_map$layer[raster_map$year == years[i] & raster_map$season == 1]
    
    prcp_seasonal[[paste0("year_", years[i])]] <- sum(prcp_monthly[[active]], na.rm = T)
    tmax_seasonal[[paste0("year_", years[i])]] <- mean(tmax_monthly[[active]], na.rm = T)
    tmin_seasonal[[paste0("year_", years[i])]] <- mean(tmin_monthly[[active]], na.rm = T)
    
    cat("\r processing seasonal climate: ", round((i/length(years))*100, 0), "% complete", sep = "")
    
  }
  
  prcp_seasonal <- rast(prcp_seasonal)
  names(prcp_seasonal) <- years
  
  tmax_seasonal <- rast(tmax_seasonal)
  names(tmax_seasonal) <- years
  
  tmin_seasonal <- rast(tmin_seasonal)
  names(tmin_seasonal) <- years
  
  cat('Writing monthly and seasonal climate data to `./data/climate/`')
  
  writeRaster(prcp_monthly, "data/environment/climate/prcp_stack_monthly.tif", overwrite = T)
  writeRaster(tmax_monthly, "data/environment/climate/tmax_stack_monthly.tif", overwrite = T)
  writeRaster(tmin_monthly, "data/environment/climate/tmin_stack_monthly.tif", overwrite = T)
  
  writeRaster(prcp_seasonal, "data/environment/climate/prcp_stack_seasonal.tif", overwrite = T)
  writeRaster(tmax_seasonal, "data/environment/climate/tmax_stack_seasonal.tif", overwrite = T)
  writeRaster(tmin_seasonal, "data/environment/climate/tmin_stack_seasonal.tif", overwrite = T)
  
  # extract site-level estimates
  
  ### creating lagged versions of the chicken data:
  
  cat('Extracting monthly climate to site-level data')
  
  data_monthly_lag <- monthly_data %>%
    mutate(date_lag1 = make_date(year, month) %m-% months(1),
           month_lag1 = month(date_lag1),
           year_lag1 = year(date_lag1),
           date_lag2 = make_date(year, month) %m-% months(2),
           month_lag2 = month(date_lag2),
           year_lag2 = year(date_lag2)) 
  
  sites <- monthly_data %>%
    select(ID, lon, lat) %>%
    unique()
  
  prcp_extract <- extract(prcp_monthly, sites, method = "bilinear", xy = T, ID = F) %>%
    bind_cols(ID = sites %>%
                select(ID, lon, lat)) %>%
    pivot_longer(
      cols = matches("\\d+-\\d{4}"),         
      names_to = c("month", "year"),         
      names_sep = "-",                       
      names_transform = list(                
        month = as.integer,
        year = as.integer),
      values_to = "prcp")
  
  tmax_extract <- extract(tmax_monthly, sites, method = "bilinear", xy = T, ID = F) %>%
    bind_cols(ID = sites %>%
                select(ID, lon, lat)) %>%
    pivot_longer(
      cols = matches("\\d+-\\d{4}"),         
      names_to = c("month", "year"),         
      names_sep = "-",                       
      names_transform = list(                
        month = as.integer,
        year = as.integer),
      values_to = "tmax")
  
  tmin_extract <- extract(tmin_monthly, sites, method = "bilinear", xy = T, ID = F) %>%
    bind_cols(ID = sites %>%
                select(ID, lon, lat)) %>%
    pivot_longer(
      cols = matches("\\d+-\\d{4}"),         
      names_to = c("month", "year"),         
      names_sep = "-",                       
      names_transform = list(                
        month = as.integer,
        year = as.integer),
      values_to = "tmin")
  
  clim_extract <- prcp_extract %>%
    bind_cols(tmax = tmax_extract$tmax, 
              tmin = tmin_extract$tmin) %>%
    st_drop_geometry() %>%
    select(-c("geometry"))
  
  data_monthly_clim <- data_monthly_lag %>%
    left_join(clim_extract, 
              by = c("ID", "month", "year"))
  
  data_monthly_clim_lag1 <- data_monthly_lag %>%
    left_join(clim_extract, 
              by = c("ID", "month_lag1" = "month", "year_lag1" = "year")) %>%
    rename(tmax_lag1 = "tmax", 
           tmin_lag1 = "tmin", 
           prcp_lag1 = "prcp") %>%
    st_drop_geometry()
  
  data_monthly_clim_lag2 <- data_monthly_lag %>%
    left_join(clim_extract, 
              by = c("ID", "lat", "lon", "month_lag2" = "month", "year_lag2" = "year")) %>%
    rename(tmax_lag2 = "tmax", 
           tmin_lag2 = "tmin", 
           prcp_lag2 = "prcp") %>%
    st_drop_geometry()
  
  data_monthly_clim_lags <- data_monthly_clim %>%
    left_join(data_monthly_clim_lag1) %>%
    left_join(data_monthly_clim_lag2) %>%
    select(county, ID, lat, lon, year, month, season, testing, wnv, eeev, 
           tmax, tmin, prcp,
           tmax_lag1, tmin_lag1, prcp_lag1,
           tmax_lag2, tmin_lag2, prcp_lag2,
           geometry) 
  
  data_monthly_clim_lags %>%
    st_drop_geometry() %>%
    filter(season == "2") %>%
    group_by(ID, season) %>%
    summarize(n = n(), 
              wnv = sum(wnv, na.rm = T)) %>%
    as.data.frame()
  
  cat('Writing monthly climate extractions as `./data/chickens/monthly/wnv_eeev_clim.rds`')
  
  write_rds(data_monthly_clim_lags, "data/chickens/monthly/wnv_eeev_clim.rds")

  ########
  
  cat('Extracting seasonal climate to site-level data')

  ### creating lagged versions of the seasonal chicken data:
  
  data_seasonal_lag <- seasonal_data %>%
    mutate(date_lag1 = make_date(year, month) %m-% months(6),
           month_lag1 = month(date_lag1),
           year_lag1 = year(date_lag1),
           date_lag2 = make_date(year, month) %m-% months(12),
           month_lag2 = month(date_lag2),
           year_lag2 = year(date_lag2)) %>%
    select(county, ID, lat, lon, 
           year, month, testing,
           year_lag1, month_lag1, 
           year_lag2, month_lag2, 
           season, wnv, eeev, geometry); data_seasonal_lag
  
  data_seasonal_clim_lag1 <- data_seasonal_lag %>%
    left_join(clim_extract, 
              by = c("ID", "lat", "lon", "month_lag1" = "month", "year_lag1" = "year")) %>%
    rename(tmax_lag1 = "tmax", 
           tmin_lag1 = "tmin", 
           prcp_lag1 = "prcp") %>%
    st_drop_geometry()
  
  data_seasonal_clim_lag2 <- data_seasonal_lag %>%
    left_join(clim_extract, 
              by = c("ID", "lat", "lon", "month_lag2" = "month", "year_lag2" = "year")) %>%
    rename(tmax_lag2 = "tmax", 
           tmin_lag2 = "tmin", 
           prcp_lag2 = "prcp") %>%
    st_drop_geometry()
  
  data_seasonal_clim_lags <- data_monthly_clim %>%
    left_join(data_seasonal_clim_lag1, by = c("county", "ID", "year", "month", "season", "testing", "wnv", "eeev")) %>%
    left_join(data_seasonal_clim_lag2, by = c("county", "ID", "year", "month", "season", "testing", "wnv", "eeev")) %>%
    rename("lat" = "lat.x", "lon" = "lon.x") %>%
    group_by(county, ID, lat, lon, year, season) %>%
    summarize(
      testing = sum(testing), 
      wnv = sum(wnv), 
      eeev = sum(eeev), 
      tmax = mean(tmax), 
      tmin = mean(tmin), 
      prcp = sum(prcp),
      tmax_lag1 = mean(tmax_lag1), 
      tmin_lag1 = mean(tmin_lag1), 
      prcp_lag1 = sum(prcp_lag1),
      tmax_lag2 = mean(tmax_lag2), 
      tmin_lag2 = mean(tmin_lag2), 
      prcp_lag2 = sum(prcp_lag2)) %>%
    ungroup() %>%
    select(county, ID, lat, lon, year, season, testing, wnv, eeev, 
           tmax, tmin, prcp,
           tmax_lag1, tmin_lag1, prcp_lag1,
           tmax_lag2, tmin_lag2, prcp_lag2,
           geometry) 
  
  data_seasonal_clim_lags %>%
    st_drop_geometry() %>%
    filter(season == "2") %>%
    group_by(ID, season) %>%
    summarize(n = n(), 
              wnv = sum(wnv, na.rm = T)) %>%
    as.data.frame()
  
  cat('Writing seasonal climate extractions as `./data/chickens/seasonal/wnv_eeev_clim.rds`')
  
  write_rds(data_seasonal_clim_lags, "data/chickens/seasonal/wnv_eeev_clim.rds")

  # return objects
  
  if(return_rasters == T & return_monthly_extract == T & return_monthly_extract == T){
    
    return <- list()
    return[["monthly_rasters"]] <- list(prcp_monthly, 
                                        tmax_monthly, 
                                        tmin_monthly)
    return[["seasonal_rasters"]] <- list(prcp_seasonal, 
                                         tmax_seasonal, 
                                         tmin_seasonal)
    return[["seasonal_extract"]] <- data_seasonal_clim_lags
    return[["monthly_extracts"]] <- data_monthly_clim_lags
    
    return(return)
    
  }
  
  if(return_rasters == T & return_seasonal_extract == F & return_monthly_extract == T){
    
    return <- list()
    return[["monthly_rasters"]] <- list(prcp_monthly, 
                                        tmax_monthly, 
                                        tmin_monthly)
    return[["seasonal_rasters"]] <- list(prcp_seasonal, 
                                         tmax_seasonal, 
                                         tmin_seasonal)
    return[["monthly_extract"]] <- data_monthly_clim_lags

    return(return)
    
  }
  
  if(return_rasters == T & return_seasonal_extract == T & return_monthly_extract == F){
    
    return <- list()
    return[["monthly_rasters"]] <- list(prcp_monthly, 
                                        tmax_monthly, 
                                        tmin_monthly)
    return[["seasonal_rasters"]] <- list(prcp_seasonal, 
                                         tmax_seasonal, 
                                         tmin_seasonal)
    return[["seasonal_extract"]] <- data_seasonal_clim_lags

    return(return)
    
  }
  
  if(return_rasters == F & return_seasonal_extract == T & return_monthly_extract == T){
    
    return <- list()
    return[["seasonal_extract"]] <- data_seasonal_clim_lags
    return[["monthly_extracts"]] <- data_monthly_clim_lags
    
    return(return)
    
  }
  
  if(return_rasters == F & return_seasonal_extract == T & return_monthly_extract == F){
    
    return <- list()
    return[["seasonal_extract"]] <- data_seasonal_clim_lags

    return(return)
    
  }
  
  if(return_rasters == F & return_seasonal_extract == F & return_monthly_extract == T){
    
    return <- list()
    return[["monthly_extracts"]] <- data_monthly_clim_lags
    
    return(return)
    
  }
  
}
  

