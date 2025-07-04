
climate_predictors <- function(
    state = "Florida"
    , start_year
    , end_year = format(Sys.Date(), "%Y")
    , lags = c(1, 2)
    , active_months = 6:12
    , crs = 32617
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
  ); conflicted::conflict_prefer_all("dplyr", quiet = T); conflicted::conflict_prefer_all("terra", quiet = T)  
  
  poly <- states() %>%
    filter(NAME == paste("Florida")) %>%
    st_transform(crs)
  
  dates <- tibble(
    date = seq(as.Date(paste0(2000, "-01-01")), as.Date(paste0(2025, "-12-31")), by = "day"),
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
    mutate(season = if_else(month %in% 6:12, 1, 2),
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
  
  return(
    list(prcp_seasonal, 
         tmax_seasonal, 
         tmin_seasonal)
    )
  
}
