
load_data <- function(
    filename, 
    month, 
    year, 
    testing, 
    wnv, 
    active_months,
    lat, 
    lon, 
    county, 
    site, 
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", 
    transform = 32617
    # ,
    # ,
    ){
  
  pacman::p_load(
    tidyverse
    , sf
    , lubridate
    # ,
    # ,
    ); conflicted::conflict_prefer_all("dplyr", quiet = T)
  
  # load raw seroconversion monitoring data: 
  raw <- read_csv(paste0("data/chickens/", filename)) 
  
  # process raw data: 
  # - filtering non-georeferenced observations
  # - establishing the active season
  # - aggregating observations to season by summing WNV seroresults and sampling effort
  
  data <- raw %>%
    rename("month" = {{month}}, 
           "year" = {{year}}, 
           "testing" = {{testing}}, 
           "wnv" = {{wnv}}, 
           "lat" = {{lat}}, 
           "lon" = {{lon}}, 
           "county" = {{county}}, 
           "ID" = {{site}}) %>%
    filter(!is.na(lat)) %>%
    filter(!is.na(lon)) %>%
    mutate(season = if_else(month %in% active_months, 1, 2),
           season_name = if_else(season == 1, "active", "inactive")) %>%
    group_by(county, ID, year, month, season, lat, lon) %>%
    summarize(
      wnv = sum(wnv, na.rm = T),
      testing = sum(testing, na.rm = T)) %>%
    ungroup() %>%
    st_as_sf(x = .,
             coords = c("lon", "lat"),
             crs = crs) %>%
    bind_cols(
      st_coordinates(.)) %>%
    st_transform(32617) %>%
    mutate(date = make_date(year = year, month = month)) %>%
    rename(lat = "Y", lon = "X") %>%
    filter(season == 1) %>%
    select(county, ID, date, year, month, lat, lon, wnv, testing, geometry)
  
  return(data)

}
