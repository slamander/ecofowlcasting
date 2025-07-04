
lulc_predictors <- function(
    template = "data/utils/raster.tif"
    , data
    , start_year
    , end_year = as.integer(format(Sys.Date(), "%Y"))
    , buffer = 2500
    , return_rasters = F
    # , 
    # , 
){
  
  pacman::p_load(
    tidyverse
    , terra
    , sf
    , FedData
    , pbapply
    # ,
    # ,
  ); 
  conflicted::conflict_prefer_all("dplyr", quiet = T); 
  conflicted::conflict_prefer_all("terra", quiet = T)  
  
  template <- rast(paste(template))
  
  nlcd_years <- c(2001, 2004, 2006, 2008, 2011, 2016, 2019)
  
  raster_mapping <- data.frame(
    year = start_year:end_year,
    nlcd_year = nlcd_years[sapply(start_year:end_year, FUN = function(x){
      which.min(abs(x-nlcd_years))
    })]) %>%
    mutate(layer = row_number())
  
  ## Download NLCD data
  
  nlcd_list <- pblapply(unique(raster_mapping$nlcd_year), 
                        FUN = function(x){
                          get_nlcd(
                            template = template, 
                            "nlcd",
                            year = x)})
  
  nlcd_raw <- rast(nlcd_list)
  names(nlcd_raw) <- unique(raster_mapping$nlcd_year)
  
  nlcd_annual <- lapply(raster_mapping$nlcd_year, 
                        FUN = function(x){
                          nlcd_raw[[paste(x)]]}) %>%
    unlist() %>%
    rast()
  names(nlcd_annual) <- raster_mapping$year
  
  cat('Writing annual landuse landcover data to `./data/lulc/`')
  
  writeRaster(nlcd_annual, "data/lulc/nlcd_annual.tif", overwrite = T)
  
  data_years <- data %>%
    st_drop_geometry() %>%
    select(year) %>%
    unique() %>%
    arrange(year) %>%
    pull(year)
  
  nlcd_extract <- data_list <- list()
  
  for(i in 1:length(data_years)){

    data_list[[paste(data_years[i])]] <- data %>%
      filter(year == as.integer(data_years[i]))  %>%
      group_by(ID) %>%
      summarize(years = n()) %>%
      st_buffer(buffer)
    
    nlcd_extract[[paste(data_years[i])]] <- terra::extract(
      nlcd_annual[[paste(data_years[i])]],
      data_list[[i]]) %>%
      table() %>%
      as.data.frame.matrix() %>%
      bind_cols(st_drop_geometry(data_list[[i]])) %>%
      dplyr::relocate(ID)
    
    cat("\r Extracting lulc at monitoring sites: ", 
        round(
          (i/length(data_years)*100), 0), 
        "% complete", 
        sep = "")
    
  }
  
  nlcd_prop <- do.call(rbind, nlcd_extract) %>%
    mutate(total = rowSums(select_if(., is.numeric), na.rm = TRUE)) %>%
    clean_names() %>%
    mutate(developed = (rowSums(
      dplyr::select(., 
                    developed_open_space,
                    developed_low_intensity,
                    developed_medium_intensity,
                    developed_high_intensity))/total) %>% round(5)) %>%
    
    mutate(cropland = (rowSums(
      dplyr::select(., 
                    pasture_hay,
                    cultivated_crops))/total) %>% round(5)) %>%
    
    mutate(natural = (rowSums(
      dplyr::select(., 
                    deciduous_forest,
                    evergreen_forest,
                    mixed_forest,
                    shrub_scrub,
                    sedge_herbaceous,
                    grassland_herbaceous))/total) %>% round(5)) %>%
    
    mutate(forest = (rowSums(
      dplyr::select(., 
                    deciduous_forest,
                    evergreen_forest,
                    mixed_forest))/total) %>% round(5)) %>%
    
    mutate(wetlands = (rowSums(
      dplyr::select(., 
                    woody_wetlands,
                    emergent_herbaceous_wetlands))/total) %>% round(5)) %>%
    
    rownames_to_column(var = "year") %>%
    mutate(year = as.numeric(year) %>% floor()) %>%
    dplyr::select(id, year, developed, cropland, natural, forest, wetlands)
  
  data_nlcd <- data %>%
    left_join(nlcd_prop, by = c("site" = "id", "year" = "year")) 
  
  # if(return_rasters == T){
  #   
  #   return(list(data_nlcd, 
  #               raster object))
  #   
  # }
  
  return(data_nlcd)
    
}