
create_mesh <- function(
    data
    , poly
    , buffer = 2
    , poly_sample = 8000
    , resolution = c(56,51)
    , convex = -0.02
    , offset_fact1 = 1
    , offset_fact2 = 20
    , maxedge_fact = (3*5)
    , cutoff_fact = 5
    , n_knots = 100
    # ,
    # ,
    ){
  
  pacman::p_load(
    tidyverse
    , terra
    , sf
    , sdmTMB
    , INLA
    , inlabru
    , sp
    # ,
    # ,
  ); conflicted::conflict_prefer_all("dplyr", quiet = T)
  
  locs <- data %>%
    st_coordinates() %>%
    as.data.frame()
  
  polygon <- read_rds(paste0("data/polygons/", poly))
  
  poly_coords <- polygon %>%
    st_buffer(buffer) %>%
    st_sample(10000) %>%
    st_coordinates() %>%
    as.data.frame() %>%
    bind_rows(locs)
  
  domain <- poly_coords %>%
      coordinates() %>%
      inla.nonconvex.hull(convex = -.02, resolution = c(56,51))
  
  max_edge <- locs %>%
    select(X) %>%
    range() %>%
    diff()
  
  mesh_domain <- fmesher::fm_mesh_2d_inla(
    loc = locs,
    boundary = domain,
    max.edge = c(5,10)*((max_edge)/maxedge_fact),
    offset = c(max_edge/offset_fact1, max_edge/offset_fact2),
    cutoff = cutoff_fact)
  
  mesh <- make_mesh(locs, c("X", "Y"),
                    n_knots = n_knots,
                    mesh = mesh_domain)
  
  return(mesh)
  
  }