
doms <- c("SEA", "AUS", "CAS", "WAS", "EAS", "AFR", "EUR", "NAM", "CAM", "SAM")

ff <- 
  map_chr(doms, function(d){
    
    str_glue("{dir_tmp}/block_max") %>% 
      list.files() %>% 
      str_subset(d) %>% 
      .[1]
    
  })

  


# import polygons into a list
l_p <- 
  
  future_map(ff %>% set_names(doms), function(f){
    
    read_ncdf(str_glue("{dir_tmp}/block_max/{f}"), 
              proxy = F) %>% 
      suppressMessages() %>% 
      suppressWarnings() %>% 
      slice(time, 1) %>%
      setNames("v") %>% 
      
      st_set_dimensions(1, values = st_get_dimension_values(., 1) %>% round(1) %>% {.-.1}) %>% 
      st_set_dimensions(2, values = st_get_dimension_values(., 2) %>% round(1) %>% {.-.1}) %>% 
      st_set_crs(4326) %>% 
      
      mutate(v = ifelse(is.na(v), NA, 1)) %>%
      st_as_sf(as_points = F, merge = T)
    
  },
  .options = furrr_options(seed = NULL))



regions_0p5 <- "regions_0p5.gpkg" %>% st_read()



covered <- 
  map(doms, function(dom){
    
    a <- st_covered_by(regions_0p5,
                       l_p %>% pluck(dom),
                       sparse = F)
    
    r <- regions_0p5 %>% 
      mutate(covered = a %>% as.numeric)
    
    return(r)
  })


covered %>% 
  bind_rows() %>%
  st_drop_geometry() %>% 
  group_by(id) %>% 
  summarise(covered = sum(covered)) %>% 
  filter(covered < 1)




