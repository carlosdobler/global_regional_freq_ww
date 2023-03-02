
library(tidyverse)
library(lubridate)
library(stars)
library(furrr)
library(units)
library(lmomRFA)


options(future.fork.enable = T)
plan(multicore)

dir_tmp <- "/mnt/pers_disk"


# load thresholds table
thresholds <- 
  str_glue("/mnt/bucket_mine/misc_data/CMIP5_model_temp_thresholds.csv") %>% 
  read_delim() %>%
  suppressMessages() %>% 
  select(1:6) %>% 
  pivot_longer(-Model, names_to = "wl") %>% 
  
  mutate(wl = str_sub(wl, 3)) %>% 
  mutate(wl = ifelse(str_length(wl) == 1, str_glue("{wl}.0"), wl))  %>%
  
  # add institutes
  mutate(Model = case_when(str_detect(Model, "HadGEM") ~ str_glue("MOHC-{Model}"),
                           str_detect(Model, "MPI") ~ str_glue("MPI-M-{Model}"),
                           str_detect(Model, "NorESM") ~ str_glue("NCC-{Model}"),
                           str_detect(Model, "GFDL") ~ str_glue("NOAA-GFDL-{Model}"),
                           str_detect(Model, "MIROC") ~ str_glue("MIROC-{Model}"),
                           TRUE ~ Model))



dom <- "NAM"

# BASELINE

# load regions 
regions_0p5 <- "regions_0p5.gpkg" %>% st_read()


# vector of files to import
ff <- 
  str_glue("{dir_tmp}/block_max") %>% 
  list.files() %>% 
  str_subset(dom)


# import files into a list
l_s <- 
  
  future_map(ff, function(f){
    
    read_ncdf(str_glue("{dir_tmp}/block_max/{f}"), 
              proxy = F) %>% 
      suppressMessages() %>% 
      suppressWarnings() %>% 
      setNames("v") %>% 
      mutate(v = set_units(v, kg/m^2/d))
    
  },
  .options = furrr_options(seed = NULL)) %>% 
  
  # fix time dim
  map(function(s){
    
    s %>%
      
      st_set_dimensions("time",
                        values = st_get_dimension_values(s, "time") %>%
                          as.character() %>%
                          str_sub(end = 4) %>% 
                          as.integer()) %>%
      
      mutate(v = set_units(v, NULL))
  })


# Verify correct import
print(str_glue("Imported:"))

walk2(l_s, ff, function(s, f){
  
  yrs <-
    s %>% 
    st_get_dimension_values("time")
  
  range_time <- 
    yrs %>% 
    range()
  
  time_steps <- 
    yrs %>% 
    length()
  
  mod <- 
    f %>% 
    str_extract("(?<=yr_)[:alnum:]*_[:graph:]*(?=\\.nc)")
  
  print(str_glue("   {mod}: \t{range_time[1]} - {range_time[2]} ({time_steps} timesteps)"))
  
})


# regularize grid

l_s <- 
  l_s %>% 
  map(function(s){
    s %>% 
      st_set_dimensions(1, values = st_get_dimension_values(s, 1) %>% round(1) %>% {.-.1}) %>% 
      st_set_dimensions(2, values = st_get_dimension_values(s, 2) %>% round(1) %>% {.-.1}) %>% 
      st_set_crs(4326)
  })





## SLICE BY WARMING LEVELS ------------------------------------------------

wls <- "0.5"

l_s_wl <- 
  
  # loop through warming levels
  map(wls, function(wl){
    
    print(str_glue("Slicing WL {wl}"))
    
    # loop through models
    map2(ff, l_s, function(f, s){
      
      # extract GCM to identify threshold year
      gcm_ <- 
        f %>% 
        str_split("_", simplify = T) %>% 
        .[,5] %>% 
        str_remove(".nc")
      
      # baseline:
      if(wl == "0.5"){
        
        s %>% 
          filter(time >= 1971,
                 time <= 2000)
        
        # other warming levels:
      } else {
        
        thres_val <-
          thresholds %>%
          filter(str_detect(Model, str_glue("{gcm_}$"))) %>% 
          filter(wl == {{wl}})
        
        s <- 
          s %>% 
          filter(time >= thres_val$value - 10,
                 time <= thres_val$value + 10)
        
        # verify correct slicing:
        print(str_glue("   {gcm_}: {thres_val$Model}: {thres_val$value}"))
        
        return(s)
      }
      
    }) %>% 
      
      # concatenate all models and form a single time dimension
      {do.call(c, c(., along = "time"))}
    
  })



s_wl <- l_s_wl[[1]]

boundary_s <- s_wl %>% 
  slice(time, 1) %>%
  mutate(v = ifelse(is.na(v), NA, 1)) %>%
  st_as_sf(as_points = F, merge = T)
  

regions <- regions_0p5



a <- st_intersects(regions,
                   boundary_s,
                   sparse = F)

regions_intersect <-
  regions[a, ] %>% 
  mutate(total_area = st_area(.) %>% set_units(NULL))
  
regions_intersect_f <- 
  
  st_intersection(regions_intersect, boundary_s) %>% 
  mutate(area_within = st_area(.) %>% set_units(NULL)) %>% 
  mutate(covered = area_within/total_area*100) %>% 
  filter(covered > 75) %>% 
  select(1)



# loop through regions

levels_baseline <- 
  future_map_dfr(seq_len(nrow(regions_intersect_f)), function(r){
    
    print(str_glue("Processing region {r}"))
    
    reg_block_max <- 
      s_wl %>% 
      st_crop(regions_intersect_f[r, ] %>% st_make_valid()) %>% 
      
    
    l_sites <-
      reg_block_max %>% 
      as_tibble() %>% 
      na.omit() %>% 
      group_by(lon, lat) %>%
      nest() 
    
    l_sites_split <- 
      l_sites %>% 
      pmap(function(data, ...){
        data %>% pull(v)
      })
    
    reg_lmom <- 
      regsamlmu(l_sites_split) # obtains l-moments for each cell
    
    reg_gev <- 
      regfit(reg_lmom, "gev") # regional gev parameters
    
    reg_precip <- 
      sitequant(0.99, reg_gev) %>% as.vector()
    
    l_sites <- 
      l_sites %>% 
      ungroup() %>% 
      select(-data) %>% 
      mutate(lev = reg_precip)
    
    return(l_sites)
    
  })


levels_baseline <- 
  levels_baseline %>% 
  mutate(lon = round(lon, 1),
         lat = round(lat, 1))


levels_baseline %>%
  mutate(lev = raster::clamp(lev, upper = quantile(lev, 0.98))) %>% 
  ggplot(aes(lon, lat, fill = lev)) +
  geom_raster() +
  colorspace::scale_fill_continuous_sequential("viridis", rev = F) +
  coord_equal() +
  theme(axis.title = element_blank())


levels_baseline <- 
  block_max %>% 
  slice(time, 1) %>%
  as_tibble() %>% 
  select(-pr) %>% 
  mutate(lon = round(lon, 1),
         lat = round(lat, 1)) %>% 
  
  left_join(levels_baseline, by = c("lon", "lat")) %>% 
  
  pull(lev) %>% 
  matrix(nrow = dim(block_max)[1], ncol = dim(block_max)[2]) %>% 
  st_as_stars()

st_dimensions(levels_baseline) <- st_dimensions(block_max %>% slice(time, 1))



# WARMING LEVELS

wl_block_max <- 
  block_max %>% 
  filter(year(time) >= 2062-10,
         year(time) <= 2062+10)


# loop through regions

levels_wl3p0 <- 
  future_map_dfr(l_sp, function(r){
    
    # print(str_glue("Processing region {r$supercells}"))
    
    reg_block_max <- 
      wl_block_max %>% 
      st_crop(r)
    
    l_sites <-
      reg_block_max %>% 
      as_tibble() %>% 
      na.omit() %>% 
      group_by(lon, lat) %>%
      nest() 
    
    l_sites_split <- 
      l_sites %>% 
      pmap(function(data, ...){
        data %>% pull(pr) %>% drop_units()
      })
    
    reg_lmom <- 
      regsamlmu(l_sites_split) # obtains l-moments for each cell
    
    reg_gev <- 
      regfit(reg_lmom, "gev") # regional gev parameters
    
    reg_precip <- 
      sitequant(0.99, reg_gev) %>% as.vector()
    
    l_sites <- 
      l_sites %>% 
      ungroup() %>% 
      select(-data) %>% 
      mutate(lon = round(lon,1),
             lat = round(lat,1))
    
    tb_baseline_lev <- 
      l_sites %>% 
      left_join(levels_baseline %>% 
                  as_tibble() %>% 
                  mutate(lon = round(lon,1),
                         lat = round(lat,1)),
                by = c("lon", "lat"))
    
    reg_prob <- 
      map_dbl(seq_len(nrow(tb_baseline_lev)), function(i){
        
        sitequant(seq(0, 1, by = 0.001), reg_gev, i) %>% 
          {which.min(abs(. - tb_baseline_lev$A1[i]))} %>% 
          names() %>% 
          as.numeric() %>% 
          {1-.}
        
      })
    
    l_sites %>% 
      mutate(lev = reg_precip,
             prob = reg_prob)
    
    
  })


# levels_wl3p0 <- 
#   levels_wl3p0 %>% 
#   mutate(lon = round(lon, 1),
#          lat = round(lat, 1))


levels_wl3p0 %>%
  mutate(lev = raster::clamp(lev, upper = quantile(lev, 0.98))) %>% 
  ggplot(aes(lon, lat, fill = lev)) +
  geom_raster() +
  colorspace::scale_fill_continuous_sequential("viridis", rev = F) +
  coord_equal() +
  theme(axis.title = element_blank())

levels_wl3p0 %>%
  mutate(prob = raster::clamp(prob, upper = quantile(prob, 0.98))) %>% 
  ggplot(aes(lon, lat, fill = prob)) +
  geom_raster() +
  colorspace::scale_fill_continuous_sequential("viridis", rev = F) +
  coord_equal() +
  theme(axis.title = element_blank())


levels_baseline <- 
  block_max %>% 
  slice(time, 1) %>%
  as_tibble() %>% 
  select(-pr) %>% 
  mutate(lon = round(lon, 1),
         lat = round(lat, 1)) %>% 
  
  left_join(levels_baseline, by = c("lon", "lat")) %>% 
  
  pull(lev) %>% 
  matrix(nrow = dim(block_max)[1], ncol = dim(block_max)[2]) %>% 
  st_as_stars()

st_dimensions(levels_baseline) <- st_dimensions(block_max %>% slice(time, 1))





