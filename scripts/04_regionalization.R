

library(tidyverse)
library(lubridate)
library(stars)
library(furrr)
library(units)
library(supercells)
library(rgeoda)
library(regional)

options(future.fork.enable = T)
plan(multicore)

dir_tmp <- "/mnt/pers_disk"



# load block max quintiles
block_max_quintiles <- 
  str_glue("{dir_tmp}/mos_blockmax_quint.nc") %>% 
  read_ncdf()


# prepare land mask

# rast_reference_0.05 <-
#   st_as_stars(st_bbox(block_max_quintiles), dx = 0.05, dy = 0.05, values = -9999)
# 
# land <-
#   "/mnt/bucket_mine/misc_data/ne_50m_land/ne_50m_land.shp" %>%
#   st_read() %>%
#   mutate(a = 1) %>%
#   select(a) %>%
#   st_rasterize(rast_reference_0.05)
# 
# land <-
#   land %>%
#   # st_warp(block_max_quintiles %>% slice(wl, 1), 
#   #         use_gdal = T, 
#   #         method = "mode") %>% 
#   st_warp(block_max_quintiles %>% slice(wl, 1)) %>% 
#   setNames("a") %>% 
#   mutate(a = ifelse(a == -9999, 0, 1))

land <- 
  "/mnt/bucket_cmip5/Probable_futures/irunde_scripts/create_a_dataset/04_rcm_buffered_ocean_mask.nc" %>% 
  read_ncdf() %>%
  st_warp(block_max_quintiles %>% slice(wl, 1)) %>% 
  setNames("a")

# mask

# block_max_quintiles[land == 0] <- -9999
# 
# block_max_quintiles <-
#   block_max_quintiles %>%
#   mutate(bp1 = ifelse(is.na(bp1), -9999, bp1),
#          bp2 = ifelse(is.na(bp2), -9999, bp2),
#          bp3 = ifelse(is.na(bp3), -9999, bp3),
#          bp4 = ifelse(is.na(bp4), -9999, bp4),
#          bp5 = ifelse(is.na(bp5), -9999, bp5))

# block_max_quintiles[is.na(land)] <- NA



# BASELINE ********************************************************************



block_max_bl <- block_max_quintiles %>% slice(wl, 1)
block_max_bl[is.na(land)] <- NA
block_max_bl <- block_max_bl %>% as("SpatRaster")

# %>% as("SpatRaster")
# block_max_bl[block_max_bl == -9999] <- NA



# step 1: micro-regions (supercells)

{
  
  param_grid_1 <- expand_grid(step_ = seq(5, 15, by = 2),
                              compactness_ = seq(10, 70, by = 20))

  tb_inh_1 <-
    param_grid_1 %>%
    mutate(inh = future_pmap(., function(step_, compactness_){

      sp <-
        supercells(block_max_bl,
                   step = step_,
                   compactness = compactness_)

      inh <-
        reg_inhomogeneity(sp,
                          block_max_bl,
                          sample_size = 0.5)

      return(inh)


    }))


  tb_inh_1 %>%
    unnest(inh) %>%
    mutate(compactness_ = factor(compactness_)) %>%

    {
      ggplot(., aes(x = compactness_, group = compactness_, fill = compactness_, y = inh)) +
        geom_boxplot(outlier.shape = NA) +
        facet_grid(~step_) +
        scale_y_continuous(limits = quantile(.$inh, c(0.1, 0.9))) +
        guides(fill = "none") +
        labs(subtitle = "step_")
    }
  
  
  
}


# block_max_bl <- block_max_quintiles %>% slice(wl, 1)
# block_max_bl[is.na(land)] <- -9999
# block_max_bl <-
#   block_max_bl %>%
#   mutate(bp1 = ifelse(is.na(bp1), -9999, bp1),
#          bp2 = ifelse(is.na(bp2), -9999, bp2),
#          bp3 = ifelse(is.na(bp3), -9999, bp3),
#          bp4 = ifelse(is.na(bp4), -9999, bp4),
#          bp5 = ifelse(is.na(bp5), -9999, bp5))
# block_max_bl <- block_max_bl %>% as("SpatRaster")


sp_bl <- 
  supercells(block_max_bl,
             step = 5,
             compactness = 10)


# sp_bl %>% select(6) %>% mapview::mapview()




# step 2: aggregation

# APPROACH 1 KMEANS

# library(tidymodels)
# 
# sp_bl %>%
#   mutate(across(bp1:bp5, ~scale(.x))) %>% 
#   .$bp1 %>% 
#   quantile()
# 
# 
# sp_bl %>% 
#   st_drop_geometry() %>%
#   select(bp1:bp5) %>% 
#   scale(center = F, scale = T) %>% 
#   .[,1] %>% 
#   quantile()
# 
# m <- sp_bl %>% 
#   st_drop_geometry() %>%
#   select(x:bp5) %>% 
#   scale(center = F, scale = T)
#   
# sp_bl_m <- 
#   sp_bl %>% 
#   select(1) %>% 
#   bind_cols(m) %>% 
#   mutate(across(x:y, ~.x*2))
#   
# 
# clus <- 
#   kmeans(sp_bl_m[2:8] %>% st_drop_geometry(), centers = 100, nstart = 200)
# 
# 
# 
# sp_bl_m$cluster <- clus$cluster
# 
# sp_bl_c <- 
#   sp_bl_m %>% 
#   group_by(cluster) %>% 
#   summarize()
# 
# sp_bl_c %>% mapview::mapview()


# APPROACH 2
# fill holes -> skater

land_p <- 
  c(land,
    block_max_quintiles %>% slice(wl, 1)) %>% 
  mutate(a = ifelse(is.na(bp1), NA, a)) %>% 
  
  select(a) %>% 
  st_as_sf(as_points = F, merge = T) %>% 
  
  mutate(id = row_number(),
         area = st_area(.),
         area = area %>% set_units(NULL) %>% {./1000000000})
  

land_p_small <- 
  land_p %>% 
  filter(area < 700) %>% 
  select(geometry)
  
land_p_large <- 
  land_p %>% 
  filter(area > 700)
  


block_max_bl <- block_max_quintiles %>% slice(wl, 1)
  

clus_assess <- 
  map(seq_len(nrow(land_p_large)), function(i){
    
    print(str_glue("super-region {i} / {nrow(land_p_large)}"))
    
    pol <- land_p_large[i, ]
    
    conv_hull <-
      pol %>% 
      st_boundary() %>% 
      st_polygonize() %>% 
      st_rasterize(block_max_quintiles %>% 
                     slice(wl, 1) %>% 
                     select(bp1) %>% 
                     mutate(bp1 = NA))
    
    b <- block_max_bl
    b[is.na(conv_hull)] <- NA
    b <- b %>% as("SpatRaster")
    
    sp <- 
      supercells(b,
                 step = 5,
                 compactness = 10)
    
    w <-
      queen_weights(sp)
    
    d <- sp[4:8]
    
    
    map(seq(200, 500, length = 5) %>% round(), function(k_fact){    # **************
      
      print(str_glue("   k factor: {k_fact}"))
      
      k <- round(pol$area/k_fact)
      
      if(k < 2){
        
        d_f <- 
          d %>% 
          summarize() %>% 
          select(geometry)
        
      } else {
        
        invisible(capture.output(sp_agg <- skater(k, w, d, cpu_threads = 1)))
        
        d <- 
          d %>% 
          mutate(cluster = sp_agg$Clusters)
        
        d_f <- 
          d %>% 
          group_by(cluster) %>% 
          summarize() %>% 
          select(geometry)
        
      }
      
      return(d_f)
      
    })
    
  })


clus_assess_f <- 
  clus_assess %>% 
  transpose %>% 
  map(function(m){
    
    m %>% 
      bind_rows() %>%
      bind_rows(land_p_small) %>% 
      mutate(id = row_number())
    
  })


# clean up
clus_assess_f <- 
  clus_assess_f %>% 
  map(function(m){
    
    m %>% 
      st_rasterize(block_max_quintiles %>% 
                     slice(wl, 1) %>% 
                     mutate(a = NA) %>% 
                     select(a)) %>% 
      st_as_sf(as_points = F, merge = T)
    
  })


block_max_bl <- block_max_quintiles %>% slice(wl, 1)
# block_max_bl[is.na(land)] <- NA
block_max_bl <- block_max_bl %>% as("SpatRaster")

plan(multicore)

clus_assess_stats <- 
  clus_assess_f %>% 
  future_map(function(m){
    
    inh <- 
      reg_inhomogeneity(m,
                        block_max_bl,
                        sample_size = 0.5)
    
    iso <- 
      reg_isolation(m,
                    block_max_bl,
                    sample_size = 0.25)
    
    list(inh = inh, iso = iso)
    
  })


filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]
  
  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}


clus_assess_stats %>% 
  transpose() %>%
  imap_dfr(function(mm, stat){
    
    mm %>% 
      map_dfr(function(m){
        
        tibble(value = m,
               n_reg = length(m),
               stat = stat)
        
      }) 
  }) %>% 
  mutate(n_reg = factor(n_reg)) %>% 
  
  group_by(stat) %>% 
  mutate(value2 = filter_lims(value)) %>% 
  
  {
    ggplot(., aes(x = n_reg, group = n_reg, y = value2)) +
      geom_boxplot(outlier.shape = NA, na.rm = T) +
      # scale_y_continuous(limits = quantile(.$value, c(0.1, 0.9), na.rm = T)) +
      guides(fill = "none") +
      facet_wrap(~stat, ncol = 1, scales = "free_y")
  }


clus_assess_f[[3]] %>% mapview::mapview()





# *****************************************************************************

land_p <- 
  c(land,
    block_max_quintiles %>% slice(wl, 1)) %>% 
  mutate(a = ifelse(is.na(bp1), NA, a)) %>% 
  
  select(a) %>% 
  st_as_sf(as_points = F, merge = T) %>% 
  
  mutate(id = row_number(),
         area = st_area(.),
         area = area %>% set_units(NULL) %>% {./1000000000})


land_p_small <- 
  land_p %>% 
  filter(area < 700) %>% 
  select(geometry)

land_p_large <- 
  land_p %>% 
  filter(area > 700)

block_max_bl <- block_max_quintiles %>% slice(wl, 1)

regions <- 
  map(seq_len(nrow(land_p_large)), function(i){
    
    print(str_glue("super-region {i} / {nrow(land_p_large)}"))
    
    pol <- land_p_large[i, ]
    
    conv_hull <-
      pol %>% 
      st_boundary() %>% 
      st_polygonize() %>% 
      st_rasterize(block_max_quintiles %>% 
                     slice(wl, 1) %>% 
                     select(bp1) %>% 
                     mutate(bp1 = NA))
    
    b <- block_max_bl
    b[is.na(conv_hull)] <- NA
    b <- b %>% as("SpatRaster")
    
    sp <- 
      supercells(b,
                 step = 5,
                 compactness = 10)
    
    w <-
      queen_weights(sp)
    
    d <- sp[4:8]
    
    k_fact = 350
    
    k <- round(pol$area/k_fact)
    
    if(k < 2){
      
      d_f <- 
        d %>% 
        summarize() %>% 
        select(geometry)
      
    } else {
      
      invisible(capture.output(sp_agg <- skater(k, w, d, cpu_threads = 1)))
      
      d <- 
        d %>% 
        mutate(cluster = sp_agg$Clusters)
      
      d_f <- 
        d %>% 
        group_by(cluster) %>% 
        summarize() %>% 
        select(geometry)
      
    }
    
    return(d_f)
      
  })


regions_f <- 
  regions %>% 
  bind_rows() %>%
  bind_rows(land_p_small) %>% 
  mutate(id = row_number())
  

# clean up
regions_f <- 
  regions_f %>% 
  st_rasterize(block_max_quintiles %>% 
                 slice(wl, 1) %>% 
                 mutate(a = NA) %>% 
                 select(a)) %>% 
  st_as_sf(as_points = F, merge = T)

st_write(regions_f, "regions_0p5.gpkg")


