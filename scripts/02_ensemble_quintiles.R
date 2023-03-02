

# SETUP -----------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(stars)
library(furrr)
library(units)

options(future.fork.enable = T)
plan(multicore)

source("scripts/functions.R")

dir_tmp <- "/mnt/pers_disk"
# dir_ensembled <- "/mnt/bucket_mine/results/global_heat_pf/02_ensembled"

doms <- c("SEA", "AUS", "CAS", "WAS", "EAS", "AFR", "EUR", "NAM", "CAM", "SAM")

wls <- c("0.5", "1.0", "1.5", "2.0", "2.5", "3.0")

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

# load table of variables
# tb_vars <-
#   read_csv("/mnt/bucket_mine/pf_variable_table.csv") %>% 
#   suppressMessages()




# DOMAIN LOOP -----------------------------------------------------------------

for(dom in doms){
  
  print(str_glue(" "))
  print(str_glue("PROCESSING {dom}"))
  
  
  
  # VARIABLE LOOP -------------------------------------------------------------
  
  # for(derived_var in derived_vars){
  #   
  #   print(str_glue(" "))
  #   print(str_glue("Processing {derived_var}"))
    
    
    
    ## IMPORT DERIVED VAR FILES -----------------------------------------------
    
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
    
    
    
    ## SLICE BY WARMING LEVELS ------------------------------------------------
    
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
    
    
    ## CALCULATE STATISTICS ---------------------------------------------------
    
    l_s_wl_stats <-
      
      # loop through warming levels
      imap(wls, function(wl, iwl){
        
        print(str_glue("Calculating stats WL {wl}"))
        
        l_s_wl %>%
          pluck(iwl) %>%
          
          st_apply(c(1,2), function(ts){
            
            # if a given grid cell is empty, propagate NAs
            if(any(is.na(ts))){
              
              # c(perc0 = NA,
              #   perc25 = NA, 
              #   perc50 = NA,
              #   perc75 = NA,
              #   perc100 = NA)
              
              c(bp1 = NA,
                bp2 = NA, 
                bp3 = NA, 
                bp4 = NA, 
                bp5 = NA)
              
              
            } else {
              
              # quantile(ts) %>% 
              #   setNames(c("perc0", "perc25", "perc50", "perc75", "perc100"))
              
              boxplot.stats(ts) %>% 
                .$stats %>% 
                setNames(c("bp1", "bp2", "bp3", "bp4", "bp5"))
              
            }
            
          },
          FUTURE = T,
          .fname = "q") %>%
          aperm(c(2,3,1)) %>%
          split("q")
        
      })
    
    
    # concatenate warming levels
    s_result <-
      l_s_wl_stats %>%
      {do.call(c, c(., along = "wl"))} %>%
      st_set_dimensions(3, values = as.numeric(wls))
    
    
    
    ## SAVE RESULT ------------------------------------------------------------
    
    print(str_glue("Saving result"))
    
    res_filename <- 
      str_glue(
        "{dir_tmp}/ensemble_quint/{dom}_quintiles_ensemble.nc"
      )
    
    fn_write_nc(s_result, res_filename, "wl")
    
    
  # }
  
}

