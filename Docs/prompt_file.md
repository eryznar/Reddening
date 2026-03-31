# SST reddening
- You are restarting the project in https://github.com/eryznar/Reddening.git
- The goal is to understand trends in SST autocorrelation in the Gulf of Alaska (GOA) and Eastern Bering Sea (EBS)
- Default choice is to work in R, although Matlab and Python can be used in cases where they are more applicable to analysis of climate data sets and atmosphere-ocean interactions

# Compare observed and modeled AR1 values for SST
- Create a file in the Scripts folder called "CESM_obs_comparison.R"
- Begin by plotting the following polygons on a map of the North Pacific to confirm these polygons for defining the two areas:
    # GOA: 
    - goa.x <- c(201, 201, 205, 208, 225, 231, 201)
    - goa.x <- ifelse(goa.x > 180, goa.x-360, goa.x)
    - goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)
    # EBS
    - ebs.x <- c(183, 183, 203, 203, 191)
    - ebs.x <- ifelse(ebs.x > 180, ebs.x-360, ebs.x)
    -ebs.y <- c(53, 65, 65, 57.5, 53)

- Produce the map for confirmation

# Add map   
- Add the coastline to the map.

# Query ERA5 data
- Look at the code at the bottom of the script after "library(ecmwfr)"
- Re-write this to download the necessary spatial domain, all years, all months.
- Save the netcdf file to "Data"

- Add a snippet to use wf_set_key() to use the API key.