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

# Calculate observed SD and AR
- Load Data/era5/sst/GOA/EBS_1950_2025.nc
- Calculate the mean SST for each month in the EBS and GOA. Use a weighted mean based on the size of cells (by latitude).
- Restrict to winter (November-March), with year corresponding to January.
- Calculate AR(1) and SD for each time series on a 15-year rolling window basis (right aligned).
- Plot the resulting time series. 

# Debug
- I get the error below at line 98. Explain why and institute a fix.
-Error in MtrxSet(x, dim, type = "POLYGON", needClosed = TRUE) : 
  polygons not (all) closed

  # Expand to North Pacific
  - Rewrite ERA5 download to include North Pacific domain: 20-66N, 110-250E
  - Rewrite time domain to 1950:2024
  - Revise the map on line 29 to show the three polygons 
  - Make the map Pacific-centered

  # Reset rate limit
  - I get the error below when attempting to download the entire ERA5 data set - institute a fix
  - Error: rate limit exceeded--- check your retry rate!rate limit exceeded--- check your retry rate!429--- check your retry rate!Rate limit exceeded. Please wait 3 seconds.--- check your retry rate!https://cds.climate.copernicus.eu/api/retrieve/v1/jobs/b6e76f3b-7a6a-44bb-9beb-6017412b197d--- check your retry rate!2d1f7c6e-1e06-4aab-99bd-00a1f328d4ec--- check your retry rate!

  # Add North Pacific to AR1/SD plots   
  - Add calculations of AR(1) and SD time series for the entire North Pacific domain
  - Add the resulting time series to the plots

  # Use cell-wise detrended monthly anomalies rather than area-weighted mean SST values
  - Calculate monthly climatology for each SST cell in the North Pacific for 1950-1979 
  - Express the time series for each cell as the normalized anomaly from that climatology (units of SD)
  - Detrend the anomaly time series for each cell
  - Plug the resulting data into the workflow for calculating and plotting AR1 and SD time series 
  - Vectorize the lm() loop

  # Debug
  - The order of operations appears to be incorrect.
  - Move from cell-wise calculation of detrended anomalies to an area-wide calculation for each of the three domains (EBS, GOA, N. Pac)
  - Specifically: 
    -- Calculate the area-wide mean SST (weighted by cell area) for each domain
    -- Calculate the montly climatology for each domain for 1950-1979
    -- Calculate the anomaly time series for the entire 1950-2025 period (in units of degrees C)
    -- Detrend the anomaly time series
    -- Use these detrended anomaly time series to calculate 15-year rolling-window AR1 and SD time series
    -- Update the plots with these new time series.

# CCF plots
- Now calculate autocorrelation values for each domain for months Add this to the current end of the script
- Use the detrended monthly SST anomaly time series for all months (January-December) in this step.
- Break the data into three eras: 1950-1988, 1989-2010, 2011-2025
- Plot the autocorrelation values over lags 0-60 months. Each spatial domain gets its own panel, and each era is plotted in a separate color (lines and points) 

# Clean up
- Remove the confidence intervals (I didn't ask for them)
- Use these eras instead of the ones originally specified: 1950-1988, 1989-2000, 2001-2025

# Time series plots
- Add a plot of the full monthly SST anomaly time series (not detrended) for the three areas
- Add plots of annual winter (November-March, with the year corresponding to January) values of non-detrended anomalies. Plot as points and lines.

# SLP analysis
- Create a new script to download monthly mean SLP values from ERA5 for 20-70N, 120-250E, 1950-2025

# Map
- Add a Pacific-centered map with winter (November-March) mean values plotted.

# SLP EOF
- Add to this script to calculate the leading EOF (Empirical Orthogonal Function) of this SLP field.
- EOF should be fit to detrended monthly anomalies.
- I only need EOF 1 and 2 and the corresponding PC time series, so if there is a way calculate only those two leading axes to reduce computing time, use that.
- EOF should be fit to the covariance matrix and the SLP data should be weighted by cell area based on latitude.
- The EOF should include SLP over land and ocean.
- Remove the references to variance explained from the script and the plots.
- Add this box to the EOF 1 plot to evaluate if it identifies the area of highest loadings: x <- c(191, 191, 208, 208, 191), y <- c(44, 55, 55, 44, 44)

# SLP time series for capturing Aleutian Low variability
- Extract the mean detrended monthly SLP anomaly for the area inside the box in the EOF1 panel.
- Calculate winter (November-March) means, with the year corresponding to January.
- Save the resulting time series to the "Output" file.

# Evaluate a driver-response relationship
- Return to CESM_obs_comparison.R
- Goal is to evaluate whether variability in the AL is related to AR1 (reddening) in SST for the GOA and EBS.
- Follow these steps:
    - Load the AL_winter_SLP_anomaly.csv file from Output
    - Calculate a right-aligned SD time series on 15-year rolling windows as was done earlier for SST
    - Make a two-panel scatter plot with AL SLP SD as the x-variable and AR1 for the GOA and EBS as the y-variables.

# Test the relationships
- Fit linear regressions with first-order autocorrelation in the residuals and add the resulting p and R2 values to each panel.

# Now add time series plots
- Make a two-panel plot with the EBS AR(1) time series and the SLP SD time series, both with points and lines; other (first) panel is for GOA SST AR1 and AL SLP SD

# Other windows
- Expand this plot to include windows of 10, 15, 20, and 25 years, for both the SLP SD and the SST AR1 results. 
- Facet with ecosystem in the rows and window length in the columns. 