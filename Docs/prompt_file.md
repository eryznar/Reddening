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

# Wavelets
- Add a wavelet plot for each spatial domain SST anomaly time series for the full time period
- Plot all three domain wavelet plots as panels on a single plot

# Debug
- I get the error below. Explain and fix.
- Error in data.frame(time = dec.year, coi = wt$coi.1, domain = dom) : 
  arguments imply differing number of rows: 912, 916, 1

# Debug again
- I now get the error below.
- Also, the colors are different, from the previous plots, the log no longer shows as base 2 in the figure legend, and the scale for the log values has changed. 
- Diagnose, explain, and fix.
- Warning messages:
1: The following aesthetics were dropped during statistical transformation: fill.
ℹ This can happen when ggplot fails to infer the correct grouping structure in the data.
ℹ Did you forget to specify a `group` aesthetic or to convert a numerical variable into a factor? 
2: The following aesthetics were dropped during statistical transformation: fill.
ℹ This can happen when ggplot fails to infer the correct grouping structure in the data.
ℹ Did you forget to specify a `group` aesthetic or to convert a numerical variable into a factor? 
3: The following aesthetics were dropped during statistical transformation: fill.
ℹ This can happen when ggplot fails to infer the correct grouping structure in the data.
ℹ Did you forget to specify a `group` aesthetic or to convert a numerical variable into a factor? 

# Debug
- Why were the log values positive in the original version (separate plots for each area) but negative now?
- Why is there the cone of uncertainty in the original version, but only a dashed asymmetrical line in the new version?
- Diagnose, explain,and fix

# Query
- For issue 1, revert to original definition of power that gave vales between ~0 and ~13.5 for log2 power 

# Wavelet exploration
Please revise my R wavelet plotting code so the wavelet power spectrum uses a positive-only power scale like a standard Torrence-style wavelet power plot, rather than plotting log2(power/variance) with negative values.
Current problem:
My code divides wt$Power by series.var
then plots fill = log2(power)
this produces negative values whenever local power is less than the series variance
I want a figure like my reference example, where the plotted color scale is positive-only
What I want you to do:
Diagnose this explicitly in the code comments.
Remove the variance normalization used only for plotting.
Plot a positive-only measure of wavelet power.
Prefer plotting raw wt$Power as the fill variable unless you think a different positive-only transform is better.
Keep the COI masking and significance contour logic intact.
Remove interpolate = TRUE from geom_raster() unless there is a strong reason to keep it.
Update legend labels and title text so they match the revised quantity being plotted.
Return the full corrected R code, not just a diff.
Briefly explain the scientific meaning of the revised color scale and how it differs from log2(power/variance).

-Why is the cone of influence not symmetrical? Diagnose and fix the plot.

- Revert to relative scale