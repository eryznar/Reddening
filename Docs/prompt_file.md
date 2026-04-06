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

# Debug
- Something is wrong with the AL SLP SD time series. Are these supposed to be z-scores? They should be - and thus should be centered on 0.
- Make this fix.

# Scatter plot
- Now make a scatter plot for each SLP SD - SST AR1 combination, with p and R2 values as above.

# Examine residuals over time
- Fit linear regressions to the AL SLP SD - GOA SST AR1 and AL SLP SD - EBS SST AR1 relationships, based on 15-year rolling windows.
- Plot the time series of residuals for each in a two-panel plot. Include a horizontal line at 0.

# Examine CESM results
- Now the analysis moves on to using CESM ensembles run with and without interannual variability in wind stress fields to evaluate the role of AL variability as a driver of reddening in GOA/EBS SST variability.
- This step depends on comparing the Fully Coupled Model (FCM) ensmeble with the Mechanistically Decoupled Model (MDM) ensemble.
- Begin by plotting the SST AR(1) time series for the GOA and EBS based on 15-year windows, following the same workflow as for the ERA5 observations.
- For both FCM and MDM, plot the ensemble mean as a heavy line and the individual ensemble members for the individual ensemble members.
- Result should be a four-panel plot with ecosystem in rows and model type (FCM or MDM) in columns.
- Create the initial script to plot only five ensemble members from each model type to test the script. 


# Debug
- For the first call to process_ensemble on line 129, I get the error below (with traceback). Diagnose, explain, and fix.
Error in sst[, , 1, ] : incorrect number of dimensions
9.
process_member(files[i])
8.
mutate(., member = i, model = model.type)
7.
process_member(files[i]) %>% mutate(member = i, model = model.type)
6.
FUN(X[[i]], ...)
5.
lapply(seq_along(files), function(i) {
cat(" member", i, "\n")
process_member(files[i]) %>% mutate(member = i, model = model.type)
})
4.
list2(...)
3.
bind_rows(.)
2.
lapply(seq_along(files), function(i) {
cat(" member", i, "\n")
process_member(files[i]) %>% mutate(member = i, model = model.type)
}) %>% bind_rows()
1.
process_ensemble(fcm.dir, "FCM", n.members)

# Expand 
- Make the following changes to the script:
    - Plot all ensemble members.
    - Plot the 95% CI as a ribbon around the mean line, rather than plotting the individual ensemble members.
    - Add the observed AR(1) time series (ERA5) to each panel.

# Clarification
- Change the CESM2 plots to winter annual values as for ERA5
- This means using the November-March means, with the year corresponding to January, as was done for CESM_obs_comparison.R
- Also limit this to years of overlap between ERA5 and CESM2.

# Debug
- I now get the error below. Diagnose, explain, and fix.
- Error in `mutate()`:
ℹ In argument: `GOA.anom = residuals(lm(GOA.anom ~ seq_along(GOA.anom)))`.
Caused by error in `lm.fit()`:
! 0 (non-NA) cases
Run `rlang::last_trace()` to see where the error occurred.

# Fix time series length
- The CESM2 time series should only be plotted for right-aligned 15-year windows that can be calculated from the time series available for ERA5 and CESM2 (individual years for 1950-2014).

# Debug again.
- I seem to be getting the same or similar error as described in line 185 of the prompt file. Diagnose, explain, and fix.
- Error in `mutate()`:
ℹ In argument: `GOA.anom = residuals(lm(GOA.anom ~ seq_along(GOA.anom)))`.
Caused by error in `lm.fit()`:
! 0 (non-NA) cases
Run `rlang::last_trace()` to see where the error occurred.

> rlang::last_trace()
<error/dplyr:::mutate_error>
Error in `mutate()`:
ℹ In argument: `GOA.anom = residuals(lm(GOA.anom ~ seq_along(GOA.anom)))`.
Caused by error in `lm.fit()`:
! 0 (non-NA) cases
---
Backtrace:
     ▆
  1. ├─global process_ensemble(fcm.dir, "FCM", n.members)
  2. │ ├─... %>% bind_rows()
  3. │ └─base::lapply(...)
  4. │   └─FUN(X[[i]], ...)
  5. │     ├─process_member(files[i]) %>% mutate(member = i, model = model.type)
  6. │     └─global process_member(files[i])
  7. │       └─monthly %>% left_join(clim, by = "month") %>% ...
  8. ├─dplyr::bind_rows(.)
  9. │ └─rlang::list2(...)
 10. ├─dplyr::mutate(., member = i, model = model.type)
 11. ├─dplyr::mutate(...)
 12. ├─dplyr:::mutate.data.frame(...)
 13. │ └─dplyr:::mutate_cols(.data, dplyr_quosures(...), by)
 14. │   ├─base::withCallingHandlers(...)
 15. │   └─dplyr:::mutate_col(dots[[i]], data, mask, new_columns)
 16. │     └─mask$eval_all_mutate(quo)
 17. │       └─dplyr (local) eval()
 18. ├─stats::residuals(lm(GOA.anom ~ seq_along(GOA.anom)))
 19. └─stats::lm(GOA.anom ~ seq_along(GOA.anom))
 20.   └─stats::lm.fit(...)
 21.     └─base::stop("0 (non-NA) cases")
Run rlang::last_trace(drop = FALSE) to see 3 hidden frames.

# Debug
- Still getting this error. Fix.
-Error in `mutate()`:
ℹ In argument: `GOA.ar1 = rollapply(...)`.
Caused by error in `na.fail.default()`:
! missing values in object
Run `rlang::last_trace()` to see where the error occurred.

> rlang::last_trace()
<error/dplyr:::mutate_error>
Error in `mutate()`:
ℹ In argument: `GOA.ar1 = rollapply(...)`.
Caused by error in `na.fail.default()`:
! missing values in object
---
Backtrace:
     ▆
  1. ├─global process_ensemble(fcm.dir, "FCM", n.members)
  2. │ ├─... %>% bind_rows()
  3. │ └─base::lapply(...)
  4. │   └─FUN(X[[i]], ...)
  5. │     ├─process_member(files[i]) %>% mutate(member = i, model = model.type)
  6. │     └─global process_member(files[i])
  7. │       └─... %>% select(year, GOA.ar1, EBS.ar1)
  8. ├─dplyr::bind_rows(.)
  9. │ └─rlang::list2(...)
 10. ├─dplyr::mutate(., member = i, model = model.type)
 11. ├─dplyr::select(., year, GOA.ar1, EBS.ar1)
 12. ├─dplyr::filter(., year >= 1964, year <= 2014)
 13. ├─dplyr::mutate(...)
 14. ├─dplyr:::mutate.data.frame(...)
 15. │ └─dplyr:::mutate_cols(.data, dplyr_quosures(...), by)
 16. │   ├─base::withCallingHandlers(...)
 17. │   └─dplyr:::mutate_col(dots[[i]], data, mask, new_columns)
 18. │     └─mask$eval_all_mutate(quo)
 19. │       └─dplyr (local) eval()
 20. ├─zoo::rollapply(...)
 21. └─zoo:::rollapply.default(...)
 22.   ├─zoo::coredata(rollapply(zoo(data), ...))
 23.   ├─zoo::rollapply(zoo(data), ...)
 24.   └─zoo:::rollapply.zoo(zoo(data), ...)
 25.     └─base::mapply(...)
 26.       └─zoo (local) `<fn>`(dots[[1L]][[15L]], dots[[2L]][[15L]], data = `<dbl>`)
 27.         └─FUN(data[posns], ...)
 28.           └─stats::acf(v, lag.max = 1, plot = FALSE)
 29.             ├─stats::na.action(as.ts(x))
 30.             └─stats:::na.fail.default(as.ts(x))
 31.               └─base::stop("missing values in object")
Run rlang::last_trace(drop = FALSE) to see 3 hidden frames.
> 
# SLP-SST correlation in observations
- Evaluate the direct relationship 
- Plot cross-correlations between detrended SLP anomalies in the AL box with detrended SST anomalies for the two systems (annual Nov-March mean for SLP and SST).
- Plot CCF for unsmoothed SLP/SST data, 2-year rolling means and three-year rolling means.

- Add a similar analysis using November-March values of N. Pac PC1 instead of AL SLP.

# Examine MLD responses
- Start a new script to download mixed layer depth fields for the N. Pacific spatial domain for 1950-2025 from ORAS5.
- Use the N. Pacific spatial domain that was queried for SST.

# Debug
- I get this error. Fix.
- Error in copernicus_marine_subset(dataset_id = "global-reanalysis-phy-001-031-grepv2-monthly",  : 
  could not find function "copernicus_marine_subset"

  # Restart
  - I am in a new session and trying to run ORAS5_MLD_download.R
  - I get this error when running the script:
  sh: copernicusmarine: command not found
Warning message:
In system(cmd) : error in running command
- I have copernicusmarine open in the terminal for this VS Code session
- But I am running the script in an R studio session.
- How can I successfully run the R script calling the MLD data from Copernicus?

# Debug
- I get this warning and this error. Diagnose, explain, and fix.
- WARNING - 2026-04-02T23:27:45Z - '--force-download' has been deprecated.
ERROR - 2026-04-02T23:27:49Z - Dataset not found: global-reanalysis-phy-001-031-grepv2-monthly Please check that the dataset exists and the input datasetID is correct.

# Summarize/confirm MLD output
- Plot maps with cell-wise mean MLD for the North Pacific spatial domain, all years.
- Create separate maps for all months and winter (November-March) only.
- Add the polygons for the EBS and GOA areas to the plots.

# Debug
- At line 16 I get the error below.
- Tell me the spatial resolution of the MLD data. Is it possible to coarsen the data too 1/12 degree to make the vector smaller?
- Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()


# Observation - CESM2 comparison
- Return to the comparison between ERA5 SST and SLP vs CESM2 FCM/MDM
- Create a new R script to do this.
- Rely on data handling and analysis approaches used earlier in this project.
- For both SST and SLP, compare the loadings (on maps) for EOF1.
- Compare the EOF loadings between observations (1 panel) and 11 ensemble members (12 panel plot, 4 columns x 3 rows)
- Four separate plots - SST observations vs FCM, SST observations vs MDM, SLP observations vs FCM, SLP observations vs MDM. 
- Confirm that the EOFs are fit to the same time domain for observations and CESM2
- Rewrite to call needed packages 

# Debug
- I get this error when running the SST EOF for ERA5: cannote allocate vector size of 38.3 Gb
- Confirm that you are using package irlba and fitting only the leading EOF axis
- Confirm that the cells are weighted by size based on latitude, and the EOF is fit on the covariance matrix.
- Confirm the EOF is fitted to detrended monthly anomaly values for each cell (1950-1979 climatology).
- The EOF should be fit to all months, not winter means.