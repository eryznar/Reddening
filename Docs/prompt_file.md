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
- This step depends on comparing the Fully Coupled Model (FCM) ensemble with the Mechanistically Decoupled Model (MDM) ensemble.
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

# MLD time series
- Continue in the same script.
- Calculate monthly area-wide means for the EBS and GOA. Use area-weighted values for each cell, with the areas based on latitude.
- Calculate the monthly anomalies for each area for the 1958-1987 reference period.
- Summarize three time series as winter (November-March) means with the year corresponding to January:
    - Raw mean MLD.
    - Mean monthly anomaly MLD.
    - Detrended monthly anomaly MLD.
    - Remember that each of these three time series is plotted as a winter mean.
    - Remember that these three time series should be plotted separately for the GOA and EBS areas only.

# Debug
- At line 16 I get the error below.
- Tell me the spatial resolution of the MLD data. Is it possible to coarsen the data too 1/12 degree to make the vector smaller?
- Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()

# Regress SLP PC1 onto cellwise MLD
- Continue adding onto the same script
- Use SLP PC1 for detrended monthly anomaly SLP values as calculated earlier. 
- Plot the regression coefficients for the regression of November-March SLP PC1 (x) onto cellwise detrended MLD monthly anomalies (November-March means, y).
- Restrict the spatial domain to the NE Pacific: 170-250 degrees E, >=30 degrees N.
- Use a regression with 1st-order autocorrelation in the residuals.
- Surround areas with p < 0.05 with a line. 
- Confirm that the SLP PC1 time series from ERA5_SLP_download.R is being used in the regressions (and is then correctly subset fo November-March mean values.)

# Redraw 15-yr windows AL SD vs SST AR1
- Return to CESM_obs_comparison.R
- At line 460, add in script that was previously in place to plot AL SD and EBS/GOA SST AR(1) for 15-year windows.
- AL SD and SST AR1 should be on two labeled y-axes, one panel for EBS, 1 for GOA.
- Save the figure to /figures as a .png.

# Draw 15-yr windows AL SD vs MLD AR1
- Use the same layout and approach as ./Figures/AL_SD_SST_AR1_15yr_dual_axis.png
- Specifically, 15-year windows, annual means of winter (November-March) values.
- Instead of SST AR(1), use MLD AR(1) for each system.-
- This is not a regressopn, rather each time series is plotted, with dual y axes, as in the reference figure.
- Add the correlation coefficient to each panel.

# Draw regression maps for AL - MLD by era
- Make a second version of the map in "Regress SLP PC1 onto cellwise MLD"
- This version will be two panels (one column, two rows)
- First panel is for the low-variability AL era (years 1998-2010).
- Second panel is for the high-variability era (years 2016-2025).
- All other details of the regression and plotting foloow the first map for all years.
- Save as a png in the Figures folder. 

# Compare era differences in cellwise SST and MLD AR(1) values in two maps.
- Think hard. There are multiple steps for this analysis, and it is critical to get them all right.
- In addition to executing the steps presented here, also develop a series of checks that can evalaute whether each step is coded correctly. Do this in terms of:
    - 1) the coding, 
    - 2) consistency with earlier accepted results on this project,
    - 3) the consistency with understanding of physical SLP-SST-MLD relationships. 
- Present this error-checking plan before writing code, and then integrate the approved code checks in the script.
- The goal here is to compare changes in AR(1) values for SST and MLD for era 2 (2003-2024) vs. era 1 (1987-2008).
- (As an aside, these are the years that are included in 15-year windows for 2017-2024 and 2017-2024.)
- Analysis will use annual winter (November-March) means values of SST and MLD for each cell.
- Use the full North Pacific domain (20-66N, 110-250E).
- For each cell, calculate AR(1) values for each cell in each era, then plot the era difference (era2 - era1).
- Make a separate map for SST AR(1) differences (left panel) and for MLD AR(1) differences (right panel).

# Check eras
- I realized the two eras overlap, which is not what we want.
- Plot the time series of mean winter SST for the EBS and GOA, save as a .png file.

# Revisions / next steps
- Return to the era difference maps for SST and MLD.
- Instead of 1987-2008 and 2003-2024 for era 1, and era 2, use 1989-2006 and 2007-2024. 
- This refers to the section of "ORAS5_MLD_summary.R" beginning on line 695.

# Second version, using average MLD
- Create a second version of this plot.
- The only difference is that instead of plotting the change in AR(1) values for MLD, plot the change in average MLD values between eras.

# Now add SST era difference in AR(1) / SLP PC1 - MLD regression maps
- Left-hand panel should be the era difference in SST AR(1) as just implemented.
- Right-hand panel should be the cellwise SLP PC1 - MLD regression map.
- The regression map should be for the same spatial domain as the SST AR(1) map, but otherwise should follow the same workflow as the earlier standalone SLP PC1 - MLD regression map.

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

# SLP PC1 instead of SLP box values
- Revise the p.dual plot on line 503 of the CESM_obs_comparison.R file.
- Instead of the 15-year rolling window Aleutian Low SLP SD values, subsitute 15-year rolling window (right aligned) values of SD in SLP PC1 values.
- Use SLP PC1 values calculated elsewhere in the project. 

# Data update to include winter 2026.
- Combine the existing download scripts for MLD, SLP, and SST into a single script called 'downloads.R'.
- Add a section of code to run an additional download of data from January-March 2026.
- Set up the additional download script so it will not require downloading data that are already downloaded for the project.
- Also set up the additional download so that the new data can be added to the existing data seamlessly.
- Check the full script carefully to ensure it will run as intended, including both the previous sections for downloads, as well as the section for the new months. 

# Clean up workflow into a single analysis script
- Start a new script called "analysis.R"
- Start the script by producing a reference map of the North Pacific that includes the following areas indicated as boxes on the map and referenced in a legend:
    - Area of SLP EOF analysis.
    - EBS and GOA study areas.
    - Area over which Aleutian Low SLP SD is calculated.
    - Save the map to /Figures as a .png

# SST monthly anomaly time series
- To the same analysis.R file, copy over the analysis from CESM_obs_comparison.R to: 
    -Calculate monthly anomalies with respect to the 1950-1979 climatology for EBS, GOA, and North Pacific (20-66N, 110-250E (area: N, W, S, E in -180/180))
    - Calculate annual winter (November-March) anomaly values with the year corresponding to January.
    - Plot the non-detrended time series for 1950-2026 for all three areas (seaparate panels) and save as a .png 

# All months anomaly
- Add an intermediate step, before plotting winter means, to plot the time series of non-detrended anomalies for all months for all three areas.

# AR(1) and SD time series
- Now find the earlier script for calculating EBS and GOA AR(1) and SD values for winter SST on 15-year, right aligned rolling windows and append it as the next step in the script.
- Make the plot of AR(1) and SD for each system (four panels) and save as a .png.
- Make each panel in the plot roughly 4:3 aspect ratio.

# AL SD 
- Now make the plot with EOF1 loadings for SLP (calculations done the same as earlier).
- Plot the loadings on an Albers projection.
- If needed, change the sign so the box holds positive values.
- Plot on a Mercator projection.

# SLP SD vs SST AR1
- Now add the analysis plotting 15-year windows for SLP SD and SST AR(1), one panel per system.

# Significance step
- Use a randomization approach to test the significance of the SLP Sd - SST AR(1) calculations.
- For each system, generate 1000 pairs of simulated SLP SD and SST AR(1) time series, using arima.sim to mimic the actual SD and AR(1) values in the time series plots.
- Be careful and think hard! We are not referencing the SLP and SST SD/AR in the simulations. Instead we are just using the SD and AR1 of the plotted lines and points on each panel!
- For each of the simulated pairs, calculate the correlation coefficient.
- Estimate the p value for each panel as the proportion of simulated correlation coefficients that are as strong as or stronger than the actual correlation coefficient.
- Plot the resulting p-values on each panel.
- Update: in line 474, the code identifies the absolute value of r. But we want the proportion of raw r values, not absolute value, that 1are less than the observed correlation. In other words, the test is one-tailed. 

# Intermediate plot
- Add a plot of all the r values from the randomization for each ecosystem (EBS and GOA) compared to the observed r value in each system.
- The distribution of randomized r values should be a barplot, and the observed r should be indicated by a vertical dashed line.
- Make in a two-panel plot (one for each system), and save as a .png.

# Make SLP box consistent with Litzow et al. (2020) PNAS
- Redefine the SLP box in analysis.R as 45-55N, 192.5-207.5E.

# Evaluate annual SST AR(1) & SD, and relationship with annual SLP SD
- A previous round of analysis on this project used annual SLP / SST data, not winter values.
- Create a new script called annual_analysis_comparison.R for this evaluation.
- Re-write all of the code in analysis.R that uses winter means amd use annual means instead.
- Give any output or figure files that are produced a name with 'annual' appended to make the different outputs (winter/annual) easy to distinguish.

# SLP SD - MLD map
- Return to the analysis.R workflow.
- RESUME USING WINTER-ONLY VALUES FOR ALL ANALYSIS (THIS IS CRITICAL).
- As employed earlier, calculate regression coefficients for winter MLD (response variable) regressed on winter SLP SD for the designated box (explanatory variable). 
- Expand the earlier analysis in space to cover the entire area of SLP values used in the EOF1 calculation.
- Also continue to use GLS regression to control for autocorrelated residuals.
- In addition, add Benjamini-Hochberg False Discovery Rate (FDR) control as described in Wilks (2016, BAMS).
- Add lines around p <= 0.05 after using these controls, as in the original version.
- Append this step to analysis.R and save the resulting figure as a .png file.
- Clarification: AL SLP SD should be annual Nov-Mar mean (sing year winter means) rather than a 15-year window of any kind. Response should be annual Nov-Mar mean MLD. Both variables should be aligned with the year corresponding to January. And the analysis should continue for the entire time series, through winter 2026.

# Revisiting CESM analysis
- Start with three tests on the FCM / MDM ensembles (not linked SST-SLP outputs from individual members)
    Test 1 — SST AR(1) trajectories: per-member 15-yr rolling winter AR(1) in GOA and EBS, then plots the 5–95% envelope + ensemble mean for FCM and MDM, with ERA5 overlaid as a single black line. Saved to Figures/CESM_ensemble_SST_AR1_trajectories.png.

    Test 2 — AL SLP SD sanity check: per-member 15-yr rolling SD of winter AL-box SLP for both ensembles. FCM should track ERA5; MDM should be ~flat by design. Saved to Figures/CESM_ensemble_AL_SLP_SD_trajectories.png.

    Test 3 — Trend distribution test: per-member slope of rolling AR(1) vs. year over 1964–2014, histogrammed by ensemble × region. ERA5's slope is overlaid as a dashed line. Writes an empirical one-sided p-value (fraction of members with slope ≥ ERA5) to Output/CESM_ensemble_AR1_trend_pvals.csv. The key attribution statement: if ERA5's reddening sits in the tail of MDM but inside FCM, wind variability is required.
- Add these in a new CESM_ensemble_tests.R script.

# Revise CESM ensemble analysis
- Test 1 and Test 2 are repeats or expansions on earlier analysis, and do not seem to track what we expect - I'm not sure how different members are initialized and how closely MDM and FCM runs should track each other and/or observations.
- Drop Test 1 and Test 2 from the script.
- Test 3 is executed inappropriately - it evaluates linear trends in SST AR(1), but that is not the pattern that we see in observations - we see a temporal pattern that would be better evaluated by a high-EDF GAM.
- Revise Test 3 to test for the EDF in a GAMM following the form SST AR(1) ~ s(year), where SST AR(1) refers to the time series of values for right-aligned 15-year rolling window.
- Replace the comparison of ensemble member and observation linear trends with the effective degrees of freedom (EDF) value for the smooth on year. 
- Ensure that the same EBS / GOA boxes are being used for ERA5, FCM, and MDM.
- Plot in the same style as the previous linear trend plot, but use a free scale on the y axis; save as a png with a new name referring to EDF, not trend. 
# Seems to be no difference in FCM/MDM yet again!

# CESM individual member analysis
- First, confirm that SST and SLP runs from the same individual members can be identified from the file name information on hand. THIS IS A CRITICAL STEP, AS GETTING THIS WRONG WILL PROVIDE SPURIOUS RESULTS FROM THE DOWNSTREAM ANALYSIS. Think hard, evaluate multiple authoritative sources, and suggest tests for evaluating that the SST/SLP information can be linked to shared member IDs.
- Start a new script called "CESM_member_regression.

# Proceed with correlation
- T1-3 all pass, confirming we have matched SLP-SST pairs from individual members.
- Now institute correlation tests for all FCM and all MDM model runs and compare the distribution of r values for each.
- The correlation is for AL SLP SD vs SST AR(1) on right-aligned 15-year windows for both the EBS and GOA, exactly as was done for the observations.
- Plot the r distributions for each system, FCM vs MDM in a four-panel plot and save as a .png.
- Free scale the y axis.
# CESM2 FCM-MDM comparison appears to be a dead end, return to observations.

# Return to observations
- Resume work on analysis.R
- Add 2-axis plots comparing AL SLP SD and EBS/GOA MLD AR(1) in the same style previously used for SLP SD - SST AR(1).
- Again use winter values on 15-year right-aligned rolling windows.

# Compare AL SLP SD with MLD anomaly 
- Make a plot in the same style and with the same data handling, but replace MLD AR(1) with MLD anomaly.
- Again, 15-year rolling windows, right-aligned, with winter values. 
- MLD anomalies should be winter means of detrended monthly anomalies.

# Regression map - MLD AR(1) on AL SLP SD
- Make a map in the same style as the SLP - MLD regression map.
- The data for these cellwise regressions should be 15-year rolling window SD of detrended AL SLP anomalies for annual winter values (x variable) and cell values of 15-year rolling window AR(1) for detrended MLD anomalies. 

# Correlation map / spatially aggregated regression map
- That pattern of regression coefficients is not spatially coherent, though there seems to be some indication of the NE Pacific - central N Pacific PDO dipole.
- Make two new maps to investigate further.
- First map is the same regression workflow, but with MLD aggregate into 1 degree x 1 degree blocks at the first step. That is, MLD values are averaged over 1 degree blocks, then anomalies are calculated, and the rest of the regression workflow is followed. Save this as a new .png plot with a clarifying name.
- Second map is a for AL SLP SD - cellwise MLD AR(1) correlation coeffients rather than regression coefficients. Otherwise this is the same as the initial regression map. Save as a .png. 
- Include the Modified Chelton adjustment for the significance test in the correlation map, in addition to the FDR control.

# Revisions
- For the 1 degree spatially aggregated regression map, change q for FDR to 0.1 
- Try an additional plot with MLD aggregated at 2 degree blocks.
# Summary - there appears to be a weak PDO-like dipole in the MLD AR(1) - SLP Sd regression, but never significant at the cell level and spatially patchier than the PDO or SST regression on SLP SD

# Era maps
- Read the previous scripts to find an era map in either MLD anomaly or MLD AR(1). Tell me the exact quantities being mapped and the eras used.
- Bring the script to create three maps into analysis.R.
- Use Era 1: 1989-2006 and Era 2: 2007-2024
- Make cellwise maps of delta SST AR(1), delta MLD AR(1), and delta MLD mean.
- Save as three panel plot with era difference for the three variables.
- Each panel should have its own scale bar.
- Indicate the EBS and GOA polygons with thin black lines, no need to identify with legend. 

- Increase map area to the full N. Pacific domain. 
- Add the EBS and GOA legends as thin black lines ON TOP OF the plotted values!
- Make the western edge of the plots at 160E to block out large-magnitude areas in the western Pacific.

# SST AR(1) regressed on AL SLP SD
- Make a map of cellwise SST AR(1) (y variable) regressed on AL SLP SD (x variable).
- Use 15 year rolling windows of winter mean values.
- Only include the area east of 160E.
- Do not add significance.
- Save as .png.
- Add significance lines to the plot as for the MLD - SLP regression map; use the p-values from GLS with AR(1) residuals and FDR.

# Seasonal sensitivity test
- This test is motivated by Newman et al. 2016, Journal of Climate.
- Add a new section to analysis.R.
- Redraw the dual-axis plots of AL SLP SD and SST AR(1) for the EBS and GOA, but instead of using November-March means for each, use November-January for the AL SLP SD and February-April for SST AR(1).
- This is to capture the lag of ocean response to AL forcing documented by Newman et al.
- Everything else in the analysis should stay the same - 15-year rolling windows, year corresponding to January, etc.
- Save the new plot with a distinguishing name as a .png.

# Confirm
- The results appear to be highly sensitive to seasonality.
- In particular, the pattern of SLP SD changes dramatically when calculated for NDJ.
- Double-check this new code is correct. Design and implement independent checks, and report back.

## Note: seasonal sensitivity verified (2026-04-28)
- Independent checks added in Scripts/analysis.R Section 14 (Checks 1-7) confirm code is correct.
- Annual mean SLP NDJ vs NDJFM: r = 0.80 (n = 77). Matches i.i.d. theoretical bound √(3/5) ≈ 0.775; not a bug.
- 15-yr rolling SD correlations (post-scale, full overlap):
    - NDJFM ~ NDJ: r = 0.012 (raw), r = 0.053 (after dropping partial 1950 year). Robust to edge effects.
    - NDJFM ~ FMA: r = 0.698. Canonical NDJFM AL-SD signal is driven by late winter (FMA).
    - NDJ ~ FMA: r = -0.568. Early- and late-winter AL volatility are anti-correlated decadally.
- Implication: the NDJ AL SD series captures a mode essentially orthogonal to the NDJFM signal used elsewhere. The Newman-lag dual-axis plot (NDJ AL SD vs FMA SST AR(1)) tests independent information, not a minor seasonal tweak.
- Diagnostic figures: Figures/AL_SD_NDJ_vs_NDJFM_check.png, Figures/AL_SD_three_seasons_check.png

# Additional seasonal sensitivity check
- This follows the finding of strongest decadal signal for AL SD in FMA, and Newman et al. (2016) finding of ~ 1-3 mo correlations with PDO ocean variability.
- Fit an additional dual-axis plot using JFM AL SLP SD and FMA for SST AR(1).
- Save the resulting plot as a .png with a distinguishing name.

# Map for additional CESM2 FCM / MDM model comparison
- Revisit the cellwise regression of SST AR(1) on AL SLP SD (using observations).
- Replicate the same analysis for coupled SLP-SST outputs for FCM and MDM.
- Continue to use GLS with AR(1) residuals.
- Run the cellwise regression for the full time domain for each ensemble member.
- Use 15-year rolling means for NDJFM values. 
- Exclude cases with incomplete data (years corresponding to January for which we do not have a full set of NDJFM observations).
- Plot the maps of mean cellwise regression values for MDM and FCM.
- Retain the ensemble SD for each cell to allow significance testing later. 
- Final plot should be three-panel, comparing observartions with FCM with MDM.
- Do not attempt to match temporal domain (years available) for observations and CESM2; instead just use the full range of years available for each, while following the directions about excluding incomplete observations (winters) above.
- Save as a .png.

# Revisions
- Put the FCM and MDM results on the same scale bar, different from the observation scale bar.
- This is because the relationships are much stronger for observations than for CESM2.

# FCM vs. MDM significance test
- Make a new plot that shows the FCM and MDM ensemble means, and a third panel that shows the difference with areas of p <= 0.05 outlined.
- FCM/MDM ensemble mean plots should be the same as in previous plot (no observations for comparison).
- Use FDR, but no correction for autotocorrelation - just simple statistical comparison based on the mean, n, and SD at each cell.
- Save as a .png.

# Comparison by system
- Compare mean regression coefficients (area-weighted) for the GOA and EBS polygon in each ensemble member.
- Conduct a t-test comparison of the area means between FCM and MDM.

# Plot revision
- Find the section of analysis.R that plots the regression of winter MLD AR(1) on AL SLP SD (15-yr rolling windows).
- Replot with 160E as the western boundary - all three of the plots at different spatial resolution (cell size).
- Add this as a revision to the relevant section. Do not change the names of the saved figures.

# MLD anomaly - AL SLP SD regression
- Add a new section after the MLD AR(1) - AL SLP SD regression map sections.
- Fit the exact same regression maps, but replace MLD AR(1) values with MLD mean anomaly values.
- All other data handling is the same - e.g., 15-year rolling windows of November-March values, GLS with AR(1) residuals, three different levels of spatial aggregation, etc.
- Save the plots as .png files.
- Renumber the subsequent sections accordingly. 

# MLD AR(1) regression on AL SLP SD compared with the PDO spatial pattern.
- Create a new section in analysis.R to compare these spatial patterns.
- Use the regression map from section 12 - original spatial resolution for ERA5.
- Also calculate and map the leading mode for winter (Nov-Mar) detrended monthly SST anomalies, using the full time domain that was used in section 12 (but no rolling windows for the EOF analysis).
- Make the EOF loadings map in the same format as the regression map.
- Suggest two-three robust statistical methods for evaluating spatial coherence between the two patterns.

# N Pacific SST plot
- Go back to the section that plots monthly anomaly SST time series for EBS, GOA, and N. Pacific and make a version that only plot N. Pacific.
- Save as a .png
- Add a dashed vertical line at the midpoint of 2014.

# Effect of variable AR(1) values on extreme event risk
- Add another section to analysis.R.
- Replot the EBS and GOA annual winter mean monthly SST anomaly time series with lowess trends fit. 
- Plot and save as a .png. 
- Lowess is over-fit. Try a GAM with k = 5.
- That shrinks to a linear fit. Change the script to fit a linear trend and save the trend for further analysis.
- Now calculate and save the SD for each time series.
- Now simulate each time series using the observed trend and SD. Vary the simulations under three conditions: AR(1) = minimum AR(1) value calculated for each system, AR(1) = 0, and AR(1) = the maximum value observed for each time series. Run 1000 simulations under each AR(1) condition for each system.
- Also calculate the climatology for 1950-1979 values (mean and SD).
- For each system under each AR(1) condition, calculate the % of the 1000 simulations > 2 SD above 1950-1979 climatological mean (heatwaves), and the % more than 2 SD below the climatological mean (cold spells).
- Plot the results in a four-panel plot - 2 systems, with the time series of % heatwaves and % cold spells in different panels. Color-code the three AR(1) conditions.

# SST variability
- Revisit the SST SSD time series analysis in Fig. 3 of the results summary. 
- Recalculate as rolling-window winter CV (SD/mean) and plot as a two-panel figure (one panel per ecosystem).
- Append to the analysis.R file at the correct place, save figure.

# MLD variabiity
- Calculate the first and second EOF modes and corresponding PCs for MLD in the area of the N. Pacific basin used elsewhere in analysis. 
- Include winter (NDJFM) only, calculate EOFs using detrended cellwise anomalies, and correct for cell area based on latitude.
- Add this section to a logical place in analysis.R, at the start of the MLD analysis.
- Plot the North et al. (1982) error bars for the first 5 eigenvalues (only calculate the first five modes).
- Plot the loadings and corresponding PC time series for modes 1-2 in a four-panel figure.

# SST regressed on AL SLP anomalies.
- Find the section of analysis that produces the map in MLD_AL_SLP_regression.png.
- This is a regression of MLD anomalies onto SLP anomalies in the AL box.
- Follow the same workflow to produce a map of SST anomaly regression onto AL SLP anomalies, using winter values.
- Plot as a .png.
- Add this to analysis on an adjoining section next to the section producing the original MLD regression map.
- Combine the two panels into a single figure - SST on top, MLD on bottom.
- Plot only the NE Pacific (East of 150E) for both panels.