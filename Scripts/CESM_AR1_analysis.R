# PURPOSE: Compare FCM and MDM CESM2 ensemble SST AR(1) for GOA and EBS
# Author: Mike Litzow

source("./Scripts/load.libs.functions.R")

# SETTINGS ----------------------------------

n.members <- NULL   # NULL = all members

fcm.dir <- "./Data/CESM2 ensemble/SST/FCM"
mdm.dir <- "./Data/CESM2 ensemble/SST/MDM"

# GOA and EBS polygon vertices (same as CESM_obs_comparison.R, -180/180)
goa.x <- c(198, 198, 203, 205, 208, 225, 231, 201)
goa.x <- ifelse(goa.x > 180, goa.x - 360, goa.x)
goa.y <- c(54, 55.5, 57.5, 59, 61, 61, 54, 54)

ebs.x <- c(183, 183, 203, 203, 191)
ebs.x <- ifelse(ebs.x > 180, ebs.x - 360, ebs.x)
ebs.y <- c(53, 65, 65, 57.5, 53)

# FUNCTION: process one ensemble member file -> winter AR1 time series ----------------------------------

process_member <- function(file) {

  nc   <- nc_open(file)
  sst  <- ncvar_get(nc, "SST", collapse_degen = FALSE)  # [lon x lat x z_t x time]
  lons <- ncvar_get(nc, "lon")
  lats <- ncvar_get(nc, "lat")
  time <- ncvar_get(nc, "time")         # days since 0000-01-01
  nc_close(nc)

  # Drop z_t dimension (size 1) -> [lon x lat x time]
  sst <- sst[, , 1, ]

  # Replace fill values with NA
  sst[sst > 1e20] <- NA_real_

  # Generate year/month sequence directly from file name date range (1850-01 to 2014-12)
  # More reliable than decoding the 365-day calendar arithmetic
  n.time <- length(time)
  years  <- rep(1850:2014, each = 12)[seq_len(n.time)]
  months <- rep(1:12, times = 165)[seq_len(n.time)]

  # Build sf masks once per unique lon/lat grid
  goa.coords <- rbind(cbind(goa.x, goa.y), c(goa.x[1], goa.y[1]))
  ebs.coords <- rbind(cbind(ebs.x, ebs.y), c(ebs.x[1], ebs.y[1]))
  goa.sf <- st_sf(geometry = st_sfc(st_polygon(list(goa.coords)), crs = 4326))
  ebs.sf <- st_sf(geometry = st_sfc(st_polygon(list(ebs.coords)), crs = 4326))

  # CESM lon is 0-360; convert to -180/180 for sf
  lons.sf <- ifelse(lons > 180, lons - 360, lons)
  grid    <- expand.grid(lon = lons.sf, lat = lats)
  grid.sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)

  in.goa <- lengths(st_within(grid.sf, goa.sf)) > 0
  in.ebs <- lengths(st_within(grid.sf, ebs.sf)) > 0
  weights <- cos(grid$lat * pi / 180)

  wmean <- function(mask) {
    w <- weights[mask]
    apply(sst, 3, function(sl) {
      v <- as.vector(sl)[mask]
      if (all(is.na(v))) NA_real_ else weighted.mean(v, w, na.rm = TRUE)
    })
  }

  monthly <- data.frame(
    year  = years,
    month = months,
    GOA   = wmean(in.goa),
    EBS   = wmean(in.ebs)
  )

  # Monthly anomaly (1950-1979 baseline) then detrend
  clim <- monthly %>%
    filter(year >= 1950, year <= 1979) %>%
    group_by(month) %>%
    summarise(GOA.clim = mean(GOA, na.rm = TRUE),
              EBS.clim = mean(EBS, na.rm = TRUE),
              .groups = "drop")

  safe_detrend <- function(x) {
    if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
    residuals(lm(x ~ seq_along(x), na.action = na.exclude))
  }

  monthly <- monthly %>%
    left_join(clim, by = "month") %>%
    mutate(GOA.anom = GOA - GOA.clim,
           EBS.anom = EBS - EBS.clim,
           GOA.anom = safe_detrend(GOA.anom),
           EBS.anom = safe_detrend(EBS.anom))

  # Winter (Nov-Mar) means; year corresponds to January (Nov/Dec shift to next year)
  winter <- monthly %>%
    filter(month %in% c(11, 12, 1, 2, 3)) %>%
    mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
    group_by(win.year) %>%
    summarise(GOA = mean(GOA.anom, na.rm = TRUE),
              EBS = mean(EBS.anom, na.rm = TRUE),
              n   = n(),
              .groups = "drop") %>%
    filter(n == 5) %>%   # keep only complete winters (all 5 months present)
    rename(year = win.year) %>%
    select(-n) %>%
    arrange(year)

  # 15-year rolling AR(1); restrict to 1964-2014 (first complete right-aligned
  # window from 1950 ends at 1964; CESM data ends at 2014)
  safe_ar1 <- function(v) {
    if (any(is.na(v))) return(NA_real_)
    acf(v, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2]
  }

  winter %>%
    mutate(
      GOA.ar1 = rollapply(GOA, width = 15, fill = NA, align = "right", FUN = safe_ar1),
      EBS.ar1 = rollapply(EBS, width = 15, fill = NA, align = "right", FUN = safe_ar1)
    ) %>%
    filter(year >= 1964, year <= 2014) %>%
    select(year, GOA.ar1, EBS.ar1)
}

# LOAD AND PROCESS ENSEMBLE FILES ----------------------------------

get_files <- function(dir, n) {
  files <- list.files(dir, pattern = "\\.nc$", full.names = TRUE)
  if (!is.null(n)) files <- files[seq_len(min(n, length(files)))]
  files
}

process_ensemble <- function(dir, model.type, n) {
  files <- get_files(dir, n)
  cat("Processing", model.type, "-", length(files), "members\n")
  lapply(seq_along(files), function(i) {
    cat("  member", i, "\n")
    process_member(files[i]) %>% mutate(member = i, model = model.type)
  }) %>% bind_rows()
}

fcm.ar1 <- process_ensemble(fcm.dir, "FCM", n.members)
mdm.ar1 <- process_ensemble(mdm.dir, "MDM", n.members)

all.ar1 <- bind_rows(fcm.ar1, mdm.ar1)

# Ensemble mean and 95% CI across members
ens.stats <- all.ar1 %>%
  pivot_longer(c(GOA.ar1, EBS.ar1), names_to = "region", values_to = "ar1") %>%
  mutate(region = sub(".ar1", "", region, fixed = TRUE)) %>%
  group_by(model, region, year) %>%
  summarise(mean.ar1 = mean(ar1, na.rm = TRUE),
            lo95     = quantile(ar1, 0.025, na.rm = TRUE),
            hi95     = quantile(ar1, 0.975, na.rm = TRUE),
            .groups  = "drop")

# Observed ERA5 AR(1) — from winter data frame computed in CESM_obs_comparison.R
# Load if not already in environment
if (!exists("winter")) {
  stop("Run CESM_obs_comparison.R first to create the 'winter' object, or source it.")
}

obs.long <- winter %>%
  select(year, GOA.ar1, EBS.ar1) %>%
  pivot_longer(-year, names_to = "region", values_to = "ar1") %>%
  mutate(region = sub(".ar1", "", region, fixed = TRUE))

# Restrict all layers to 1964-2014: right-aligned 15-yr windows from 1950 base,
# capped at CESM end year
overlap.years <- 1964:2014

# PLOT ----------------------------------

cesm.plot <- ggplot() +
  geom_ribbon(data = ens.stats %>% filter(!is.na(mean.ar1), year %in% overlap.years),
              aes(x = year, ymin = lo95, ymax = hi95, fill = model),
              alpha = 0.25) +
  geom_line(data = ens.stats %>% filter(!is.na(mean.ar1), year %in% overlap.years),
            aes(x = year, y = mean.ar1, color = model),
            linewidth = 1.0) +
  geom_line(data = obs.long %>% filter(!is.na(ar1), year %in% overlap.years),
            aes(x = year, y = ar1),
            color = "black", linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("FCM" = "steelblue4", "MDM" = "darkorange4")) +
  scale_fill_manual(values  = c("FCM" = "steelblue",  "MDM" = "darkorange")) +
  facet_grid(region ~ model) +
  labs(title = "CESM2 SST AR(1) — 15-year rolling window (winter)",
       subtitle = "Ribbon = 95% CI across members; dashed = ERA5 observed",
       x = "Year", y = "AR(1)",
       color = "Model", fill = "Model") +
  theme_bw() +
  theme(legend.position = "bottom")

print(cesm.plot)
