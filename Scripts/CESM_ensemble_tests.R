# PURPOSE: FCM vs MDM CESM2 ensemble-level tests for attributing SST reddening
#          to Aleutian Low variability. These tests do NOT rely on matching
#          individual SLP and SST members to each other — the FCM/MDM contrast
#          (wind variability vs. no wind variability) is the experimental
#          treatment, and the comparison is between ensemble response
#          distributions.
#
# Tests implemented:
#   1. Per-member 15-yr rolling winter SST AR(1) for GOA and EBS.
#      Ribbon = ensemble 5-95% envelope, line = ensemble mean,
#      ERA5 overlaid as single draw.
#   2. Per-member 15-yr rolling winter AL-box SLP SD.
#      Sanity check: FCM should track ERA5 magnitude; MDM should be ~flat.
#   3. Distribution of per-member linear trends in rolling AR(1) (1964-2014).
#      ERA5 trend overlaid as a single draw. Compares whether reddening
#      (positive AR(1) trend) is a forced response, internal variability,
#      or requires AL wind variability.
#
# Output files:
#   Output/CESM_ensemble_sst_ar1.csv
#   Output/CESM_ensemble_al_slp_sd.csv
#   Figures/CESM_ensemble_SST_AR1_trajectories.png
#   Figures/CESM_ensemble_AL_SLP_SD_trajectories.png
#   Figures/CESM_ensemble_AR1_trend_distributions.png

source("./Scripts/load.libs.functions.R")

library(ggplot2)
library(dplyr)
library(tidyr)
library(ncdf4)
library(zoo)
library(sf)
library(patchwork)

# ============================================================
# CONFIG
# ============================================================

fcm.sst.dir <- "./Data/CESM2 ensemble/SST/FCM"
mdm.sst.dir <- "./Data/CESM2 ensemble/SST/MDM"
fcm.slp.dir <- "./Data/CESM2 ensemble/SLP/FCM"
mdm.slp.dir <- "./Data/CESM2 ensemble/SLP/MDM"

sst.cache   <- "./Output/CESM_ensemble_sst_ar1.csv"
slp.cache   <- "./Output/CESM_ensemble_al_slp_sd.csv"

# Overlap period with ERA5 and 15-yr right-aligned windows:
#   CESM2 is 1850-2014; ERA5 baseline starts 1950; first complete 15-yr
#   window right-aligned on 1964 -> window length up to 2014 = 51 values.
win.width <- 15
yr.min    <- 1964
yr.max    <- 2014

# GOA/EBS polygons (matches CESM_AR1_analysis.R, -180/180 space)
goa.x.180 <- ifelse(c(198, 198, 203, 205, 208, 225, 231, 201) > 180,
                    c(198, 198, 203, 205, 208, 225, 231, 201) - 360,
                    c(198, 198, 203, 205, 208, 225, 231, 201))
goa.y     <- c(54, 55.5, 57.5, 59, 61, 61, 54, 54)
ebs.x.180 <- ifelse(c(183, 183, 203, 203, 191) > 180,
                    c(183, 183, 203, 203, 191) - 360,
                    c(183, 183, 203, 203, 191))
ebs.y     <- c(53, 65, 65, 57.5, 53)

# AL box (Litzow et al. 2020 PNAS): 192.5-207.5 E, 45-55 N (0-360 space)
al.lon.min <- 192.5
al.lon.max <- 207.5
al.lat.min <- 45
al.lat.max <- 55

# ============================================================
# HELPERS
# ============================================================

safe_detrend <- function(x) {
  if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
  residuals(lm(x ~ seq_along(x), na.action = na.exclude))
}

safe_ar1 <- function(v) {
  if (any(is.na(v))) return(NA_real_)
  acf(v, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2]
}

roll_ar1 <- function(v, width = win.width) {
  rollapply(v, width = width, fill = NA, align = "right", FUN = safe_ar1)
}

roll_sd <- function(v, width = win.width) {
  rollapply(v, width = width, fill = NA, align = "right", FUN = sd)
}

# Extract member ID from filename (FCM: LE2-####.###; MDM: cesm2.MD.##)
member_id <- function(file) {
  bn <- basename(file)
  m <- regmatches(bn, regexpr("LE2-[0-9]{4}\\.[0-9]{3}", bn))
  if (length(m) == 1 && nzchar(m)) return(m)
  m2 <- regmatches(bn, regexpr("MD\\.[0-9]{2}", bn))
  if (length(m2) == 1 && nzchar(m2)) return(m2)
  tools::file_path_sans_ext(bn)
}

# ============================================================
# SST: per-member winter GOA/EBS AR(1)
# ============================================================

process_sst_member <- function(file) {
  nc   <- nc_open(file)
  sst  <- ncvar_get(nc, "SST", collapse_degen = FALSE)
  lons <- ncvar_get(nc, "lon")
  lats <- ncvar_get(nc, "lat")
  n.t  <- dim(sst)[length(dim(sst))]
  nc_close(nc)

  if (length(dim(sst)) == 4) sst <- sst[, , 1, ]
  sst[sst > 1e20] <- NA_real_

  years.m  <- rep(1850:2014, each = 12)[seq_len(n.t)]
  months.m <- rep(1:12, times = 165)[seq_len(n.t)]

  goa.coords <- rbind(cbind(goa.x.180, goa.y), c(goa.x.180[1], goa.y[1]))
  ebs.coords <- rbind(cbind(ebs.x.180, ebs.y), c(ebs.x.180[1], ebs.y[1]))
  goa.sf <- st_sf(geometry = st_sfc(st_polygon(list(goa.coords)), crs = 4326))
  ebs.sf <- st_sf(geometry = st_sfc(st_polygon(list(ebs.coords)), crs = 4326))

  lons.sf <- ifelse(lons > 180, lons - 360, lons)
  grid    <- expand.grid(lon = lons.sf, lat = lats)
  grid.sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)

  in.goa  <- lengths(st_within(grid.sf, goa.sf)) > 0
  in.ebs  <- lengths(st_within(grid.sf, ebs.sf)) > 0
  weights <- cos(grid$lat * pi / 180)

  wmean <- function(mask) {
    w <- weights[mask]
    apply(sst, 3, function(sl) {
      v <- as.vector(sl)[mask]
      if (all(is.na(v))) NA_real_ else weighted.mean(v, w, na.rm = TRUE)
    })
  }

  monthly <- data.frame(
    year  = years.m,
    month = months.m,
    GOA   = wmean(in.goa),
    EBS   = wmean(in.ebs)
  )

  clim <- monthly %>%
    filter(year >= 1950, year <= 1979) %>%
    group_by(month) %>%
    summarise(GOA.c = mean(GOA, na.rm = TRUE),
              EBS.c = mean(EBS, na.rm = TRUE),
              .groups = "drop")

  monthly <- monthly %>%
    left_join(clim, by = "month") %>%
    mutate(GOA.a = safe_detrend(GOA - GOA.c),
           EBS.a = safe_detrend(EBS - EBS.c))

  winter <- monthly %>%
    filter(month %in% c(11, 12, 1, 2, 3)) %>%
    mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
    group_by(win.year) %>%
    summarise(GOA = mean(GOA.a, na.rm = TRUE),
              EBS = mean(EBS.a, na.rm = TRUE),
              n   = n(), .groups = "drop") %>%
    filter(n == 5) %>%
    rename(year = win.year) %>% select(-n) %>% arrange(year)

  winter %>%
    mutate(GOA.ar1 = roll_ar1(GOA),
           EBS.ar1 = roll_ar1(EBS)) %>%
    filter(year >= yr.min, year <= yr.max) %>%
    select(year, GOA.ar1, EBS.ar1)
}

# ============================================================
# SLP: per-member winter AL-box SLP SD (15-yr rolling)
# ============================================================

process_slp_member <- function(file) {
  nc   <- nc_open(file)
  psl  <- ncvar_get(nc, "PSL", collapse_degen = FALSE)  # Pa; [lon x lat x time]
  lons <- ncvar_get(nc, "lon")
  lats <- ncvar_get(nc, "lat")
  n.t  <- dim(psl)[length(dim(psl))]
  nc_close(nc)

  if (length(dim(psl)) == 4) psl <- psl[, , 1, ]
  psl[psl > 1e20] <- NA_real_
  psl <- psl / 100    # Pa -> hPa

  years.m  <- rep(1850:2014, each = 12)[seq_len(n.t)]
  months.m <- rep(1:12, times = 165)[seq_len(n.t)]

  # AL-box mask (lons are 0-360 in CESM)
  lon.idx <- which(lons >= al.lon.min & lons <= al.lon.max)
  lat.idx <- which(lats >= al.lat.min & lats <= al.lat.max)
  if (length(lon.idx) == 0 || length(lat.idx) == 0)
    stop("AL box empty for ", basename(file))

  w <- cos(lats[lat.idx] * pi / 180)
  w.mat <- matrix(w, nrow = length(lon.idx), ncol = length(lat.idx), byrow = TRUE)

  al.monthly <- apply(psl[lon.idx, lat.idx, , drop = FALSE], 3, function(sl) {
    if (all(is.na(sl))) NA_real_ else weighted.mean(sl, w.mat, na.rm = TRUE)
  })

  df <- data.frame(year = years.m, month = months.m, SLP = al.monthly)

  clim <- df %>% filter(year >= 1950, year <= 1979) %>%
    group_by(month) %>%
    summarise(c = mean(SLP, na.rm = TRUE), .groups = "drop")

  df <- df %>% left_join(clim, by = "month") %>%
    mutate(anom = safe_detrend(SLP - c))

  winter <- df %>%
    filter(month %in% c(11, 12, 1, 2, 3)) %>%
    mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
    group_by(win.year) %>%
    summarise(SLP = mean(anom, na.rm = TRUE), n = n(), .groups = "drop") %>%
    filter(n == 5) %>% rename(year = win.year) %>% select(-n) %>% arrange(year)

  winter %>%
    mutate(AL.sd = roll_sd(SLP),
           AL.sd = as.numeric(scale(AL.sd))) %>%
    filter(year >= yr.min, year <= yr.max) %>%
    select(year, AL.sd)
}

# ============================================================
# RUN (with caching)
# ============================================================

run_ensemble <- function(dir, model.type, fn) {
  files <- list.files(dir, pattern = "\\.nc$", full.names = TRUE)
  cat(model.type, ":", length(files), "files\n")
  lapply(files, function(f) {
    cat("  ", basename(f), "\n")
    out <- tryCatch(fn(f), error = function(e) { cat("    ERROR:", e$message, "\n"); NULL })
    if (!is.null(out)) out %>% mutate(member = member_id(f), model = model.type)
  }) %>% bind_rows()
}

if (file.exists(sst.cache)) {
  message("Loading cached SST AR(1): ", sst.cache)
  sst.df <- read.csv(sst.cache)
} else {
  message("Processing CESM2 SST ensembles...")
  sst.df <- bind_rows(
    run_ensemble(fcm.sst.dir, "FCM", process_sst_member),
    run_ensemble(mdm.sst.dir, "MDM", process_sst_member)
  )
  write.csv(sst.df, sst.cache, row.names = FALSE)
  message("Saved: ", sst.cache)
}

if (file.exists(slp.cache)) {
  message("Loading cached SLP SD: ", slp.cache)
  slp.df <- read.csv(slp.cache)
} else {
  message("Processing CESM2 SLP ensembles...")
  slp.df <- bind_rows(
    run_ensemble(fcm.slp.dir, "FCM", process_slp_member),
    run_ensemble(mdm.slp.dir, "MDM", process_slp_member)
  )
  write.csv(slp.df, slp.cache, row.names = FALSE)
  message("Saved: ", slp.cache)
}

# ============================================================
# ERA5 OBSERVED SINGLE-DRAW COMPARATORS
# ============================================================
# Use analysis.R outputs if already produced; otherwise compute inline here
# would duplicate analysis.R. Rely on caches produced by analysis.R.

era.sst.ar1 <- tryCatch({
  # winter.roll from analysis.R stores GOA.ar1 / EBS.ar1 by year; re-derive if needed
  if (exists("winter.roll")) {
    winter.roll %>% select(year, GOA.ar1, EBS.ar1) %>%
      filter(year >= yr.min, year <= yr.max)
  } else {
    source("./Scripts/analysis.R", local = TRUE, echo = FALSE)
    winter.roll %>% select(year, GOA.ar1, EBS.ar1) %>%
      filter(year >= yr.min, year <= yr.max)
  }
}, error = function(e) { warning("ERA5 SST AR(1) unavailable: ", e$message); NULL })

era.al.sd <- tryCatch({
  if (exists("al")) {
    al %>% select(year, AL.sd) %>% filter(year >= yr.min, year <= yr.max)
  } else {
    NULL
  }
}, error = function(e) NULL)

# ============================================================
# TEST 1: SST AR(1) ENSEMBLE TRAJECTORIES
# ============================================================

sst.long <- sst.df %>%
  pivot_longer(c(GOA.ar1, EBS.ar1), names_to = "region", values_to = "ar1") %>%
  mutate(region = recode(region, "GOA.ar1" = "GOA", "EBS.ar1" = "EBS"))

env.sst <- sst.long %>% group_by(model, region, year) %>%
  summarise(ar1.mean = mean(ar1, na.rm = TRUE),
            lo = quantile(ar1, 0.05, na.rm = TRUE),
            hi = quantile(ar1, 0.95, na.rm = TRUE),
            .groups = "drop")

era.sst.long <- if (!is.null(era.sst.ar1)) {
  era.sst.ar1 %>%
    pivot_longer(c(GOA.ar1, EBS.ar1), names_to = "region", values_to = "ar1") %>%
    mutate(region = recode(region, "GOA.ar1" = "GOA", "EBS.ar1" = "EBS"))
} else NULL

p.sst <- ggplot(env.sst, aes(x = year)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = 0.25) +
  geom_line(aes(y = ar1.mean, color = model), linewidth = 0.9) +
  { if (!is.null(era.sst.long))
      list(geom_line(data = era.sst.long, aes(y = ar1), color = "black", linewidth = 0.7),
           geom_point(data = era.sst.long, aes(y = ar1), color = "black", size = 1.2))
    else NULL } +
  facet_grid(region ~ model) +
  scale_color_manual(values = c("FCM" = "darkorange4", "MDM" = "steelblue4")) +
  scale_fill_manual(values  = c("FCM" = "darkorange",  "MDM" = "steelblue")) +
  labs(title    = "CESM2 ensemble SST AR(1), 15-yr rolling winter means",
       subtitle = "Ribbon = 5-95% ensemble envelope; line = ensemble mean; black = ERA5",
       x = "Year", y = "AR(1)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"))

ggsave("./Figures/CESM_ensemble_SST_AR1_trajectories.png",
       plot = p.sst, width = 9, height = 6, dpi = 300)
message("Saved: Figures/CESM_ensemble_SST_AR1_trajectories.png")

# ============================================================
# TEST 2: AL SLP SD ENSEMBLE TRAJECTORIES (sanity check)
# ============================================================

env.slp <- slp.df %>% group_by(model, year) %>%
  summarise(sd.mean = mean(AL.sd, na.rm = TRUE),
            lo = quantile(AL.sd, 0.05, na.rm = TRUE),
            hi = quantile(AL.sd, 0.95, na.rm = TRUE),
            .groups = "drop")

p.slp <- ggplot(env.slp, aes(x = year)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = 0.25) +
  geom_line(aes(y = sd.mean, color = model), linewidth = 0.9) +
  { if (!is.null(era.al.sd))
      list(geom_line(data = era.al.sd, aes(y = AL.sd), color = "black", linewidth = 0.7),
           geom_point(data = era.al.sd, aes(y = AL.sd), color = "black", size = 1.2))
    else NULL } +
  facet_wrap(~ model, ncol = 2) +
  scale_color_manual(values = c("FCM" = "darkorange4", "MDM" = "steelblue4")) +
  scale_fill_manual(values  = c("FCM" = "darkorange",  "MDM" = "steelblue")) +
  labs(title    = "CESM2 ensemble AL SLP SD, 15-yr rolling winter (z-scored)",
       subtitle = "MDM by design has no interannual AL variability",
       x = "Year", y = "AL SLP SD (z)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"))

ggsave("./Figures/CESM_ensemble_AL_SLP_SD_trajectories.png",
       plot = p.slp, width = 9, height = 4, dpi = 300)
message("Saved: Figures/CESM_ensemble_AL_SLP_SD_trajectories.png")

# ============================================================
# TEST 3: DISTRIBUTION OF PER-MEMBER AR(1) TRENDS
# ============================================================
# For each ensemble member, fit lm(ar1 ~ year) across yr.min:yr.max
# and collect the slope. Gives an internal-variability null distribution
# of reddening trends for each ensemble.

fit_trend <- function(df, yvar) {
  df %>% filter(!is.na(.data[[yvar]])) %>%
    group_by(model, member) %>%
    summarise(slope = {
      if (n() >= 10) coef(lm(.data[[yvar]] ~ year))[2] else NA_real_
    }, .groups = "drop")
}

trend.goa <- fit_trend(sst.df, "GOA.ar1") %>% mutate(region = "GOA")
trend.ebs <- fit_trend(sst.df, "EBS.ar1") %>% mutate(region = "EBS")
trends    <- bind_rows(trend.goa, trend.ebs) %>%
  mutate(region = factor(region, levels = c("GOA", "EBS")))

era.trend <- if (!is.null(era.sst.ar1)) {
  data.frame(
    region = c("GOA", "EBS"),
    slope  = c(
      coef(lm(GOA.ar1 ~ year, data = era.sst.ar1 %>% filter(!is.na(GOA.ar1))))[2],
      coef(lm(EBS.ar1 ~ year, data = era.sst.ar1 %>% filter(!is.na(EBS.ar1))))[2]
    )
  ) %>% mutate(region = factor(region, levels = c("GOA", "EBS")))
} else NULL

# One-sided empirical p: fraction of ensemble members with slope >= ERA5 slope
pvals <- if (!is.null(era.trend)) {
  trends %>% left_join(era.trend %>% rename(era.slope = slope), by = "region") %>%
    group_by(model, region) %>%
    summarise(p = mean(slope >= era.slope, na.rm = TRUE), .groups = "drop")
} else NULL

if (!is.null(pvals)) {
  print(pvals)
  write.csv(pvals, "./Output/CESM_ensemble_AR1_trend_pvals.csv", row.names = FALSE)
}

p.trend <- ggplot(trends, aes(x = slope, fill = model)) +
  geom_histogram(bins = 20, color = "gray30", alpha = 0.7, position = "identity") +
  { if (!is.null(era.trend))
      geom_vline(data = era.trend, aes(xintercept = slope),
                 color = "black", linewidth = 0.9, linetype = "dashed")
    else NULL } +
  facet_grid(region ~ model, scales = "free_y") +
  scale_fill_manual(values = c("FCM" = "darkorange", "MDM" = "steelblue")) +
  labs(title    = "Per-member linear trend in 15-yr rolling SST AR(1) (1964-2014)",
       subtitle = "Dashed line = ERA5 observed trend. One-sided p = fraction of members with slope \u2265 ERA5.",
       x = "AR(1) trend (per year)", y = "Number of ensemble members") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/CESM_ensemble_AR1_trend_distributions.png",
       plot = p.trend, width = 9, height = 6, dpi = 300)
message("Saved: Figures/CESM_ensemble_AR1_trend_distributions.png")

message("All CESM ensemble tests complete.")
