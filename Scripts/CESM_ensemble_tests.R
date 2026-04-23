# PURPOSE: CESM2 FCM vs MDM ensemble-level attribution test for SST reddening.
#          Fits a GAMM of the form  AR(1) ~ s(year)  (with AR(1) residuals)
#          to each member's 15-yr right-aligned rolling winter SST AR(1) series,
#          and extracts the effective degrees of freedom (EDF) of the smooth.
#          EDF > 1 indicates nonlinear temporal structure in the reddening
#          trajectory; ERA5 observations show a clearly nonlinear pattern,
#          so EDF is a more appropriate descriptor than a single linear slope.
#
#          Compares: per-member EDF distributions (FCM and MDM, GOA and EBS)
#          against ERA5 as a single draw. If ERA5 EDF sits in the tail of MDM
#          but inside the FCM distribution, wind variability is required to
#          reproduce the observed nonlinear reddening pattern.
#
#          GOA and EBS boxes match the ERA5 polygons in analysis.R exactly.
#
# Outputs:
#   Output/CESM_ensemble_sst_ar1.csv             — per-member rolling AR(1) series
#   Output/CESM_ensemble_AR1_EDF.csv             — per-member EDF values
#   Figures/CESM_ensemble_AR1_EDF_distributions.png

source("./Scripts/load.libs.functions.R")

library(ggplot2)
library(dplyr)
library(tidyr)
library(ncdf4)
library(zoo)
library(sf)
library(mgcv)

# ============================================================
# CONFIG
# ============================================================

fcm.sst.dir <- "./Data/CESM2 ensemble/SST/FCM"
mdm.sst.dir <- "./Data/CESM2 ensemble/SST/MDM"

sst.cache <- "./Output/CESM_ensemble_sst_ar1.csv"
edf.cache <- "./Output/CESM_ensemble_AR1_EDF.csv"

win.width <- 15
yr.min    <- 1964   # first complete 15-yr right-aligned window from 1950
yr.max    <- 2014   # CESM2 historical runs end 2014

# GOA / EBS polygons — match analysis.R exactly (0-360 space)
goa.x.360 <- c(191, 191, 203, 205, 208, 223, 234)
goa.y     <- c(50,  53,  57.5, 59,  61,  61,  50)
ebs.x.360 <- c(183, 183, 203, 203, 191)
ebs.y     <- c(53,  65,  65,  57.5, 53)

# Convert to -180/180 for sf masking
goa.x <- ifelse(goa.x.360 > 180, goa.x.360 - 360, goa.x.360)
ebs.x <- ifelse(ebs.x.360 > 180, ebs.x.360 - 360, ebs.x.360)

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

member_id <- function(file) {
  bn <- basename(file)
  m <- regmatches(bn, regexpr("LE2-[0-9]{4}\\.[0-9]{3}", bn))
  if (length(m) == 1 && nzchar(m)) return(m)
  m2 <- regmatches(bn, regexpr("MD\\.[0-9]{2}", bn))
  if (length(m2) == 1 && nzchar(m2)) return(m2)
  tools::file_path_sans_ext(bn)
}

# ============================================================
# SST: per-member winter GOA/EBS rolling AR(1)
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

  goa.coords <- rbind(cbind(goa.x, goa.y), c(goa.x[1], goa.y[1]))
  ebs.coords <- rbind(cbind(ebs.x, ebs.y), c(ebs.x[1], ebs.y[1]))
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

# ============================================================
# ERA5 OBSERVED SINGLE-DRAW COMPARATOR
# ============================================================
# Source analysis.R if winter.roll is not already in the environment.

era.sst.ar1 <- tryCatch({
  if (!exists("winter.roll")) source("./Scripts/analysis.R", local = TRUE, echo = FALSE)
  winter.roll %>% select(year, GOA.ar1, EBS.ar1) %>%
    filter(year >= yr.min, year <= yr.max)
}, error = function(e) {
  warning("ERA5 SST AR(1) unavailable: ", e$message); NULL
})

# ============================================================
# TEST: GAMM(AR(1) ~ s(year)) — extract effective degrees of freedom
# ============================================================
# Rolling AR(1) values from overlapping 15-yr windows are themselves
# strongly autocorrelated, so we fit gamm() with corAR1() residuals.
# The smooth EDF measures how "wiggly" the reddening trajectory is:
# EDF = 1 is linear; higher EDF indicates nonlinear temporal structure.

gamm_edf <- function(df, yvar) {
  df2 <- df %>% filter(!is.na(.data[[yvar]]))
  if (nrow(df2) < 20) return(NA_real_)
  fit <- tryCatch(
    gamm(as.formula(paste0(yvar, " ~ s(year)")),
         correlation = corAR1(form = ~ 1),
         data = df2, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  s.tab <- summary(fit$gam)$s.table
  if (is.null(s.tab) || nrow(s.tab) < 1) return(NA_real_)
  s.tab[1, "edf"]
}

if (file.exists(edf.cache)) {
  message("Loading cached EDF values: ", edf.cache)
  edf.df <- read.csv(edf.cache)
} else {
  message("Fitting GAMM per member (AR(1) ~ s(year))...")
  groups <- sst.df %>% group_by(model, member) %>% group_split()
  edf.df <- lapply(groups, function(d) {
    cat("  ", d$model[1], "-", d$member[1], "\n")
    data.frame(
      model  = d$model[1],
      member = d$member[1],
      GOA    = gamm_edf(d, "GOA.ar1"),
      EBS    = gamm_edf(d, "EBS.ar1")
    )
  }) %>% bind_rows() %>%
    pivot_longer(c(GOA, EBS), names_to = "region", values_to = "edf")
  write.csv(edf.df, edf.cache, row.names = FALSE)
  message("Saved: ", edf.cache)
}

edf.df <- edf.df %>% mutate(region = factor(region, levels = c("GOA", "EBS")))

# ERA5 EDF (single draw)
era.edf <- if (!is.null(era.sst.ar1)) {
  data.frame(
    region = c("GOA", "EBS"),
    edf    = c(gamm_edf(era.sst.ar1, "GOA.ar1"),
               gamm_edf(era.sst.ar1, "EBS.ar1"))
  ) %>% mutate(region = factor(region, levels = c("GOA", "EBS")))
} else NULL

# One-sided empirical p: fraction of members with EDF >= ERA5 EDF
pvals <- if (!is.null(era.edf)) {
  edf.df %>% left_join(era.edf %>% rename(era.edf = edf), by = "region") %>%
    group_by(model, region) %>%
    summarise(n     = sum(!is.na(edf)),
              p.ge  = mean(edf >= era.edf, na.rm = TRUE),
              .groups = "drop")
} else NULL

if (!is.null(pvals)) {
  print(pvals)
  write.csv(pvals, "./Output/CESM_ensemble_AR1_EDF_pvals.csv", row.names = FALSE)
}

# ============================================================
# PLOT: per-member EDF distributions with ERA5 overlay
# ============================================================

p.edf <- ggplot(edf.df, aes(x = edf, fill = model)) +
  geom_histogram(bins = 20, color = "gray30", alpha = 0.7, position = "identity") +
  { if (!is.null(era.edf))
      geom_vline(data = era.edf, aes(xintercept = edf),
                 color = "black", linewidth = 0.9, linetype = "dashed")
    else NULL } +
  facet_grid(region ~ model, scales = "free") +
  scale_fill_manual(values = c("FCM" = "darkorange", "MDM" = "steelblue")) +
  labs(title    = "Effective degrees of freedom for GAMM smooth on year",
       subtitle = paste0("AR(1) ~ s(year) with corAR1() residuals; 15-yr rolling winter AR(1); ",
                         yr.min, "-", yr.max,
                         ". Dashed line = ERA5; EDF = 1 is linear."),
       x = "EDF  (smooth on year)", y = "Number of ensemble members") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/CESM_ensemble_AR1_EDF_distributions.png",
       plot = p.edf, width = 9, height = 6, dpi = 300)
message("Saved: Figures/CESM_ensemble_AR1_EDF_distributions.png")

message("CESM ensemble EDF test complete.")
