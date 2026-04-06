# PURPOSE: Summarize ORAS5 MLD output — cell-wise mean maps
# Author: Mike Litzow
#
# Produces two maps saved to Figures/:
#   1. Annual mean MLD  (all months, all years)
#   2. Winter mean MLD  (November-March, all years)

source("./Scripts/load.libs.functions.R")

# ---- PRECOMPUTE MEANS VIA CDO (avoids loading full 2GB array into R) ----
# Data is on the ORCA025 tripolar grid (~1/4 degree resolution).
# CDO streams the file without loading all time steps at once.

cdo <- "/usr/local/bin/cdo"
src <- "./Data/oras5_mld_NP_1958_2025.nc"

if (!file.exists("./Data/oras5_mld_NP_mean_all.nc")) {
  message("Computing all-month mean via CDO...")
  system(paste(cdo, "timmean", src, "./Data/oras5_mld_NP_mean_all.nc"))
}

if (!file.exists("./Data/oras5_mld_NP_mean_win.nc")) {
  message("Computing winter (Nov-Mar) mean via CDO...")
  system(paste(cdo, "timmean -selmon,11,12,1,2,3", src,
               "./Data/oras5_mld_NP_mean_win.nc"))
}

# ---- LOAD PRECOMPUTED MEANS (small files, single time step each) ----
load_mean_nc <- function(file) {
  nc   <- nc_open(file)
  raw  <- ncvar_get(nc, "somxl030")
  mld  <- if (length(dim(raw)) == 3) raw[, , 1] else raw
  lons <- ncvar_get(nc, "nav_lon")
  lats <- ncvar_get(nc, "nav_lat")
  nc_close(nc)
  mld[mld > 1e10] <- NA
  list(mld = mld, lons = lons, lats = lats)
}

d_all <- load_mean_nc("./Data/oras5_mld_NP_mean_all.nc")
d_win <- load_mean_nc("./Data/oras5_mld_NP_mean_win.nc")

mean_all <- d_all$mld;  lons <- d_all$lons;  lats <- d_all$lats
mean_win <- d_win$mld

# ---- BUILD DATA FRAMES ----
# Flatten 2D grid
df_all <- data.frame(
  lon = as.vector(lons),
  lat = as.vector(lats),
  mld = as.vector(mean_all)
) %>%
  filter(!is.na(mld), mld > 0) %>%
  mutate(lon = ifelse(lon < 0, lon + 360, lon)) %>%
  filter(lat >= 20, lat <= 66, lon >= 110, lon <= 250)

df_win <- data.frame(
  lon = as.vector(lons),
  lat = as.vector(lats),
  mld = as.vector(mean_win)
) %>%
  filter(!is.na(mld), mld > 0) %>%
  mutate(lon = ifelse(lon < 0, lon + 360, lon)) %>%
  filter(lat >= 20, lat <= 66, lon >= 110, lon <= 250)

# ---- MAP SETUP ----
mapWorld <- map_data("world", wrap = c(20, 380))

# Study region polygons in 0-360 longitude space
goa.poly <- data.frame(
  lon = c(198, 198, 203, 205, 208, 225, 231, 201),
  lat = c(54, 55.5, 57.5, 59, 61, 61, 54, 54),
  region = "GOA"
)
ebs.poly <- data.frame(
  lon = c(183, 183, 203, 203, 191),
  lat = c(53, 65, 65, 57.5, 53),
  region = "EBS"
)
polys <- bind_rows(goa.poly, ebs.poly)

mld_map <- function(df, title) {
  ggplot() +
    geom_point(data = df, aes(x = lon, y = lat, color = mld),
               shape = 15, size = 0.4) +
    geom_polygon(data = mapWorld,
                 aes(x = long, y = lat, group = group),
                 fill = "gray80", color = "gray50", linewidth = 0.2) +
    geom_polygon(data = polys,
                 aes(x = lon, y = lat, group = region),
                 fill = NA, color = "black", linewidth = 0.7) +
    scale_color_gradientn(
      colors = oceColorsPalette(64),
      name   = "MLD (m)"
    ) +
    coord_cartesian(xlim = c(110, 250), ylim = c(20, 66)) +
    labs(title = title,
         x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      plot.title      = element_text(hjust = 0.5),
      legend.position = "right"
    )
}

# ---- PLOTS ----
dir.create("./Figures", showWarnings = FALSE)

p_all <- mld_map(df_all,
  "ORAS5 MLD — Cell-wise Mean, All Months (1958-2025)")

p_win <- mld_map(df_win,
  "ORAS5 MLD — Cell-wise Mean, Winter (Nov-Mar, 1958-2025)")

ggsave("./Figures/ORAS5_MLD_mean_all_months.png",
       plot = p_all, width = 10, height = 5, dpi = 150)
ggsave("./Figures/ORAS5_MLD_mean_winter.png",
       plot = p_win, width = 10, height = 5, dpi = 150)

message("Saved: Figures/ORAS5_MLD_mean_all_months.png")
message("Saved: Figures/ORAS5_MLD_mean_winter.png")

# ---- TIME SERIES — GOA AND EBS ----
# Strategy: use CDO sellonlatbox to extract small bounding boxes (manageable in R),
# then apply polygon masking inside R before computing weighted spatial means.
#
# Polygons are defined in 0-360 space; ORAS5 nav_lon is -180/180.
# Bounding boxes converted accordingly:
#   GOA: lon 198-231E = -162 to -129W, lat 54-61N
#   EBS: lon 183-203E = -177 to -157W, lat 53-65N

# CDO bounding box extraction (run once, cached)
goa.box.file <- "./Data/oras5_mld_GOA_box.nc"
ebs.box.file <- "./Data/oras5_mld_EBS_box.nc"

if (!file.exists(goa.box.file)) {
  message("Extracting GOA bounding box via CDO...")
  system(paste(cdo, "sellonlatbox,-162,-129,54,61", src, goa.box.file))
}
if (!file.exists(ebs.box.file)) {
  message("Extracting EBS bounding box via CDO...")
  system(paste(cdo, "sellonlatbox,-177,-157,53,65", src, ebs.box.file))
}

# Load a bounding-box file, return list with mld [x,y,t], nav_lon, nav_lat, dates
load_box <- function(file) {
  nc     <- nc_open(file)
  mld    <- ncvar_get(nc, "somxl030")
  lons   <- ncvar_get(nc, "nav_lon")
  lats   <- ncvar_get(nc, "nav_lat")
  time   <- ncvar_get(nc, "time_counter")
  tunits <- ncatt_get(nc, "time_counter", "units")$value
  nc_close(nc)
  mld[mld > 1e10] <- NA
  origin <- sub(".*since ", "", tunits)
  dates  <- as.Date(as.POSIXct(time, origin = origin, tz = "UTC"))
  list(mld = mld, lons = lons, lats = lats, dates = dates)
}

goa.d <- load_box(goa.box.file)
ebs.d <- load_box(ebs.box.file)

# Build sf polygon objects for masking (coordinates in -180/180)
make_sf_poly <- function(lon360, lat) {
  lon180 <- ifelse(lon360 > 180, lon360 - 360, lon360)
  coords  <- rbind(cbind(lon180, lat), c(lon180[1], lat[1]))
  st_sf(geometry = st_sfc(st_polygon(list(coords)), crs = 4326))
}

goa.sf <- make_sf_poly(
  c(198, 198, 203, 205, 208, 225, 231, 201),
  c(54,  55.5, 57.5, 59,  61,  61,  54,  54)
)
ebs.sf <- make_sf_poly(
  c(183, 183, 203, 203, 191),
  c(53,  65,  65,  57.5, 53)
)

# Compute latitude-weighted spatial mean time series for one region
region_ts <- function(d, poly.sf) {
  lons <- as.vector(d$lons)
  lats <- as.vector(d$lats)

  # Point-in-polygon mask
  grid.sf <- st_as_sf(data.frame(lon = lons, lat = lats),
                      coords = c("lon", "lat"), crs = 4326)
  in.poly <- lengths(st_within(grid.sf, poly.sf)) > 0

  if (sum(in.poly) == 0) stop("No grid cells inside polygon")

  weights <- cos(lats[in.poly] * pi / 180)
  nt      <- dim(d$mld)[3]

  # Flatten spatial dims: mld is [x, y, t] or [xy, t] after matrix()
  mat <- matrix(d$mld, nrow = prod(dim(d$mld)[1:2]), ncol = nt)

  sapply(seq_len(nt), function(t) {
    vals <- mat[in.poly, t]
    if (all(is.na(vals))) return(NA_real_)
    weighted.mean(vals, w = weights, na.rm = TRUE)
  })
}

message("Computing GOA spatial mean time series...")
goa.ts <- region_ts(goa.d, goa.sf)
message("Computing EBS spatial mean time series...")
ebs.ts <- region_ts(ebs.d, ebs.sf)

dates  <- goa.d$dates
years  <- as.integer(format(dates, "%Y"))
months <- as.integer(format(dates, "%m"))

# Build monthly data frame, compute anomaly and detrended anomaly
make_monthly <- function(ts.vec, yrs, mons) {
  df <- data.frame(year = yrs, month = mons, mld = ts.vec)
  clim <- df %>%
    group_by(month) %>%
    summarise(clim = mean(mld, na.rm = TRUE), .groups = "drop")
  df <- df %>%
    left_join(clim, by = "month") %>%
    mutate(anom      = mld - clim,
           anom.detr = residuals(lm(anom ~ seq_along(anom))))
  df
}

goa.monthly <- make_monthly(goa.ts, years, months)
ebs.monthly <- make_monthly(ebs.ts, years, months)

# Winter means (Nov-Mar, year = January year)
make_winter <- function(df, region) {
  df %>%
    filter(month %in% c(11, 12, 1, 2, 3)) %>%
    mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
    group_by(win.year) %>%
    summarise(raw       = mean(mld,       na.rm = TRUE),
              anom      = mean(anom,      na.rm = TRUE),
              anom.detr = mean(anom.detr, na.rm = TRUE),
              .groups   = "drop") %>%
    rename(year = win.year) %>%
    mutate(region = region)
}

winter <- bind_rows(
  make_winter(goa.monthly, "GOA"),
  make_winter(ebs.monthly, "EBS")
)

region.colors <- c("GOA" = "steelblue4", "EBS" = "darkorange4")

ts_long <- winter %>%
  pivot_longer(c(raw, anom, anom.detr),
               names_to = "series", values_to = "value") %>%
  mutate(series = factor(series,
    levels = c("raw", "anom", "anom.detr"),
    labels = c("Raw MLD (m)",
               "Monthly anomaly (m)",
               "Detrended anomaly (m)")))

p_ts <- ggplot(ts_long, aes(x = year, y = value, color = region)) +
  geom_line() +
  geom_point(size = 1.2) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray50") +
  scale_color_manual(values = region.colors) +
  facet_wrap(~ series, ncol = 1, scales = "free_y") +
  labs(title = "ORAS5 MLD — GOA and EBS Winter (Nov-Mar) Time Series",
       x = "Year", y = NULL, color = "Region") +
  theme_bw() +
  theme(plot.title       = element_text(hjust = 0.5),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom")

ggsave("./Figures/ORAS5_MLD_timeseries_winter.png",
       plot = p_ts, width = 10, height = 8, dpi = 150)
message("Saved: Figures/ORAS5_MLD_timeseries_winter.png")

# ---- CELLWISE REGRESSION: SLP PC1 -> MLD ----
# Regress winter-mean detrended MLD anomaly (y) on winter SLP PC1 (x)
# at each grid cell. Restricted to NE Pacific: 170-250E, >=30N.
# Uses GLS with AR(1) residuals (nlme::gls + corAR1).
# Contour line drawn around cells with p < 0.05.
#
# CDO pipeline (run once, cached):
#   1. ymonmean   -> monthly climatology
#   2. ymonsub    -> monthly anomalies
#   3. detrend    -> detrended monthly anomalies
#   4. selmon     -> keep Nov-Mar only

mld.winmon.file <- "./Data/oras5_mld_NP_detr_anom_winmon.nc"

if (!file.exists(mld.winmon.file)) {
  message("CDO: computing detrended anomalies and selecting winter months...")
  system(paste(
    cdo, "selmon,11,12,1,2,3",
    "-detrend",
    "-ymonsub", src,
    paste0("-ymonmean ", src),
    mld.winmon.file
  ))
}

# GLS with AR(1) residuals for each cell — skip everything if cached results exist
reg.cache <- "./Output/oras5_mld_slp_pc1_regression.csv"
if (file.exists(reg.cache)) {
  message("Loading cached regression results...")
  reg.df <- read.csv(reg.cache)
} else {
  # Pre-crop to NE Pacific bounding box via CDO before loading into R
  mld.ne.file <- "./Data/oras5_mld_NEP_detr_anom_winmon.nc"
  if (!file.exists(mld.ne.file)) {
    message("CDO: cropping to NE Pacific bounding box...")
    system(paste(cdo, "sellonlatbox,-190,-110,30,66",
                 mld.winmon.file, mld.ne.file))
  }

  message("Loading cellwise winter-month MLD anomalies (NE Pacific)...")
  nc     <- nc_open(mld.ne.file)
  mld.w  <- ncvar_get(nc, "somxl030")
  lons.w <- ncvar_get(nc, "nav_lon")
  lats.w <- ncvar_get(nc, "nav_lat")
  time.w <- ncvar_get(nc, "time_counter")
  tun.w  <- ncatt_get(nc, "time_counter", "units")$value
  nc_close(nc)

  mld.w[mld.w > 1e10] <- NA

  origin.w   <- sub(".*since ", "", tun.w)
  dates.w    <- as.Date(as.POSIXct(time.w, origin = origin.w, tz = "UTC"))
  years.w    <- as.integer(format(dates.w, "%Y"))
  months.w   <- as.integer(format(dates.w, "%m"))
  winyears.w <- ifelse(months.w %in% c(11, 12), years.w + 1L, years.w)

  lons.w360  <- ifelse(lons.w < 0, lons.w + 360, lons.w)
  ne.mask    <- lons.w360 >= 170 & lons.w360 <= 250 & lats.w >= 30

  pc1        <- read.csv("./Output/NP_PC1_winter.csv")
  common.yrs <- sort(intersect(unique(winyears.w), pc1$year))
  pc1.vec    <- pc1$PC1[match(common.yrs, pc1$year)]
  n.yrs      <- length(common.yrs)

  nx   <- dim(mld.w)[1]; ny <- dim(mld.w)[2]; nt <- dim(mld.w)[3]
  mat  <- matrix(mld.w, nrow = nx * ny, ncol = nt)
  good.cells <- which(ne.mask & rowSums(!is.na(mat)) > nt * 0.5)
  rm(mld.w); gc()

  message("Computing cellwise winter means...")
  win.mat <- matrix(NA_real_, nrow = length(good.cells), ncol = n.yrs)
  for (k in seq_along(common.yrs)) {
    ti <- which(winyears.w == common.yrs[k])
    win.mat[, k] <- rowMeans(mat[good.cells, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(mat); gc()

  message("Fitting GLS AR(1) regressions (this may take several minutes)...")
beta.vec <- rep(NA_real_, length(good.cells))
pval.vec <- rep(NA_real_, length(good.cells))

df.reg <- data.frame(y = NA_real_, x = pc1.vec)

for (k in seq_along(good.cells)) {
  y <- win.mat[k, ]
  if (sum(!is.na(y)) < 10) next
  df.reg$y <- y
  fit <- tryCatch(
    gls(y ~ x, data = df.reg,
        correlation = corAR1(form = ~ 1),
        method = "ML", na.action = na.omit),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  s <- summary(fit)$tTable
  if (nrow(s) < 2) next
  beta.vec[k] <- s[2, 1]
  pval.vec[k] <- s[2, 4]
}

# Assemble full result grid
beta.full <- rep(NA_real_, nx * ny)
pval.full <- rep(NA_real_, nx * ny)
beta.full[good.cells] <- beta.vec
pval.full[good.cells] <- pval.vec

reg.df <- data.frame(
  lon  = as.vector(lons.w360),
  lat  = as.vector(lats.w),
  beta = beta.full,
  pval = pval.full
) %>%
  filter(!is.na(beta)) %>%
  filter(lon >= 170, lon <= 250, lat >= 30)

# Cache regression results to avoid recomputing
  write.csv(reg.df, "./Output/oras5_mld_slp_pc1_regression.csv", row.names = FALSE)
  message("Regression results saved to Output/oras5_mld_slp_pc1_regression.csv")
} # end else (regression computed fresh)

col.lim <- max(abs(reg.df$beta), na.rm = TRUE)

# Significance contour: bin to regular grid then geom_contour
reg.grid <- reg.df %>%
  mutate(lon = round(lon / 0.5) * 0.5,
         lat = round(lat / 0.5) * 0.5) %>%
  group_by(lon, lat) %>%
  summarise(beta = mean(beta, na.rm = TRUE),
            pval = mean(pval, na.rm = TRUE),
            .groups = "drop") %>%
  complete(lon = seq(170, 250, by = 0.5),
           lat = seq(30,  66,  by = 0.5))

p_reg <- ggplot() +
  geom_point(data = filter(reg.df, !is.na(beta)),
             aes(x = lon, y = lat, color = beta),
             shape = 15, size = 0.5) +
  geom_contour(data = filter(reg.grid, !is.na(pval)),
               aes(x = lon, y = lat, z = pval),
               breaks = 0.05, color = "black", linewidth = 0.4) +
  geom_polygon(data = mapWorld,
               aes(x = long, y = lat, group = group),
               fill = "gray80", color = "gray50", linewidth = 0.2) +
  geom_polygon(data = polys,
               aes(x = lon, y = lat, group = region),
               fill = NA, color = "black", linewidth = 0.7) +
  scale_color_gradientn(
    colors   = rev(oceColorsPalette(64)),
    limits   = c(-col.lim, col.lim),
    values   = scales::rescale(c(-col.lim, 0, col.lim)),
    na.value = "white",
    name     = "β (m / PC1)"
  ) +
  coord_cartesian(xlim = c(170, 250), ylim = c(30, 66)) +
  labs(
    title    = "Regression: SLP PC1 → Cellwise Winter MLD Anomaly (detrended)",
    subtitle = paste0("GLS AR(1); contour = p < 0.05; ",
                      min(common.yrs), "–", max(common.yrs)),
    x = "Longitude", y = "Latitude"
  ) +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 8),
        legend.position = "right")

ggsave("./Figures/ORAS5_MLD_SLP_PC1_regression.png",
       plot = p_reg, width = 9, height = 5, dpi = 150)
message("Saved: Figures/ORAS5_MLD_SLP_PC1_regression.png")
