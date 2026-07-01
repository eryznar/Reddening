# PURPOSE: Reddening project — main analysis script
# Author: Mike Litzow
#
# Sections:
#   1. Reference map — North Pacific study areas

library(ggplot2)
library(dplyr)
library(maps)
library(ncdf4)
library(zoo)
library(tidyr)
library(sf)
library(patchwork)

# ============================================================
# SECTION 1: REFERENCE MAP — North Pacific study areas
# ============================================================
# Polygons in 0-360 longitude space:
#   GOA polygon : irregular (from CESM_obs_comparison.R)
#   EBS polygon : irregular (from CESM_obs_comparison.R)

mapWorld <- map_data("world", wrap = c(20, 380))

# --- Build polygon data frames (all in 0-360 lon space) ---

# GOA: vertices follow the irregular polygon used for masking
goa.poly <- data.frame(
  lon  = c(191, 191, 203, 205, 208, 223, 234),
  lat  = c(50,  53,  57.5, 59,  61,  61,  50),
  area = "GOA"
)

# EBS: vertices follow the irregular polygon used for masking
ebs.poly <- data.frame(
  lon  = c(183, 183, 203, 203, 191),
  lat  = c(53,  65,  65,  57.5, 53),
  area = "EBS"
)

all.areas <- bind_rows(goa.poly, ebs.poly) %>%
  mutate(area = factor(area, levels = c("GOA", "EBS")))

# --- Axis label helper: convert 0-360 lon to ±180 notation ---
lon_label <- function(x) {
  ifelse(x > 180, paste0(360 - x, "\u00b0W"),
  ifelse(x == 180, "180\u00b0",
         paste0(x, "\u00b0E")))
}

area.colors    <- c("GOA" = "darkorange4",
                    "EBS" = "steelblue4")

area.linetypes <- c("GOA" = "solid",
                    "EBS" = "solid")

area.linewidths <- c("GOA" = 0.9,
                     "EBS" = 0.9)

area.labels <- c(
  "GOA" = "Gulf of Alaska (GOA)",
  "EBS" = "Eastern Bering Sea (EBS)"
)

region.fills <- c("GOA" = "darkorange", "EBS" = "steelblue")

ref.map <- ggplot() +
  geom_polygon(data = mapWorld,
               aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray55", linewidth = 0.2) +
  geom_polygon(data = filter(all.areas, area %in% c("GOA", "EBS")),
               aes(x = lon, y = lat, group = area, fill = area),
               color = NA, alpha = 0.25) +
  scale_fill_manual(values = region.fills, guide = "none") +
  geom_polygon(data = all.areas,
               aes(x = lon, y = lat, group = area,
                   color = area, linetype = area, linewidth = area),
               fill = NA) +
  scale_color_manual(name = NULL, values = area.colors, labels = area.labels) +
  scale_linetype_manual(name = NULL, values = area.linetypes, labels = area.labels) +
  scale_linewidth_manual(name = NULL, values = area.linewidths, labels = area.labels) +
  coord_map(
    projection = "rectangular",
    parameters = 55,
    xlim = c(120, 250),
    ylim = c(20, 66)
  ) +
  scale_x_continuous(
    breaks = seq(120, 250, 30),
    labels = lon_label
  ) +
  scale_y_continuous(
    breaks = seq(20, 60, 10),
    labels = function(y) paste0(y, "\u00b0N")
  ) +
  labs(title = "North Pacific Study Regions", x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.direction = "vertical",
    legend.text      = element_text(size = 9),
    legend.key.width = unit(1.8, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3)
  ) +
  guides(color    = guide_legend(override.aes = list(linewidth = 1)),
         linetype = guide_legend(override.aes = list(linewidth = 1)),
         linewidth = "none")

ggsave("./Figures/NP_reference_map.png",
       plot = ref.map, width = 8, height = 6, dpi = 300)
message("Saved: Figures/NP_reference_map.png")

# ============================================================
# SECTION 2: ERA5 SST WINTER ANOMALY TIME SERIES
# ============================================================
# Calculates monthly SST anomalies (1950-1979 climatology) for GOA, EBS,
# and full North Pacific (20-66N, 110-250E). Summarises as annual winter
# (Nov-Mar) means and plots non-detrended time series for 1950-2026.

# --- Load ERA5 SST (use 2026 file if available, otherwise 2025) ---
sst.file <- if (file.exists("./Data/era5_sst_NP_1950_2026.nc")) {
  "./Data/era5_sst_NP_1950_2026.nc"
} else {
  "./Data/era5_sst_NP_1950_2025.nc"
}
message("Loading SST from: ", sst.file)

nc   <- nc_open(sst.file)
sst  <- ncvar_get(nc, "sst")
lons <- ncvar_get(nc, "longitude")
lats <- ncvar_get(nc, "latitude")
time <- ncvar_get(nc, "valid_time")
nc_close(nc)

dates  <- as.Date(as.POSIXct(time, origin = "1970-01-01", tz = "UTC"))
months <- as.integer(format(dates, "%m"))
years  <- as.integer(format(dates, "%Y"))
sst    <- sst - 273.15   # K → °C

# --- Study-area polygons (same vertices as CESM_obs_comparison.R) ---
goa.x <- c(191, 191, 203, 205, 208, 223, 234)
goa.y <- c(50,  53,  57.5, 59,  61,  61,  50)
ebs.x <- c(183, 183, 203, 203, 191)
ebs.y <- c(53,  65,  65,  57.5, 53)

# Convert 0-360 → -180/180 for sf
to180 <- function(x) ifelse(x > 180, x - 360, x)

make_sf_poly <- function(lon360, lat) {
  coords <- rbind(cbind(to180(lon360), lat), c(to180(lon360[1]), lat[1]))
  st_sf(geometry = st_sfc(st_polygon(list(coords)), crs = 4326))
}

goa.sf <- make_sf_poly(goa.x, goa.y)
ebs.sf <- make_sf_poly(ebs.x, ebs.y)

# ERA5 lons are already in -180/180 space
grid    <- expand.grid(lon = lons, lat = lats)
grid.sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)

in.goa <- lengths(st_within(grid.sf, goa.sf)) > 0
in.ebs <- lengths(st_within(grid.sf, ebs.sf)) > 0
in.np  <- !is.na(as.vector(sst[,,1]))   # all non-land ocean cells

weights <- cos(grid$lat * pi / 180)

weighted_region_mean <- function(mask) {
  w <- weights[mask]
  apply(sst, 3, function(slice) {
    vals <- as.vector(slice)[mask]
    if (all(is.na(vals))) return(NA_real_)
    weighted.mean(vals, w, na.rm = TRUE)
  })
}

message("Computing area-weighted SST means...")
goa.sst <- weighted_region_mean(in.goa)
ebs.sst <- weighted_region_mean(in.ebs)
np.sst  <- weighted_region_mean(in.np)

monthly.sst <- data.frame(year = years, month = months,
                           GOA = goa.sst, EBS = ebs.sst, NP = np.sst)

# --- Monthly climatology (1950-1979) ---
clim.sst <- monthly.sst %>%
  filter(year >= 1950, year <= 1979) %>%
  group_by(month) %>%
  summarise(GOA.clim = mean(GOA, na.rm = TRUE),
            EBS.clim = mean(EBS, na.rm = TRUE),
            NP.clim  = mean(NP,  na.rm = TRUE),
            .groups = "drop")

# --- Anomalies (non-detrended, stored for plotting) ---
monthly.anom <- monthly.sst %>%
  left_join(clim.sst, by = "month") %>%
  mutate(GOA.anom = GOA - GOA.clim,
         EBS.anom = EBS - EBS.clim,
         NP.anom  = NP  - NP.clim,
         date     = as.Date(paste(year, month, "15", sep = "-"))) %>%
  select(date, year, month, GOA.anom, EBS.anom, NP.anom)

region.colors <- c("GOA" = "darkorange4", "EBS" = "steelblue4", "NP" = "forestgreen")
region.labels <- c("GOA" = "Gulf of Alaska", "EBS" = "Eastern Bering Sea",
                   "NP"  = "North Pacific")

# --- All-months anomaly time series plot ---
end.mo <- max(monthly.anom$date, na.rm = TRUE)

p.monthly.anom <- monthly.anom %>%
  pivot_longer(cols = c(GOA.anom, EBS.anom, NP.anom),
               names_to = "region", values_to = "anom") %>%
  mutate(region = recode(region,
                         "GOA.anom" = "GOA",
                         "EBS.anom" = "EBS",
                         "NP.anom"  = "NP"),
         region = factor(region, levels = c("GOA", "EBS", "NP"))) %>%
  ggplot(aes(x = date, y = anom, color = region)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_line(linewidth = 0.35) +
  scale_color_manual(values = region.colors, labels = region.labels) +
  facet_wrap(~ region, ncol = 1, scales = "free_y",
             labeller = labeller(region = region.labels)) +
  labs(title   = paste0("Monthly SST Anomaly, 1950\u2013", format(end.mo, "%Y")),
       x = "Year", y = "SST Anomaly (\u00b0C)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"))

ggsave("./Figures/ERA5_SST_monthly_anom_ts.png",
       plot = p.monthly.anom, width = 7, height = 8, dpi = 300)
message("Saved: Figures/ERA5_SST_monthly_anom_ts.png")

# --- N. Pacific only version (single panel) ---
p.monthly.anom.np <- monthly.anom %>%
  ggplot(aes(x = date, y = NP.anom)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_line(linewidth = 0.4, color = region.colors[["NP"]]) +
  labs(title = paste0("N. Pacific Monthly SST Anomaly, 1950\u2013", format(end.mo, "%Y")),
       x = "Year", y = "SST Anomaly (\u00b0C)") +
  theme_bw(base_size = 11)

ggsave("./Figures/ERA5_SST_monthly_anom_NP_only.png",
       plot = p.monthly.anom.np, width = 8, height = 3.5, dpi = 300)
message("Saved: Figures/ERA5_SST_monthly_anom_NP_only.png")

# --- Annual winter (Nov-Mar) means, year = January year ---
winter.anom <- monthly.anom %>%
  filter(month %in% c(11, 12, 1, 2, 3)) %>%
  mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
  group_by(win.year) %>%
  summarise(GOA = mean(GOA.anom, na.rm = TRUE),
            EBS = mean(EBS.anom, na.rm = TRUE),
            NP  = mean(NP.anom,  na.rm = TRUE),
            .groups = "drop") %>%
  rename(year = win.year) %>%
  arrange(year)

# --- Plot ---
end.yr <- max(winter.anom$year, na.rm = TRUE)

p.win.anom <- winter.anom %>%
  pivot_longer(-year, names_to = "region", values_to = "anom") %>%
  mutate(region = factor(region, levels = c("GOA", "EBS", "NP"))) %>%
  ggplot(aes(x = year, y = anom, color = region)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.4) +
  scale_color_manual(values = region.colors, labels = region.labels) +
  facet_wrap(~ region, ncol = 1, scales = "free_y",
             labeller = labeller(region = region.labels)) +
  scale_x_continuous(breaks = seq(1950, end.yr, 10)) +
  labs(title   = paste0("Winter (Nov\u2013Mar) SST Anomaly, 1950\u2013", end.yr),
       x = "Year", y = "SST Anomaly (\u00b0C)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"))

ggsave("./Figures/ERA5_SST_winter_anom_ts.png",
       plot = p.win.anom, width = 7, height = 8, dpi = 300)
message("Saved: Figures/ERA5_SST_winter_anom_ts.png")

# ============================================================
# SECTION 3: 15-YEAR ROLLING AR(1) AND SD — GOA AND EBS
# ============================================================
# Detrend the winter anomaly series for each region, then compute
# right-aligned 15-year rolling AR(1) and SD. Four-panel plot: AR(1)
# and SD for GOA (top row) and EBS (bottom row).

detrend <- function(x) residuals(lm(x ~ seq_along(x)))

roll_ar1 <- function(x, width = 15) {
  rollapply(x, width = width, fill = NA, align = "right",
            FUN = function(v) acf(v, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2])
}
roll_sd <- function(x, width = 15) {
  rollapply(x, width = width, fill = NA, align = "right", FUN = sd)
}

winter.detr <- winter.anom %>%
  mutate(GOA = detrend(GOA),
         EBS = detrend(EBS),
         NP  = detrend(NP))

winter.roll <- winter.detr %>%
  mutate(GOA.ar1 = roll_ar1(GOA),
         EBS.ar1 = roll_ar1(EBS),
         GOA.sd  = roll_sd(GOA),
         EBS.sd  = roll_sd(EBS))

# Four panels: GOA AR1, EBS AR1, GOA SD, EBS SD
p.goa.ar1 <- ggplot(winter.roll, aes(x = year, y = GOA.ar1)) +
  geom_line(color = "darkorange4") +
  geom_point(color = "darkorange4", size = 1.4) +
  labs(title = "Gulf of Alaska", x = NULL, y = "AR(1)") +
  theme_bw(base_size = 11)

p.ebs.ar1 <- ggplot(winter.roll, aes(x = year, y = EBS.ar1)) +
  geom_line(color = "steelblue4") +
  geom_point(color = "steelblue4", size = 1.4) +
  labs(title = "Eastern Bering Sea", x = NULL, y = "AR(1)") +
  theme_bw(base_size = 11)

p.goa.sd <- ggplot(winter.roll, aes(x = year, y = GOA.sd)) +
  geom_line(color = "darkorange4") +
  geom_point(color = "darkorange4", size = 1.4) +
  labs(x = "Year", y = "SD (\u00b0C)") +
  theme_bw(base_size = 11)

p.ebs.sd <- ggplot(winter.roll, aes(x = year, y = EBS.sd)) +
  geom_line(color = "steelblue4") +
  geom_point(color = "steelblue4", size = 1.4) +
  labs(x = "Year", y = "SD (\u00b0C)") +
  theme_bw(base_size = 11)

p.roll <- (p.goa.ar1 | p.ebs.ar1) / (p.goa.sd | p.ebs.sd) +
  plot_annotation(
    title   = "15-yr Rolling Winter SST AR(1) and SD (right-aligned)",
    theme   = theme(plot.title = element_text(size = 12, hjust = 0.5))
  )

ggsave("./Figures/ERA5_SST_AR1_SD_15yr_rolling.png",
       plot = p.roll, width = 10, height = 6.5, dpi = 300)
message("Saved: Figures/ERA5_SST_AR1_SD_15yr_rolling.png")

# ============================================================
# SECTION 3B: 15-YEAR ROLLING WINTER SST CV (SD / mean)
# ============================================================
# Revisits the Section 3 SD series as a coefficient of variation.
# CV requires a meaningful (non-zero) mean, so it is computed on the
# RAW winter-mean SST (deg C) -- NOT the detrended anomaly used in
# Section 3 (whose mean is ~0 and would make SD/mean undefined).
# Per right-aligned 15-yr window: mean and SD of the raw winter SST,
# CV = SD / mean. Two-panel plot, one ecosystem per panel.
#
# NOTE: SST in deg C is an interval (not ratio) scale, so the absolute
# CV magnitude depends on the deg-C zero point; read each panel as
# relative change through time within a system, not as a cross-system
# absolute ratio. The EBS winter mean can run cold (near 0 deg C), which
# can inflate/destabilise its CV -- inspect before over-interpreting.

# Raw winter (Nov-Mar, year = Jan) mean SST from monthly.sst (deg C)
winter.raw <- monthly.sst %>%
  filter(month %in% c(11, 12, 1, 2, 3)) %>%
  mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
  group_by(win.year) %>%
  summarise(GOA = mean(GOA, na.rm = TRUE),
            EBS = mean(EBS, na.rm = TRUE),
            .groups = "drop") %>%
  rename(year = win.year) %>%
  arrange(year)

roll_mean <- function(x, width = 15) {
  rollapply(x, width = width, fill = NA, align = "right", FUN = mean)
}
# roll_sd() is defined in Section 3 and reused here.

winter.cv <- winter.raw %>%
  mutate(GOA.cv = roll_sd(GOA) / roll_mean(GOA),
         EBS.cv = roll_sd(EBS) / roll_mean(EBS))

p.goa.cv <- ggplot(winter.cv, aes(x = year, y = GOA.cv)) +
  geom_line(color = "darkorange4") +
  geom_point(color = "darkorange4", size = 1.4) +
  labs(title = "Gulf of Alaska", x = "Year", y = "Winter SST CV (SD / mean)") +
  theme_bw(base_size = 11)

p.ebs.cv <- ggplot(winter.cv, aes(x = year, y = EBS.cv)) +
  geom_line(color = "steelblue4") +
  geom_point(color = "steelblue4", size = 1.4) +
  labs(title = "Eastern Bering Sea", x = "Year", y = "Winter SST CV (SD / mean)") +
  theme_bw(base_size = 11)

p.cv <- (p.goa.cv / p.ebs.cv) +
  plot_annotation(
    title = "15-yr Rolling Winter SST CV (right-aligned)",
    theme = theme(plot.title = element_text(size = 12, hjust = 0.5))
  )

ggsave("./Figures/ERA5_SST_CV_15yr_rolling.png",
       plot = p.cv, width = 7, height = 6.5, dpi = 300)
message("Saved: Figures/ERA5_SST_CV_15yr_rolling.png")

# ============================================================
# SECTION 4: SLP EOF1 LOADINGS MAP — MERCATOR PROJECTION
# ============================================================
# Reproduces the EOF1 analysis from ERA5_SLP_download.R and plots
# the spatial loadings on a Mercator projection, with the Aleutian
# Low box (192.5-207.5E, 45-55N; Litzow et al. 2020 PNAS) overlaid.
# Sign is flipped if needed so the AL box region holds positive loadings.

library(irlba)

slp.file <- if (file.exists("./Data/era5_slp_NP_1950_2026.nc")) {
  "./Data/era5_slp_NP_1950_2026.nc"
} else {
  "./Data/era5_slp_NP_1950_2025.nc"
}
message("Loading SLP from: ", slp.file)

nc.slp   <- nc_open(slp.file)
slp      <- ncvar_get(nc.slp, "msl")
lons.slp <- ncvar_get(nc.slp, "longitude")
lats.slp <- ncvar_get(nc.slp, "latitude")
time.slp <- ncvar_get(nc.slp, "valid_time")
nc_close(nc.slp)

slp        <- slp / 100   # Pa → hPa
dates.slp  <- as.Date(as.POSIXct(time.slp, origin = "1970-01-01", tz = "UTC"))
months.slp <- as.integer(format(dates.slp, "%m"))
years.slp  <- as.integer(format(dates.slp, "%Y"))
n.time.slp <- dim(slp)[3]

# --- Step 1: monthly anomalies (1950-1979 climatology) ---
slp.anom <- array(NA_real_, dim = dim(slp))
base.slp <- years.slp >= 1950 & years.slp <= 1979

for (m in 1:12) {
  t.base <- which(base.slp & months.slp == m)
  t.all  <- which(months.slp == m)
  clim.m <- apply(slp[, , t.base, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  for (t in t.all) slp.anom[, , t] <- slp[, , t] - clim.m
}

# --- Step 2: detrend each cell ---
X.slp    <- cbind(1, seq_len(n.time.slp))
IH.slp   <- diag(n.time.slp) - X.slp %*% solve(t(X.slp) %*% X.slp) %*% t(X.slp)
anom.mat.slp  <- matrix(slp.anom,
                        nrow = length(lons.slp) * length(lats.slp),
                        ncol = n.time.slp)
detrend.slp <- anom.mat.slp %*% t(IH.slp)

# --- Step 3: area-weight and truncated SVD for EOF1 only ---
grid.slp <- expand.grid(lon = lons.slp, lat = lats.slp)
w.slp    <- sqrt(cos(grid.slp$lat * pi / 180))
sv.slp   <- irlba(sweep(detrend.slp, 1, w.slp, "*"), nv = 1)
eof1.slp <- sv.slp$u[, 1] / w.slp

# Flip sign so AL box (192.5-207.5E, 45-55N; Litzow et al. 2020 PNAS)
# has positive mean loading.
# grid.slp$lon is in 0-360 space (ERA5 returns 120-250E for this domain).
al.box.mask <- grid.slp$lon >= 192.5 & grid.slp$lon <= 207.5 &
               grid.slp$lat >=  45   & grid.slp$lat <=  55
al.box.mean <- mean(eof1.slp[al.box.mask], na.rm = TRUE)
if (!is.na(al.box.mean) && al.box.mean < 0) eof1.slp <- -eof1.slp

# Convert lons to 0-360 for Pacific-centered plot
eof.df.slp <- grid.slp %>%
  mutate(lon  = ifelse(lon < 0, lon + 360, lon),
         EOF1 = eof1.slp)

# AL box in 0-360 space (Litzow et al. 2020 PNAS: 192.5-207.5E, 45-55N)
al.box.df <- data.frame(
  x = c(192.5, 207.5, 207.5, 192.5, 192.5),
  y = c(45,    45,    55,    55,    45)
)

mapWorld.slp <- map_data("world", wrap = c(20, 380))

eof1.map <- ggplot() +
  geom_raster(data = eof.df.slp, aes(x = lon, y = lat, fill = EOF1)) +
  geom_polygon(data = mapWorld.slp,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "gray30", linewidth = 0.2) +
  geom_path(data = al.box.df, aes(x = x, y = y),
            color = "black", linewidth = 0.9) +
  scale_fill_distiller(palette = "RdBu", direction = -1, name = "Loading") +
  coord_map(projection = "rectangular",
            parameters = 55,
            xlim = c(120, 250), ylim = c(20, 66)) +
  scale_x_continuous(breaks = seq(120, 250, 30),
                     labels = lon_label) +
  scale_y_continuous(breaks = seq(20, 60, 10),
                     labels = function(y) paste0(y, "\u00b0N")) +
  labs(title = "ERA5 SLP EOF1 Loadings",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3))

ggsave("./Figures/ERA5_SLP_EOF1_loadings.png",
       plot = eof1.map, width = 8, height = 6, dpi = 300)
message("Saved: Figures/ERA5_SLP_EOF1_loadings.png")

# ============================================================
# SECTION 5: AL SLP SD vs SST AR(1) — 15-yr dual-axis plot
# ============================================================
# Computes AL winter SLP anomaly inline from the detrended monthly
# anomaly matrix built in Section 4 (al.box.mask = Litzow et al. 2020
# PNAS box, 192.5-207.5E / 45-55N). Calculates 15-yr rolling SD
# (z-scored), joins with SST AR(1) from Section 3, and plots
# dual-axis time series for GOA and EBS.

al.monthly <- data.frame(
  year  = years.slp,
  month = months.slp,
  SLP   = apply(detrend.slp[al.box.mask, ], 2, mean, na.rm = TRUE)
)

al <- al.monthly %>%
  filter(month %in% c(11, 12, 1, 2, 3)) %>%
  mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
  group_by(win.year) %>%
  summarise(SLP = mean(SLP, na.rm = TRUE), .groups = "drop") %>%
  rename(year = win.year) %>%
  arrange(year) %>%
  mutate(AL.sd = roll_sd(SLP),
         AL.sd = as.numeric(scale(AL.sd)))

al.sst <- winter.roll %>%
  select(year, GOA.ar1, EBS.ar1) %>%
  left_join(al %>% select(year, AL.sd), by = "year") %>%
  filter(!is.na(AL.sd))

al.sd.range   <- range(al.sst$AL.sd,   na.rm = TRUE)
goa.ar1.range <- range(al.sst$GOA.ar1, na.rm = TRUE)
ebs.ar1.range <- range(al.sst$EBS.ar1, na.rm = TRUE)

# --- Modified Chelton method (Pyper & Peterman 1998) ---
# Effective sample size accounts for autocorrelation in both series:
#   1/N* ≈ 1/N + (2/N) * Σ_{j=1}^{N/5} (1 - j/N) * ρ_xx(j) * ρ_yy(j)
# t = r * sqrt(N* - 2) / sqrt(1 - r²), one-sided p from t-distribution.
chelton_test <- function(x, y) {
  n       <- length(x)
  r.obs   <- cor(x, y, use = "complete.obs")
  max.lag <- floor(n / 5)

  acf.x <- acf(x, lag.max = max.lag, plot = FALSE, na.action = na.pass)$acf[-1]
  acf.y <- acf(y, lag.max = max.lag, plot = FALSE, na.action = na.pass)$acf[-1]

  j     <- seq_len(max.lag)
  n.eff <- n / (1 + 2 * sum((1 - j / n) * acf.x * acf.y))
  n.eff <- max(3, min(n, n.eff))   # clamp to [3, n]

  t.obs <- r.obs * sqrt(n.eff - 2) / sqrt(1 - r.obs^2)
  df    <- n.eff - 2
  pval  <- pt(t.obs, df = df, lower.tail = FALSE)   # one-sided

  list(r.obs = r.obs, t.obs = t.obs, n.eff = n.eff, df = df, pval = pval)
}

dat.goa  <- al.sst %>% filter(!is.na(GOA.ar1), !is.na(AL.sd))
dat.ebs  <- al.sst %>% filter(!is.na(EBS.ar1), !is.na(AL.sd))

chel.goa <- chelton_test(dat.goa$GOA.ar1, dat.goa$AL.sd)
chel.ebs <- chelton_test(dat.ebs$EBS.ar1, dat.ebs$AL.sd)

message(sprintf("GOA Chelton: r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
                chel.goa$r.obs, chel.goa$n.eff, chel.goa$t.obs, chel.goa$pval))
message(sprintf("EBS Chelton: r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
                chel.ebs$r.obs, chel.ebs$n.eff, chel.ebs$t.obs, chel.ebs$pval))

make_dual_axis_plot <- function(dat, ar1.col, ar1.color, region.label, ar1.range, chel) {
  k <- diff(ar1.range) / diff(al.sd.range)
  b <- ar1.range[1] - k * al.sd.range[1]

  dat <- dat %>% filter(!is.na(.data[[ar1.col]]), !is.na(AL.sd))

  lbl <- sprintf("r = %.2f, N* = %.1f\nt = %.2f, p = %.3f",
                 chel$r.obs, chel$n.eff, chel$t.obs, chel$pval)

  ggplot(dat, aes(x = year)) +
    geom_line(aes(y = .data[[ar1.col]], color = "SST AR(1)"), linewidth = 0.7) +
    geom_point(aes(y = .data[[ar1.col]], color = "SST AR(1)"), size = 1.5) +
    geom_line(aes(y = k * AL.sd + b, color = "AL SLP SD"),
              linewidth = 0.7, linetype = "dashed") +
    geom_point(aes(y = k * AL.sd + b, color = "AL SLP SD"), size = 1.5, shape = 1) +
    annotate("text", x = Inf, y = -Inf, label = lbl,
             hjust = 1.1, vjust = -0.4, size = 3.5, lineheight = 1.3) +
    scale_y_continuous(
      name     = "SST AR(1)",
      sec.axis = sec_axis(~ (. - b) / k, name = "AL SLP SD (z-scored)")
    ) +
    scale_color_manual(values = c("SST AR(1)" = ar1.color, "AL SLP SD" = "gray30")) +
    labs(title = region.label, x = "Year", color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position    = "bottom",
          axis.title.y.left  = element_text(color = ar1.color),
          axis.text.y.left   = element_text(color = ar1.color),
          axis.title.y.right = element_text(color = "gray30"),
          axis.text.y.right  = element_text(color = "gray30"))
}

p.goa.dual <- make_dual_axis_plot(al.sst, "GOA.ar1", "darkorange4", "GOA", goa.ar1.range, chel.goa)
p.ebs.dual <- make_dual_axis_plot(al.sst, "EBS.ar1", "steelblue4", "EBS", ebs.ar1.range, chel.ebs)

p.dual <- p.goa.dual / p.ebs.dual +
  plot_annotation(
    title = "15-yr Rolling AL SLP SD vs. SST AR(1) \u2014 Winter (Nov\u2013Mar)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 12))
  )

ggsave("./Figures/AL_SD_SST_AR1_15yr_dual_axis.png",
       plot = p.dual, width = 7, height = 7, dpi = 300)
message("Saved: Figures/AL_SD_SST_AR1_15yr_dual_axis.png")

make_chelton_panel <- function(t.obs, df, n.eff, pval, region.label) {
  t.range <- seq(-4, max(4, t.obs + 0.5), length.out = 500)
  lbl <- sprintf("p = %.4f\nN* = %.1f", pval, n.eff)
  ggplot(data.frame(t = t.range, d = dt(t.range, df = df)), aes(x = t, y = d)) +
    geom_area(fill = "gray80", color = "gray50", linewidth = 0.3) +
    geom_vline(xintercept = t.obs, linetype = "dashed", linewidth = 0.9) +
    annotate("text", x = Inf, y = Inf, label = lbl,
             hjust = 1.1, vjust = 1.5, size = 3.5, lineheight = 1.3) +
    labs(title = paste(region.label, "— Chelton method"),
         x = sprintf("t  (df = %.1f)", df), y = "Density") +
    theme_bw(base_size = 11)
}

p.goa.chel <- make_chelton_panel(chel.goa$t.obs, chel.goa$df, chel.goa$n.eff, chel.goa$pval, "GOA")
p.ebs.chel <- make_chelton_panel(chel.ebs$t.obs, chel.ebs$df, chel.ebs$n.eff, chel.ebs$pval, "EBS")

p.chelton <- p.goa.chel / p.ebs.chel +
  plot_annotation(
    title = "Chelton significance test: AL SLP SD vs. SST AR(1)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 12))
  )

ggsave("./Figures/AL_SD_SST_AR1_chelton_dist.png",
       plot = p.chelton, width = 5, height = 7, dpi = 300)
message("Saved: Figures/AL_SD_SST_AR1_chelton_dist.png")

# ============================================================
# SECTION 5B: MLD EOF MODES 1-2 — WINTER (NDJFM)
# ============================================================
# Leading two EOF modes of winter (NDJFM) mixed-layer depth over the
# N. Pacific basin domain used elsewhere (110-250E, 20-66N; the SLP
# EOF1 / MLD-regression domain). Built from the same detrended cellwise
# winter-month MLD anomalies used in the regression sections
# (oras5_mld_NP_detr_anom_winmon.nc, already detrended per cell), here
# averaged to NDJFM winter-year means (year = January year). Cells are
# area-weighted by sqrt(cos(lat)); only the first five modes are
# extracted (truncated SVD, irlba). Outputs:
#   - North et al. (1982) eigenvalue error-bar plot, modes 1-5
#   - 4-panel figure: EOF1/EOF2 loadings + PC1/PC2 time series
# NOTE: domain matches Section 6 (includes the western Pacific). If
# EOF1 is dominated by the high-variance NW Pacific, trim to >=160E.

mld.winmon.file    <- "./Data/oras5_mld_NP_detr_anom_winmon.nc"
mld.eof.domain     <- list(lon = c(110, 250), lat = c(20, 66))
mld.eof.eig.cache  <- "./Output/mld_eof_eigenvalues.csv"
mld.eof.load.cache <- "./Output/mld_eof_loadings.csv"
mld.eof.pc.cache   <- "./Output/mld_eof_PCs.csv"

if (all(file.exists(c(mld.eof.eig.cache, mld.eof.load.cache, mld.eof.pc.cache)))) {
  message("Loading cached MLD EOF results...")
  north.df <- read.csv(mld.eof.eig.cache)
  load.df  <- read.csv(mld.eof.load.cache)
  pc.df    <- read.csv(mld.eof.pc.cache)
} else {
  message("Loading detrended winter-month MLD anomalies for EOF...")
  nc.e   <- nc_open(mld.winmon.file)
  mld.e  <- ncvar_get(nc.e, "somxl030")
  lons.e <- ncvar_get(nc.e, "nav_lon")
  lats.e <- ncvar_get(nc.e, "nav_lat")
  time.e <- ncvar_get(nc.e, "time_counter")
  tun.e  <- ncatt_get(nc.e, "time_counter", "units")$value
  nc_close(nc.e)

  mld.e[mld.e > 1e10] <- NA
  origin.e <- sub(".*since ", "", tun.e)
  dates.e  <- as.Date(as.POSIXct(time.e, origin = origin.e, tz = "UTC"))
  years.e  <- as.integer(format(dates.e, "%Y"))
  months.e <- as.integer(format(dates.e, "%m"))
  winyr.e  <- ifelse(months.e %in% c(11, 12), years.e + 1L, years.e)
  lon.e360 <- ifelse(lons.e < 0, lons.e + 360, lons.e)

  nx.e <- dim(mld.e)[1]; ny.e <- dim(mld.e)[2]; nt.e <- dim(mld.e)[3]
  mat.e <- matrix(mld.e, nrow = nx.e * ny.e, ncol = nt.e)
  rm(mld.e); gc()

  # --- NDJFM winter-year means (year = January year) ---
  uwy   <- sort(unique(winyr.e))
  win.e <- matrix(NA_real_, nrow = nx.e * ny.e, ncol = length(uwy))
  for (k in seq_along(uwy)) {
    ti <- which(winyr.e == uwy[k])
    win.e[, k] <- rowMeans(mat.e[, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(mat.e); gc()

  # Keep only winters with full Nov-Mar coverage so the time axis is balanced
  mo.per.wy <- tapply(months.e, winyr.e, function(m) length(unique(m)))
  full.wy   <- as.integer(names(mo.per.wy)[mo.per.wy >= 5])
  keep.t    <- which(uwy %in% full.wy)
  win.e     <- win.e[, keep.t, drop = FALSE]
  eof.years <- uwy[keep.t]
  n.time.e  <- length(eof.years)

  # --- Restrict to basin domain and fully-observed ocean cells (SVD needs no NA) ---
  coord.lon <- as.vector(lon.e360)
  coord.lat <- as.vector(lats.e)
  in.dom    <- coord.lon >= mld.eof.domain$lon[1] & coord.lon <= mld.eof.domain$lon[2] &
               coord.lat >= mld.eof.domain$lat[1] & coord.lat <= mld.eof.domain$lat[2]
  ocean     <- rowSums(is.na(win.e)) == 0
  use.cell  <- which(in.dom & ocean)
  message(sprintf("MLD EOF: %d ocean cells x %d winters (%d-%d)",
                  length(use.cell), n.time.e, min(eof.years), max(eof.years)))

  # --- Center each cell in time, area-weight, truncated SVD (5 modes) ---
  A    <- win.e[use.cell, , drop = FALSE]
  A    <- A - rowMeans(A)                          # covariance EOF (center in time)
  w.e  <- sqrt(cos(coord.lat[use.cell] * pi / 180))
  Aw   <- sweep(A, 1, w.e, "*")
  sv.e <- irlba(Aw, nv = 5, nu = 5)

  # Eigenvalues, % variance, North et al. (1982) sampling error
  lambda  <- sv.e$d^2 / (n.time.e - 1)
  tot.var <- sum(Aw^2) / (n.time.e - 1)            # trace of weighted covariance
  pct.var <- 100 * lambda / tot.var
  N.eff   <- n.time.e                              # # winters (independent realizations)
  pct.err <- pct.var * sqrt(2 / N.eff)             # dlambda ~ lambda * sqrt(2/N)

  north.df <- data.frame(mode = 1:5, pct = pct.var, err = pct.err, N = N.eff)

  # --- De-weighted loadings + PC time series for modes 1-2 (fix arbitrary sign) ---
  flip <- function(load, pc) {
    i <- which.max(abs(load))
    if (load[i] < 0) list(load = -load, pc = -pc) else list(load = load, pc = pc)
  }
  m1 <- flip(sv.e$u[, 1] / w.e, sv.e$d[1] * sv.e$v[, 1])
  m2 <- flip(sv.e$u[, 2] / w.e, sv.e$d[2] * sv.e$v[, 2])

  # Loadings to 0.5-deg grid for raster (ORAS5 grid is irregular)
  load.df <- data.frame(lon = coord.lon[use.cell], lat = coord.lat[use.cell],
                        EOF1 = m1$load, EOF2 = m2$load) %>%
    mutate(lon = round(lon / 0.5) * 0.5, lat = round(lat / 0.5) * 0.5) %>%
    group_by(lon, lat) %>%
    summarise(EOF1 = mean(EOF1, na.rm = TRUE),
              EOF2 = mean(EOF2, na.rm = TRUE), .groups = "drop")

  pc.df <- data.frame(year = eof.years, PC1 = m1$pc, PC2 = m2$pc)

  write.csv(north.df, mld.eof.eig.cache,  row.names = FALSE)
  write.csv(load.df,  mld.eof.load.cache, row.names = FALSE)
  write.csv(pc.df,    mld.eof.pc.cache,   row.names = FALSE)
  message("Saved: Output/mld_eof_{eigenvalues,loadings,PCs}.csv")
}

# Common plotting variables
pct.var   <- north.df$pct
eof.years <- pc.df$year
N.eff     <- north.df$N[1]

# --- North et al. (1982) eigenvalue error-bar plot (modes 1-5) ---
p.north <- ggplot(north.df, aes(x = factor(mode), y = pct)) +
  geom_point(size = 2.4, color = "steelblue4") +
  geom_errorbar(aes(ymin = pct - err, ymax = pct + err),
                width = 0.18, color = "steelblue4") +
  labs(title = "MLD EOF eigenvalue spectrum — North et al. (1982) error bars",
       subtitle = sprintf("Winter (NDJFM) %d–%d, N = %d winters",
                          min(eof.years), max(eof.years), N.eff),
       x = "EOF mode", y = "Variance explained (%)") +
  theme_bw(base_size = 11)

ggsave("./Figures/MLD_EOF_North_errorbars.png",
       plot = p.north, width = 6, height = 4.5, dpi = 300)
message("Saved: Figures/MLD_EOF_North_errorbars.png")

# --- 4-panel: EOF1/EOF2 loadings + PC1/PC2 time series ---
goa.path   <- data.frame(x = c(goa.x, goa.x[1]), y = c(goa.y, goa.y[1]))
ebs.path   <- data.frame(x = c(ebs.x, ebs.x[1]), y = c(ebs.y, ebs.y[1]))
mapWorld.e <- map_data("world2")

loadings_map <- function(col, ttl, pctv) {
  lim <- max(abs(load.df[[col]]), na.rm = TRUE)
  ggplot() +
    geom_raster(data = load.df, aes(x = lon, y = lat, fill = .data[[col]])) +
    geom_polygon(data = mapWorld.e, aes(x = long, y = lat, group = group),
                 fill = "gray85", color = "gray40", linewidth = 0.2) +
    geom_path(data = goa.path, aes(x = x, y = y), color = "black", linewidth = 0.4) +
    geom_path(data = ebs.path, aes(x = x, y = y), color = "black", linewidth = 0.4) +
    scale_fill_gradient2(low = "steelblue4", mid = "white", high = "darkorange4",
                         midpoint = 0, limits = c(-lim, lim), name = "Loading") +
    coord_map(projection = "rectangular", parameters = 55,
              xlim = mld.eof.domain$lon, ylim = mld.eof.domain$lat) +
    scale_x_continuous(breaks = seq(120, 240, 30), labels = lon_label) +
    scale_y_continuous(breaks = seq(20, 60, 10),
                       labels = function(y) paste0(y, "°N")) +
    labs(title = sprintf("%s (%.1f%% var)", ttl, pctv), x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())
}

pc_ts_plot <- function(col, ttl, color) {
  ggplot(pc.df, aes(x = year, y = .data[[col]])) +
    geom_hline(yintercept = 0, color = "gray60", linewidth = 0.3) +
    geom_line(color = color, linewidth = 0.6) +
    geom_point(color = color, size = 1.2) +
    labs(title = ttl, x = "Year", y = "PC (std. units)") +
    theme_bw(base_size = 10)
}

p.mld.eof <- (loadings_map("EOF1", "MLD EOF1", pct.var[1]) |
              loadings_map("EOF2", "MLD EOF2", pct.var[2])) /
             (pc_ts_plot("PC1", "MLD PC1", "darkorange4") |
              pc_ts_plot("PC2", "MLD PC2", "steelblue4")) +
  plot_annotation(
    title = sprintf("Winter (NDJFM) MLD EOF modes 1–2, %d–%d",
                    min(eof.years), max(eof.years)),
    theme = theme(plot.title = element_text(size = 12, hjust = 0.5))
  )

ggsave("./Figures/MLD_EOF1-2_loadings_PCs.png",
       plot = p.mld.eof, width = 11, height = 7.5, dpi = 300)
message("Saved: Figures/MLD_EOF1-2_loadings_PCs.png")

# ============================================================
# SECTION 6: WINTER MLD REGRESSED ON WINTER AL SLP (full EOF domain)
# ============================================================
# Cellwise regression of detrended winter-mean MLD anomaly (y) on the
# annual Nov-Mar mean AL-box SLP anomaly (x; al$SLP from Section 5,
# single-year winter means, year = January year). Spatial domain matches
# the SLP EOF1 domain (20-66N, 110-250E). Full series through winter 2026.
# GLS with AR(1) residuals (nlme::gls + corAR1), then Benjamini-Hochberg
# FDR (Wilks 2016, BAMS). Contour drawn around cells with q <= 0.05.

library(nlme)

mld.winmon.file <- "./Data/oras5_mld_NP_detr_anom_winmon.nc"
mld.reg.cache   <- "./Output/mld_al_slp_regression.csv"

if (file.exists(mld.reg.cache)) {
  message("Loading cached MLD ~ AL SLP SD regression from: ", mld.reg.cache)
  mld.reg.df <- read.csv(mld.reg.cache)
} else {
  message("Loading detrended winter-month MLD anomalies (full NP domain)...")
  nc.m   <- nc_open(mld.winmon.file)
  mld.w  <- ncvar_get(nc.m, "somxl030")
  lons.m <- ncvar_get(nc.m, "nav_lon")
  lats.m <- ncvar_get(nc.m, "nav_lat")
  time.m <- ncvar_get(nc.m, "time_counter")
  tun.m  <- ncatt_get(nc.m, "time_counter", "units")$value
  nc_close(nc.m)

  mld.w[mld.w > 1e10] <- NA
  origin.m  <- sub(".*since ", "", tun.m)
  dates.m   <- as.Date(as.POSIXct(time.m, origin = origin.m, tz = "UTC"))
  years.m   <- as.integer(format(dates.m, "%Y"))
  months.m  <- as.integer(format(dates.m, "%m"))
  winyrs.m  <- ifelse(months.m %in% c(11, 12), years.m + 1L, years.m)
  lons.m360 <- ifelse(lons.m < 0, lons.m + 360, lons.m)

  nx <- dim(mld.w)[1]; ny <- dim(mld.w)[2]; nt <- dim(mld.w)[3]
  mat <- matrix(mld.w, nrow = nx * ny, ncol = nt)
  rm(mld.w); gc()

  uy <- sort(unique(winyrs.m))
  win.mld <- matrix(NA_real_, nrow = nx * ny, ncol = length(uy))
  for (k in seq_along(uy)) {
    ti <- which(winyrs.m == uy[k])
    win.mld[, k] <- rowMeans(mat[, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(mat); gc()

  # Align years with al$SLP (annual Nov-Mar mean AL-box SLP anomaly, single-year)
  al.df   <- al %>% select(year, SLP) %>% filter(!is.na(SLP))
  yr.use  <- intersect(uy, al.df$year)
  x.vec   <- al.df$SLP[match(yr.use, al.df$year)]
  win.mld <- win.mld[, match(yr.use, uy), drop = FALSE]
  n.yrs   <- length(yr.use)
  message(sprintf("Regression years: %d-%d (%d years)",
                  min(yr.use), max(yr.use), n.yrs))

  good.cells <- which(rowSums(!is.na(win.mld)) >= max(10, 0.5 * n.yrs))
  message(sprintf("Fitting GLS AR(1) at %d cells over %d years (slow)...",
                  length(good.cells), n.yrs))

  beta.vec <- rep(NA_real_, length(good.cells))
  pval.vec <- rep(NA_real_, length(good.cells))
  df.reg   <- data.frame(y = NA_real_, x = x.vec)

  for (k in seq_along(good.cells)) {
    y <- win.mld[good.cells[k], ]
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

  beta.full <- rep(NA_real_, nx * ny)
  pval.full <- rep(NA_real_, nx * ny)
  beta.full[good.cells] <- beta.vec
  pval.full[good.cells] <- pval.vec

  mld.reg.df <- data.frame(
    lon  = as.vector(lons.m360),
    lat  = as.vector(lats.m),
    beta = beta.full,
    pval = pval.full
  ) %>% filter(!is.na(beta))

  write.csv(mld.reg.df, mld.reg.cache, row.names = FALSE)
  message("Saved: ", mld.reg.cache)
}

# Restrict to SLP EOF1 domain
mld.reg.df <- mld.reg.df %>%
  filter(lon >= 110, lon <= 250, lat >= 20, lat <= 66)

# --- Benjamini-Hochberg FDR (Wilks 2016, BAMS) ---
mld.reg.df$qval <- p.adjust(mld.reg.df$pval, method = "BH")

n.sig.raw <- sum(mld.reg.df$pval <= 0.05, na.rm = TRUE)
n.sig.fdr <- sum(mld.reg.df$qval <= 0.05, na.rm = TRUE)
message(sprintf("Raw p<=0.05: %d cells | BH-FDR q<=0.05: %d cells (of %d)",
                n.sig.raw, n.sig.fdr, nrow(mld.reg.df)))

# Regular 0.5-deg grid for continuous raster + FDR contouring (ORAS5 is irregular)
mld.reg.grid <- mld.reg.df %>%
  mutate(lon = round(lon / 0.5) * 0.5,
         lat = round(lat / 0.5) * 0.5) %>%
  group_by(lon, lat) %>%
  summarise(beta = mean(beta, na.rm = TRUE),
            qval = mean(qval, na.rm = TRUE),
            .groups = "drop") %>%
  complete(lon = seq(110, 250, by = 0.5),
           lat = seq(20,  66,  by = 0.5))

col.lim.mld <- max(abs(mld.reg.df$beta), na.rm = TRUE)

# Clean 0-360 world map (avoids horizontal wrap artifacts from map_data("world", wrap=...))
mapWorld.clean <- map_data("world2")

p.mld.reg <- ggplot() +
  geom_raster(data = filter(mld.reg.grid, !is.na(beta)),
              aes(x = lon, y = lat, fill = beta)) +
  geom_contour(data = filter(mld.reg.grid, !is.na(qval)),
               aes(x = lon, y = lat, z = qval),
               breaks = 0.05, color = "black", linewidth = 0.4) +
  geom_polygon(data = mapWorld.clean,
               aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray30", linewidth = 0.25) +
  geom_path(data = al.box.df, aes(x = x, y = y),
            color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue4", mid = "white",
                       high = "darkorange4", midpoint = 0,
                       limits = c(-col.lim.mld, col.lim.mld),
                       name = "\u03b2 (m / hPa)") +
  coord_map(projection = "rectangular", parameters = 55,
            xlim = c(110, 250), ylim = c(20, 66)) +
  scale_x_continuous(breaks = seq(120, 250, 30), labels = lon_label) +
  scale_y_continuous(breaks = seq(20, 60, 10),
                     labels = function(y) paste0(y, "\u00b0N")) +
  labs(title    = "Winter MLD regressed on winter AL SLP (annual Nov\u2013Mar means)",
       subtitle = "GLS AR(1); contour = Benjamini-Hochberg FDR q \u2264 0.05 (Wilks 2016)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/MLD_AL_SLP_regression.png",
       plot = p.mld.reg, width = 9, height = 6, dpi = 300)
message("Saved: Figures/MLD_AL_SLP_regression.png")

# ============================================================
# SECTION 6B: WINTER SST REGRESSED ON WINTER AL SLP (full EOF domain)
# ============================================================
# Mirror of Section 6 with cellwise winter-mean SST anomaly (y) as the
# response instead of MLD. x is the same annual Nov-Mar mean AL-box SLP
# anomaly (al$SLP from Section 5). SST cells: monthly anomalies vs the
# 1950-1979 climatology, averaged to NDJFM winter-year means (year =
# January), then detrended per cell (matching the MLD pre-detrending and
# the Section 3 SST convention). GLS with AR(1) residuals; BH-FDR contour
# (Wilks 2016). Reuses the in-memory ERA5 `sst` array from Section 2.

sst.reg.cache <- "./Output/sst_al_slp_regression.csv"

if (file.exists(sst.reg.cache)) {
  message("Loading cached SST ~ AL SLP regression from: ", sst.reg.cache)
  sst.reg.df <- read.csv(sst.reg.cache)
} else {
  message("Computing cellwise winter-mean SST anomalies...")
  nx.s <- dim(sst)[1]; ny.s <- dim(sst)[2]
  base.s  <- years >= 1950 & years <= 1979
  winyr.s <- ifelse(months %in% c(11, 12), years + 1L, years)

  # Monthly climatology (1950-1979) per cell
  clim.s <- vector("list", 12)
  for (m in 1:12) {
    tb <- which(base.s & months == m)
    clim.s[[m]] <- apply(sst[, , tb, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  }

  # NDJFM winter-year mean anomaly per cell (built slice-wise to limit memory)
  uy.s    <- sort(unique(winyr.s))
  win.sst <- matrix(NA_real_, nrow = nx.s * ny.s, ncol = length(uy.s))
  for (k in seq_along(uy.s)) {
    ti <- which(winyr.s == uy.s[k] & months %in% c(11, 12, 1, 2, 3))
    if (!length(ti)) next
    acc <- matrix(0, nx.s, ny.s); cnt <- matrix(0, nx.s, ny.s)
    for (t in ti) {
      a    <- sst[, , t] - clim.s[[months[t]]]
      good <- is.finite(a)
      acc[good] <- acc[good] + a[good]
      cnt[good] <- cnt[good] + 1
    }
    wm <- acc / cnt; wm[cnt == 0] <- NA
    win.sst[, k] <- as.vector(wm)
  }

  # Align years with al$SLP (annual Nov-Mar mean AL-box SLP anomaly)
  al.df   <- al %>% select(year, SLP) %>% filter(!is.na(SLP))
  yr.use  <- intersect(uy.s, al.df$year)
  x.vec   <- al.df$SLP[match(yr.use, al.df$year)]
  win.sst <- win.sst[, match(yr.use, uy.s), drop = FALSE]
  n.yrs   <- length(yr.use)
  message(sprintf("Regression years: %d-%d (%d years)",
                  min(yr.use), max(yr.use), n.yrs))

  good.cells <- which(rowSums(!is.na(win.sst)) >= max(10, 0.5 * n.yrs))
  message(sprintf("Fitting GLS AR(1) at %d cells over %d years (slow)...",
                  length(good.cells), n.yrs))

  beta.vec <- rep(NA_real_, length(good.cells))
  pval.vec <- rep(NA_real_, length(good.cells))
  df.reg   <- data.frame(y = NA_real_, x = x.vec)
  tt       <- seq_len(n.yrs)

  for (k in seq_along(good.cells)) {
    y  <- win.sst[good.cells[k], ]
    ok <- is.finite(y)
    if (sum(ok) < 10) next
    yd     <- rep(NA_real_, n.yrs)            # detrend the winter series per cell
    yd[ok] <- residuals(lm(y[ok] ~ tt[ok]))
    df.reg$y <- yd
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

  beta.full <- rep(NA_real_, nx.s * ny.s)
  pval.full <- rep(NA_real_, nx.s * ny.s)
  beta.full[good.cells] <- beta.vec
  pval.full[good.cells] <- pval.vec

  grid.s   <- expand.grid(lon = lons, lat = lats)
  lon.s360 <- ifelse(grid.s$lon < 0, grid.s$lon + 360, grid.s$lon)

  sst.reg.df <- data.frame(
    lon  = lon.s360,
    lat  = grid.s$lat,
    beta = beta.full,
    pval = pval.full
  ) %>% filter(!is.na(beta))

  write.csv(sst.reg.df, sst.reg.cache, row.names = FALSE)
  message("Saved: ", sst.reg.cache)
}

# Restrict to SLP EOF1 domain
sst.reg.df <- sst.reg.df %>%
  filter(lon >= 110, lon <= 250, lat >= 20, lat <= 66)

# --- Benjamini-Hochberg FDR (Wilks 2016, BAMS) ---
sst.reg.df$qval <- p.adjust(sst.reg.df$pval, method = "BH")

n.sig.raw.s <- sum(sst.reg.df$pval <= 0.05, na.rm = TRUE)
n.sig.fdr.s <- sum(sst.reg.df$qval <= 0.05, na.rm = TRUE)
message(sprintf("Raw p<=0.05: %d cells | BH-FDR q<=0.05: %d cells (of %d)",
                n.sig.raw.s, n.sig.fdr.s, nrow(sst.reg.df)))

# Regular 0.5-deg grid for continuous raster + FDR contouring
sst.reg.grid <- sst.reg.df %>%
  mutate(lon = round(lon / 0.5) * 0.5,
         lat = round(lat / 0.5) * 0.5) %>%
  group_by(lon, lat) %>%
  summarise(beta = mean(beta, na.rm = TRUE),
            qval = mean(qval, na.rm = TRUE),
            .groups = "drop") %>%
  complete(lon = seq(110, 250, by = 0.5),
           lat = seq(20,  66,  by = 0.5))

col.lim.sst <- max(abs(sst.reg.df$beta), na.rm = TRUE)

p.sst.reg <- ggplot() +
  geom_raster(data = filter(sst.reg.grid, !is.na(beta)),
              aes(x = lon, y = lat, fill = beta)) +
  geom_contour(data = filter(sst.reg.grid, !is.na(qval)),
               aes(x = lon, y = lat, z = qval),
               breaks = 0.05, color = "black", linewidth = 0.4) +
  geom_polygon(data = mapWorld.clean,
               aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray30", linewidth = 0.25) +
  geom_path(data = al.box.df, aes(x = x, y = y),
            color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue4", mid = "white",
                       high = "darkorange4", midpoint = 0,
                       limits = c(-col.lim.sst, col.lim.sst),
                       name = "β (°C / hPa)") +
  coord_map(projection = "rectangular", parameters = 55,
            xlim = c(110, 250), ylim = c(20, 66)) +
  scale_x_continuous(breaks = seq(120, 250, 30), labels = lon_label) +
  scale_y_continuous(breaks = seq(20, 60, 10),
                     labels = function(y) paste0(y, "°N")) +
  labs(title    = "Winter SST regressed on winter AL SLP (annual Nov–Mar means)",
       subtitle = "GLS AR(1); contour = Benjamini-Hochberg FDR q ≤ 0.05 (Wilks 2016)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/SST_AL_SLP_regression.png",
       plot = p.sst.reg, width = 9, height = 6, dpi = 300)
message("Saved: Figures/SST_AL_SLP_regression.png")

# ============================================================
# SECTION 6C: COMBINED SST + MLD ~ AL SLP REGRESSION (NE Pacific)
# ============================================================
# Stacks the Section 6B SST regression (top) and Section 6 MLD regression
# (bottom) into one figure, restricted to the NE Pacific (east of 150E).
# Each panel keeps its own diverging scale (units differ: °C/hPa vs m/hPa)
# and colour limits computed from the shown (east-of-150E) cells only.
# Reuses sst.reg.df (Section 6B) and mld.reg.df (Section 6); both carry
# cellwise beta + BH-FDR qval. Same GLS AR(1) / FDR workflow throughout.

nep_reg_panel <- function(df, fill.lab, ttl) {
  d <- df %>% filter(lon >= 150, lon <= 250, lat >= 20, lat <= 66)
  g <- d %>%
    mutate(lon = round(lon / 0.5) * 0.5, lat = round(lat / 0.5) * 0.5) %>%
    group_by(lon, lat) %>%
    summarise(beta = mean(beta, na.rm = TRUE),
              qval = mean(qval, na.rm = TRUE), .groups = "drop") %>%
    complete(lon = seq(150, 250, by = 0.5), lat = seq(20, 66, by = 0.5))
  lim <- max(abs(d$beta), na.rm = TRUE)

  ggplot() +
    geom_raster(data = filter(g, !is.na(beta)), aes(x = lon, y = lat, fill = beta)) +
    geom_contour(data = filter(g, !is.na(qval)), aes(x = lon, y = lat, z = qval),
                 breaks = 0.05, color = "black", linewidth = 0.4) +
    geom_polygon(data = mapWorld.clean, aes(x = long, y = lat, group = group),
                 fill = "gray85", color = "gray30", linewidth = 0.25) +
    geom_path(data = al.box.df, aes(x = x, y = y),
              color = "black", linewidth = 0.7, linetype = "dashed") +
    scale_fill_gradient2(low = "steelblue4", mid = "white", high = "darkorange4",
                         midpoint = 0, limits = c(-lim, lim), name = fill.lab) +
    coord_map(projection = "rectangular", parameters = 55,
              xlim = c(150, 250), ylim = c(20, 66)) +
    scale_x_continuous(breaks = seq(150, 250, 20), labels = lon_label) +
    scale_y_continuous(breaks = seq(20, 60, 10),
                       labels = function(y) paste0(y, "°N")) +
    labs(title = ttl, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "gray90", linewidth = 0.3))
}

p.sst.nep <- nep_reg_panel(sst.reg.df, "β (°C / hPa)",
                           "Winter SST regressed on winter AL SLP")
p.mld.nep <- nep_reg_panel(mld.reg.df, "β (m / hPa)",
                           "Winter MLD regressed on winter AL SLP")

p.sst.mld.nep <- (p.sst.nep / p.mld.nep) +
  plot_annotation(
    title    = "NE Pacific (east of 150°E): SST (top) and MLD (bottom) ~ winter AL SLP",
    subtitle = "GLS AR(1); contour = Benjamini-Hochberg FDR q ≤ 0.05 (Wilks 2016)",
    theme    = theme(plot.title    = element_text(size = 12, hjust = 0.5),
                     plot.subtitle = element_text(size = 9,  hjust = 0.5))
  )

ggsave("./Figures/SST_MLD_AL_SLP_regression_NEP.png",
       plot = p.sst.mld.nep, width = 8, height = 9, dpi = 300)
message("Saved: Figures/SST_MLD_AL_SLP_regression_NEP.png")

# ============================================================
# SECTION 7: 15-YR ROLLING AL SLP SD vs. MLD AR(1) — dual axis
# ============================================================
# Mirrors Section 5 (SST AR(1) dual-axis) but uses winter-mean MLD
# anomaly (detrended) for GOA and EBS. Loads pre-extracted ORAS5
# MLD box NetCDFs, masks to the GOA/EBS polygons (goa.sf / ebs.sf
# from Section 2), computes area-weighted monthly series, detrends
# monthly anomalies, then rolls 15-yr right-aligned AR(1). AL SLP
# SD is the same 15-yr rolling z-scored series from Section 5.

mld.ts.cache <- "./Output/oras5_mld_GOA_EBS_winter_anom.csv"

if (file.exists(mld.ts.cache)) {
  message("Loading cached GOA/EBS winter MLD anomaly from: ", mld.ts.cache)
  mld.win <- read.csv(mld.ts.cache)
} else {
  load_mld_box <- function(file) {
    nc   <- nc_open(file)
    m    <- ncvar_get(nc, "somxl030")
    lons <- ncvar_get(nc, "nav_lon")
    lats <- ncvar_get(nc, "nav_lat")
    tv   <- ncvar_get(nc, "time_counter")
    tu   <- ncatt_get(nc, "time_counter", "units")$value
    nc_close(nc)
    m[m > 1e10] <- NA
    origin <- sub(".*since ", "", tu)
    dates  <- as.Date(as.POSIXct(tv, origin = origin, tz = "UTC"))
    list(mld = m, lons = lons, lats = lats, dates = dates)
  }

  region_ts_mld <- function(d, poly.sf) {
    lons.v <- as.vector(d$lons); lats.v <- as.vector(d$lats)
    pts <- st_as_sf(data.frame(lon = lons.v, lat = lats.v),
                    coords = c("lon", "lat"), crs = 4326)
    in.poly <- lengths(st_within(pts, poly.sf)) > 0
    if (!any(in.poly)) stop("no MLD cells inside polygon")
    w  <- cos(lats.v[in.poly] * pi / 180)
    nt <- dim(d$mld)[3]
    mat <- matrix(d$mld, nrow = prod(dim(d$mld)[1:2]), ncol = nt)
    sapply(seq_len(nt), function(t) {
      v <- mat[in.poly, t]
      if (all(is.na(v))) NA_real_ else weighted.mean(v, w = w, na.rm = TRUE)
    })
  }

  message("Loading GOA / EBS MLD boxes and masking to polygons...")
  goa.d <- load_mld_box("./Data/oras5_mld_GOA_box.nc")
  ebs.d <- load_mld_box("./Data/oras5_mld_EBS_box.nc")

  goa.mld.ts <- region_ts_mld(goa.d, goa.sf)
  ebs.mld.ts <- region_ts_mld(ebs.d, ebs.sf)

  dates.mld  <- goa.d$dates
  years.mld  <- as.integer(format(dates.mld, "%Y"))
  months.mld <- as.integer(format(dates.mld, "%m"))

  monthly_anom <- function(ts) {
    df <- data.frame(year = years.mld, month = months.mld, mld = ts)
    clim <- df %>% group_by(month) %>%
      summarise(clim = mean(mld, na.rm = TRUE), .groups = "drop")
    df %>% left_join(clim, by = "month") %>%
      mutate(anom      = mld - clim,
             anom.detr = residuals(lm(anom ~ seq_along(anom))))
  }

  winter_from_monthly <- function(df) {
    df %>%
      filter(month %in% c(11, 12, 1, 2, 3)) %>%
      mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
      group_by(win.year) %>%
      summarise(anom.detr = mean(anom.detr, na.rm = TRUE), .groups = "drop") %>%
      rename(year = win.year) %>%
      arrange(year)
  }

  goa.mld.win <- winter_from_monthly(monthly_anom(goa.mld.ts)) %>%
    rename(GOA = anom.detr)
  ebs.mld.win <- winter_from_monthly(monthly_anom(ebs.mld.ts)) %>%
    rename(EBS = anom.detr)
  mld.win <- inner_join(goa.mld.win, ebs.mld.win, by = "year")

  write.csv(mld.win, mld.ts.cache, row.names = FALSE)
  message("Saved: ", mld.ts.cache)
}

mld.roll <- mld.win %>%
  arrange(year) %>%
  mutate(GOA.ar1 = roll_ar1(GOA),
         EBS.ar1 = roll_ar1(EBS))

al.mld <- mld.roll %>%
  select(year, GOA.ar1, EBS.ar1) %>%
  left_join(al %>% select(year, AL.sd), by = "year") %>%
  filter(!is.na(AL.sd))

mld.al.sd.range   <- range(al.mld$AL.sd,   na.rm = TRUE)
goa.mld.ar1.range <- range(al.mld$GOA.ar1, na.rm = TRUE)
ebs.mld.ar1.range <- range(al.mld$EBS.ar1, na.rm = TRUE)

dat.goa.mld <- al.mld %>% filter(!is.na(GOA.ar1), !is.na(AL.sd))
dat.ebs.mld <- al.mld %>% filter(!is.na(EBS.ar1), !is.na(AL.sd))

chel.goa.mld <- chelton_test(dat.goa.mld$GOA.ar1, dat.goa.mld$AL.sd)
chel.ebs.mld <- chelton_test(dat.ebs.mld$EBS.ar1, dat.ebs.mld$AL.sd)

message(sprintf("MLD GOA Chelton: r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
                chel.goa.mld$r.obs, chel.goa.mld$n.eff,
                chel.goa.mld$t.obs, chel.goa.mld$pval))
message(sprintf("MLD EBS Chelton: r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
                chel.ebs.mld$r.obs, chel.ebs.mld$n.eff,
                chel.ebs.mld$t.obs, chel.ebs.mld$pval))

make_mld_dual_axis_plot <- function(dat, ar1.col, ar1.color, region.label,
                                    ar1.range, chel) {
  k <- diff(ar1.range) / diff(mld.al.sd.range)
  b <- ar1.range[1] - k * mld.al.sd.range[1]

  dat <- dat %>% filter(!is.na(.data[[ar1.col]]), !is.na(AL.sd))

  lbl <- sprintf("r = %.2f, N* = %.1f\nt = %.2f, p = %.3f",
                 chel$r.obs, chel$n.eff, chel$t.obs, chel$pval)

  ggplot(dat, aes(x = year)) +
    geom_line(aes(y = .data[[ar1.col]], color = "MLD AR(1)"), linewidth = 0.7) +
    geom_point(aes(y = .data[[ar1.col]], color = "MLD AR(1)"), size = 1.5) +
    geom_line(aes(y = k * AL.sd + b, color = "AL SLP SD"),
              linewidth = 0.7, linetype = "dashed") +
    geom_point(aes(y = k * AL.sd + b, color = "AL SLP SD"), size = 1.5, shape = 1) +
    annotate("text", x = Inf, y = -Inf, label = lbl,
             hjust = 1.1, vjust = -0.4, size = 3.5, lineheight = 1.3) +
    scale_y_continuous(
      name     = "MLD AR(1)",
      sec.axis = sec_axis(~ (. - b) / k, name = "AL SLP SD (z-scored)")
    ) +
    scale_color_manual(values = c("MLD AR(1)" = ar1.color, "AL SLP SD" = "gray30")) +
    labs(title = region.label, x = "Year", color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position    = "bottom",
          axis.title.y.left  = element_text(color = ar1.color),
          axis.text.y.left   = element_text(color = ar1.color),
          axis.title.y.right = element_text(color = "gray30"),
          axis.text.y.right  = element_text(color = "gray30"))
}

p.goa.mld.dual <- make_mld_dual_axis_plot(al.mld, "GOA.ar1", "darkorange4", "GOA",
                                          goa.mld.ar1.range, chel.goa.mld)
p.ebs.mld.dual <- make_mld_dual_axis_plot(al.mld, "EBS.ar1", "steelblue4",  "EBS",
                                          ebs.mld.ar1.range, chel.ebs.mld)

p.mld.dual <- p.goa.mld.dual / p.ebs.mld.dual +
  plot_annotation(
    title = "15-yr Rolling AL SLP SD vs. MLD AR(1) \u2014 Winter (Nov\u2013Mar)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 12))
  )

ggsave("./Figures/AL_SD_MLD_AR1_15yr_dual_axis.png",
       plot = p.mld.dual, width = 7, height = 7, dpi = 300)
message("Saved: Figures/AL_SD_MLD_AR1_15yr_dual_axis.png")

# ============================================================
# SECTION 8: 15-YR ROLLING AL SLP SD vs. MLD ANOMALY — dual axis
# ============================================================
# Same style as Section 7, but the response is the 15-yr right-
# aligned rolling MEAN of winter detrended MLD anomaly (instead of
# AR(1)). AL SLP SD is the same 15-yr rolling z-score from Section 5.
# Inputs reuse mld.win (winter detrended MLD anomaly) built above.

roll_mean_15 <- function(x, width = 15) {
  rollapply(x, width = width, fill = NA, align = "right",
            FUN = mean, na.rm = TRUE)
}

mld.anom.roll <- mld.win %>%
  arrange(year) %>%
  mutate(GOA.anom = roll_mean_15(GOA),
         EBS.anom = roll_mean_15(EBS))

al.mld.anom <- mld.anom.roll %>%
  select(year, GOA.anom, EBS.anom) %>%
  left_join(al %>% select(year, AL.sd), by = "year") %>%
  filter(!is.na(AL.sd))

mld.anom.al.sd.range <- range(al.mld.anom$AL.sd,    na.rm = TRUE)
goa.mld.anom.range   <- range(al.mld.anom$GOA.anom, na.rm = TRUE)
ebs.mld.anom.range   <- range(al.mld.anom$EBS.anom, na.rm = TRUE)

dat.goa.manom <- al.mld.anom %>% filter(!is.na(GOA.anom), !is.na(AL.sd))
dat.ebs.manom <- al.mld.anom %>% filter(!is.na(EBS.anom), !is.na(AL.sd))

chel.goa.manom <- chelton_test(dat.goa.manom$GOA.anom, dat.goa.manom$AL.sd)
chel.ebs.manom <- chelton_test(dat.ebs.manom$EBS.anom, dat.ebs.manom$AL.sd)

message(sprintf("MLD anom GOA Chelton: r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
                chel.goa.manom$r.obs, chel.goa.manom$n.eff,
                chel.goa.manom$t.obs, chel.goa.manom$pval))
message(sprintf("MLD anom EBS Chelton: r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
                chel.ebs.manom$r.obs, chel.ebs.manom$n.eff,
                chel.ebs.manom$t.obs, chel.ebs.manom$pval))

make_mld_anom_dual_plot <- function(dat, anom.col, anom.color, region.label,
                                    anom.range, chel) {
  k <- diff(anom.range) / diff(mld.anom.al.sd.range)
  b <- anom.range[1] - k * mld.anom.al.sd.range[1]

  dat <- dat %>% filter(!is.na(.data[[anom.col]]), !is.na(AL.sd))

  lbl <- sprintf("r = %.2f, N* = %.1f\nt = %.2f, p = %.3f",
                 chel$r.obs, chel$n.eff, chel$t.obs, chel$pval)

  ggplot(dat, aes(x = year)) +
    geom_line(aes(y = .data[[anom.col]], color = "MLD anomaly"), linewidth = 0.7) +
    geom_point(aes(y = .data[[anom.col]], color = "MLD anomaly"), size = 1.5) +
    geom_line(aes(y = k * AL.sd + b, color = "AL SLP SD"),
              linewidth = 0.7, linetype = "dashed") +
    geom_point(aes(y = k * AL.sd + b, color = "AL SLP SD"), size = 1.5, shape = 1) +
    annotate("text", x = Inf, y = -Inf, label = lbl,
             hjust = 1.1, vjust = -0.4, size = 3.5, lineheight = 1.3) +
    scale_y_continuous(
      name     = "MLD anomaly (m)",
      sec.axis = sec_axis(~ (. - b) / k, name = "AL SLP SD (z-scored)")
    ) +
    scale_color_manual(values = c("MLD anomaly" = anom.color,
                                  "AL SLP SD"    = "gray30")) +
    labs(title = region.label, x = "Year", color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position    = "bottom",
          axis.title.y.left  = element_text(color = anom.color),
          axis.text.y.left   = element_text(color = anom.color),
          axis.title.y.right = element_text(color = "gray30"),
          axis.text.y.right  = element_text(color = "gray30"))
}

p.goa.manom.dual <- make_mld_anom_dual_plot(al.mld.anom, "GOA.anom", "darkorange4",
                                            "GOA", goa.mld.anom.range, chel.goa.manom)
p.ebs.manom.dual <- make_mld_anom_dual_plot(al.mld.anom, "EBS.anom", "steelblue4",
                                            "EBS", ebs.mld.anom.range, chel.ebs.manom)

p.manom.dual <- p.goa.manom.dual / p.ebs.manom.dual +
  plot_annotation(
    title = "15-yr Rolling AL SLP SD vs. MLD Anomaly \u2014 Winter (Nov\u2013Mar)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 12))
  )

ggsave("./Figures/AL_SD_MLD_anom_15yr_dual_axis.png",
       plot = p.manom.dual, width = 7, height = 7, dpi = 300)
message("Saved: Figures/AL_SD_MLD_anom_15yr_dual_axis.png")

# ============================================================
# SECTION 9: Cellwise regression of MLD AR(1) on AL SLP SD
# ============================================================
# Same style as Section 6 MLD ~ AL SLP regression map, but both
# variables are 15-yr right-aligned rolling windows of winter values.
# x variable: al$AL.sd  (rolling SD of winter AL SLP anomaly, z-scored)
# y variable: cellwise 15-yr rolling AR(1) of winter-mean detrended
#             MLD anomaly. GLS with corAR1(), BH-FDR contour at q<=0.05.

mld.ar1.reg.cache <- "./Output/mld_ar1_al_slp_sd_regression.csv"

if (file.exists(mld.ar1.reg.cache)) {
  message("Loading cached MLD AR(1) ~ AL SLP SD regression from: ",
          mld.ar1.reg.cache)
  mld.ar1.reg.df <- read.csv(mld.ar1.reg.cache)
} else {
  message("Loading detrended winter-month MLD anomalies (full NP domain)...")
  nc.m    <- nc_open(mld.winmon.file)
  mld.w   <- ncvar_get(nc.m, "somxl030")
  lons.m  <- ncvar_get(nc.m, "nav_lon")
  lats.m  <- ncvar_get(nc.m, "nav_lat")
  time.m  <- ncvar_get(nc.m, "time_counter")
  tun.m   <- ncatt_get(nc.m, "time_counter", "units")$value
  nc_close(nc.m)

  mld.w[mld.w > 1e10] <- NA
  origin.m  <- sub(".*since ", "", tun.m)
  dates.m   <- as.Date(as.POSIXct(time.m, origin = origin.m, tz = "UTC"))
  years.m   <- as.integer(format(dates.m, "%Y"))
  months.m  <- as.integer(format(dates.m, "%m"))
  winyrs.m  <- ifelse(months.m %in% c(11, 12), years.m + 1L, years.m)
  lons.m360 <- ifelse(lons.m < 0, lons.m + 360, lons.m)

  nx <- dim(mld.w)[1]; ny <- dim(mld.w)[2]; nt <- dim(mld.w)[3]
  mat <- matrix(mld.w, nrow = nx * ny, ncol = nt)
  rm(mld.w); gc()

  uy <- sort(unique(winyrs.m))
  win.mld <- matrix(NA_real_, nrow = nx * ny, ncol = length(uy))
  for (k in seq_along(uy)) {
    ti <- which(winyrs.m == uy[k])
    win.mld[, k] <- rowMeans(mat[, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(mat); gc()

  safe_ar1 <- function(v) {
    if (any(is.na(v))) return(NA_real_)
    suppressWarnings(acf(v, lag.max = 1, plot = FALSE)$acf[2])
  }
  roll_ar1_cell <- function(x, width = 15) {
    rollapply(x, width = width, fill = NA, align = "right", FUN = safe_ar1)
  }

  # AL SLP SD from Section 5 (already 15-yr rolling, z-scored)
  al.df   <- al %>% select(year, AL.sd) %>% filter(!is.na(AL.sd))
  yr.use  <- intersect(uy, al.df$year)
  x.vec   <- al.df$AL.sd[match(yr.use, al.df$year)]
  col.idx <- match(yr.use, uy)
  n.yrs   <- length(yr.use)
  message(sprintf("Regression years: %d-%d (%d years)",
                  min(yr.use), max(yr.use), n.yrs))

  # Cells with enough winter data (need a long-enough unbroken series for rolling AR(1))
  good.cells <- which(rowSums(!is.na(win.mld)) >= max(25, 0.5 * length(uy)))
  message(sprintf("Computing cellwise 15-yr rolling AR(1) at %d cells...",
                  length(good.cells)))

  # Per-cell AR(1) series aligned on uy, then subset to yr.use
  ar1.mat <- matrix(NA_real_, nrow = length(good.cells), ncol = length(yr.use))
  for (k in seq_along(good.cells)) {
    a <- roll_ar1_cell(win.mld[good.cells[k], ])
    ar1.mat[k, ] <- a[col.idx]
  }
  rm(win.mld); gc()

  message(sprintf("Fitting GLS AR(1) on %d-cell rolling AR(1) series vs AL SLP SD...",
                  length(good.cells)))

  beta.vec <- rep(NA_real_, length(good.cells))
  pval.vec <- rep(NA_real_, length(good.cells))
  df.reg   <- data.frame(y = NA_real_, x = x.vec)

  for (k in seq_along(good.cells)) {
    y <- ar1.mat[k, ]
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

  beta.full <- rep(NA_real_, nx * ny)
  pval.full <- rep(NA_real_, nx * ny)
  beta.full[good.cells] <- beta.vec
  pval.full[good.cells] <- pval.vec

  mld.ar1.reg.df <- data.frame(
    lon  = as.vector(lons.m360),
    lat  = as.vector(lats.m),
    beta = beta.full,
    pval = pval.full
  ) %>% filter(!is.na(beta))

  write.csv(mld.ar1.reg.df, mld.ar1.reg.cache, row.names = FALSE)
  message("Saved: ", mld.ar1.reg.cache)
}

mld.ar1.reg.df <- mld.ar1.reg.df %>%
  filter(lon >= 160, lon <= 250, lat >= 20, lat <= 66)

mld.ar1.reg.df$qval <- p.adjust(mld.ar1.reg.df$pval, method = "BH")

n.sig.raw <- sum(mld.ar1.reg.df$pval <= 0.05, na.rm = TRUE)
n.sig.fdr <- sum(mld.ar1.reg.df$qval <= 0.05, na.rm = TRUE)
message(sprintf("Raw p<=0.05: %d cells | BH-FDR q<=0.05: %d cells (of %d)",
                n.sig.raw, n.sig.fdr, nrow(mld.ar1.reg.df)))

mld.ar1.reg.grid <- mld.ar1.reg.df %>%
  mutate(lon = round(lon / 0.5) * 0.5,
         lat = round(lat / 0.5) * 0.5) %>%
  group_by(lon, lat) %>%
  summarise(beta = mean(beta, na.rm = TRUE),
            qval = mean(qval, na.rm = TRUE),
            .groups = "drop") %>%
  complete(lon = seq(160, 250, by = 0.5),
           lat = seq(20,  66,  by = 0.5))

col.lim.ar1 <- max(abs(mld.ar1.reg.df$beta), na.rm = TRUE)

p.mld.ar1.reg <- ggplot() +
  geom_raster(data = filter(mld.ar1.reg.grid, !is.na(beta)),
              aes(x = lon, y = lat, fill = beta)) +
  geom_contour(data = filter(mld.ar1.reg.grid, !is.na(qval)),
               aes(x = lon, y = lat, z = qval),
               breaks = 0.05, color = "black", linewidth = 0.4) +
  geom_polygon(data = mapWorld.clean,
               aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray30", linewidth = 0.25) +
  geom_path(data = al.box.df, aes(x = x, y = y),
            color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue4", mid = "white",
                       high = "darkorange4", midpoint = 0,
                       limits = c(-col.lim.ar1, col.lim.ar1),
                       name = "\u03b2 (AR(1) / z)") +
  # Revision: trim western boundary to 160E. Use coord_cartesian
  # (NOT coord_map) so world polygons don't get mangled and overdraw
  # the data area; matches Section 17 SST AR(1) treatment.
  coord_cartesian(xlim = c(160, 250), ylim = c(20, 66)) +
  scale_x_continuous(breaks = seq(160, 250, 20), labels = lon_label) +
  scale_y_continuous(breaks = seq(20, 60, 10),
                     labels = function(y) paste0(y, "\u00b0N")) +
  labs(title    = "Winter MLD AR(1) regressed on AL SLP SD (15-yr rolling windows)",
       subtitle = "GLS AR(1); contour = Benjamini-Hochberg FDR q \u2264 0.05 (Wilks 2016)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/MLD_AR1_AL_SLP_SD_regression.png",
       plot = p.mld.ar1.reg, width = 9, height = 6, dpi = 300)
message("Saved: Figures/MLD_AR1_AL_SLP_SD_regression.png")

# ============================================================
# SECTION 10: 1-DEG AGGREGATED MLD AR(1) regressed on AL SLP SD
# ============================================================
# Same workflow as Section 9, but raw MLD is first averaged into
# 1 deg x 1 deg spatial blocks BEFORE computing climatology,
# anomaly, detrending, winter means, rolling AR(1), and GLS fit.

mld.1deg.reg.cache <- "./Output/mld_ar1_al_slp_sd_regression_1deg.csv"
raw.mld.file       <- "./Data/oras5_mld_NP_1958_2026.nc"

if (file.exists(mld.1deg.reg.cache)) {
  message("Loading cached 1-deg MLD AR(1) ~ AL SLP SD regression from: ",
          mld.1deg.reg.cache)
  mld.1deg.reg.df <- read.csv(mld.1deg.reg.cache)
} else {
  message("Loading raw MLD NC (streaming, aggregating to 1-deg) ...")
  nc.r   <- nc_open(raw.mld.file)
  lons.r <- ncvar_get(nc.r, "nav_lon")
  lats.r <- ncvar_get(nc.r, "nav_lat")
  time.r <- ncvar_get(nc.r, "time_counter")
  tun.r  <- ncatt_get(nc.r, "time_counter", "units")$value
  origin.r <- sub(".*since ", "", tun.r)
  dates.r  <- as.Date(as.POSIXct(time.r, origin = origin.r, tz = "UTC"))
  years.r  <- as.integer(format(dates.r, "%Y"))
  months.r <- as.integer(format(dates.r, "%m"))

  lons.r360 <- ifelse(as.vector(lons.r) < 0, as.vector(lons.r) + 360,
                      as.vector(lons.r))
  lats.v    <- as.vector(lats.r)

  # Restrict to NP domain (20-66N, 110-250E) and build 1-deg bin IDs
  in.domain <- lons.r360 >= 110 & lons.r360 <= 250 &
               lats.v    >= 20  & lats.v    <= 66
  bin.lon <- floor(lons.r360) + 0.5
  bin.lat <- floor(lats.v)    + 0.5
  bin.key <- ifelse(in.domain, paste(bin.lon, bin.lat, sep = "_"), NA)

  uniq.keys <- sort(unique(stats::na.omit(bin.key)))
  bin.idx   <- match(bin.key, uniq.keys)   # NA for out-of-domain cells
  nbin      <- length(uniq.keys)

  nx <- dim(lons.r)[1]; ny <- dim(lons.r)[2]; nt <- length(time.r)
  bin.coords <- do.call(rbind, strsplit(uniq.keys, "_"))
  bin.lon.v  <- as.numeric(bin.coords[, 1])
  bin.lat.v  <- as.numeric(bin.coords[, 2])

  # Stream in chunks to limit memory; accumulate 1-deg means (unweighted)
  chunk.size <- 120
  starts     <- seq(1, nt, by = chunk.size)
  mat.bin    <- matrix(NA_real_, nrow = nbin, ncol = nt)

  for (s in starts) {
    n.read <- min(chunk.size, nt - s + 1)
    cat(sprintf("  chunk %d-%d / %d\n", s, s + n.read - 1, nt))
    chunk <- ncvar_get(nc.r, "somxl030",
                       start = c(1, 1, s), count = c(-1, -1, n.read))
    chunk[chunk > 1e10] <- NA
    chunk.mat <- matrix(chunk, nrow = nx * ny, ncol = n.read)
    rm(chunk)
    # Aggregate each timestep by 1-deg bin
    valid <- !is.na(bin.idx)
    for (tt in seq_len(n.read)) {
      v  <- chunk.mat[, tt]
      ok <- valid & !is.na(v)
      if (!any(ok)) next
      sums   <- tapply(v[ok], bin.idx[ok], sum)
      counts <- tapply(v[ok], bin.idx[ok], length)
      keys.t <- as.integer(names(sums))
      mat.bin[keys.t, s + tt - 1] <- sums / counts
    }
    rm(chunk.mat); gc()
  }
  nc_close(nc.r)

  # Climatology + monthly anomaly + detrend per 1-deg cell
  message("Computing climatology, anomalies, detrending, winter means per 1-deg cell...")
  clim.mat <- matrix(NA_real_, nrow = nbin, ncol = 12)
  for (m in 1:12) {
    ti <- which(months.r == m)
    clim.mat[, m] <- rowMeans(mat.bin[, ti, drop = FALSE], na.rm = TRUE)
  }
  anom.mat <- mat.bin - clim.mat[, months.r]
  rm(mat.bin); gc()

  # Detrend each row (NA-safe)
  tvec <- seq_len(nt)
  for (k in seq_len(nbin)) {
    y <- anom.mat[k, ]
    ok <- !is.na(y)
    if (sum(ok) < 24) { anom.mat[k, ] <- NA_real_; next }
    f  <- lm(y ~ tvec, subset = ok)
    anom.mat[k, ok] <- residuals(f)
  }

  # Winter annual means per 1-deg cell
  winyrs.r <- ifelse(months.r %in% c(11, 12), years.r + 1L, years.r)
  in.win   <- months.r %in% c(11, 12, 1, 2, 3)
  uy       <- sort(unique(winyrs.r[in.win]))
  win.mld.1 <- matrix(NA_real_, nrow = nbin, ncol = length(uy))
  for (k in seq_along(uy)) {
    ti <- which(winyrs.r == uy[k] & in.win)
    if (length(ti) == 0) next
    win.mld.1[, k] <- rowMeans(anom.mat[, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(anom.mat); gc()

  safe_ar1 <- function(v) {
    if (any(is.na(v))) return(NA_real_)
    suppressWarnings(acf(v, lag.max = 1, plot = FALSE)$acf[2])
  }
  roll_ar1_cell <- function(x, width = 15) {
    rollapply(x, width = width, fill = NA, align = "right", FUN = safe_ar1)
  }

  al.df   <- al %>% select(year, AL.sd) %>% filter(!is.na(AL.sd))
  yr.use  <- intersect(uy, al.df$year)
  x.vec   <- al.df$AL.sd[match(yr.use, al.df$year)]
  col.idx <- match(yr.use, uy)

  good.cells <- which(rowSums(!is.na(win.mld.1)) >= max(25, 0.5 * length(uy)))
  message(sprintf("Computing 15-yr rolling AR(1) + GLS at %d one-degree bins...",
                  length(good.cells)))

  ar1.mat <- matrix(NA_real_, nrow = length(good.cells), ncol = length(yr.use))
  for (k in seq_along(good.cells)) {
    a <- roll_ar1_cell(win.mld.1[good.cells[k], ])
    ar1.mat[k, ] <- a[col.idx]
  }
  rm(win.mld.1); gc()

  beta.vec <- rep(NA_real_, length(good.cells))
  pval.vec <- rep(NA_real_, length(good.cells))
  df.reg   <- data.frame(y = NA_real_, x = x.vec)
  for (k in seq_along(good.cells)) {
    y <- ar1.mat[k, ]
    if (sum(!is.na(y)) < 10) next
    df.reg$y <- y
    fit <- tryCatch(
      gls(y ~ x, data = df.reg, correlation = corAR1(form = ~ 1),
          method = "ML", na.action = na.omit),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    s <- summary(fit)$tTable
    if (nrow(s) < 2) next
    beta.vec[k] <- s[2, 1]
    pval.vec[k] <- s[2, 4]
  }

  beta.full <- rep(NA_real_, nbin); pval.full <- rep(NA_real_, nbin)
  beta.full[good.cells] <- beta.vec
  pval.full[good.cells] <- pval.vec

  mld.1deg.reg.df <- data.frame(
    lon  = bin.lon.v,
    lat  = bin.lat.v,
    beta = beta.full,
    pval = pval.full
  ) %>% filter(!is.na(beta))

  write.csv(mld.1deg.reg.df, mld.1deg.reg.cache, row.names = FALSE)
  message("Saved: ", mld.1deg.reg.cache)
}

mld.1deg.reg.df <- mld.1deg.reg.df %>%
  filter(lon >= 110, lon <= 250, lat >= 20, lat <= 66)
mld.1deg.reg.df$qval <- p.adjust(mld.1deg.reg.df$pval, method = "BH")

n.sig.raw <- sum(mld.1deg.reg.df$pval <= 0.05, na.rm = TRUE)
n.sig.fdr <- sum(mld.1deg.reg.df$qval <= 0.10, na.rm = TRUE)
message(sprintf("1-deg regression: raw p<=0.05: %d | BH-FDR q<=0.10: %d (of %d)",
                n.sig.raw, n.sig.fdr, nrow(mld.1deg.reg.df)))

mld.1deg.grid <- mld.1deg.reg.df %>%
  complete(lon = seq(110.5, 249.5, by = 1),
           lat = seq(20.5,  65.5,  by = 1))

col.lim.1deg <- max(abs(mld.1deg.reg.df$beta), na.rm = TRUE)

p.mld.1deg <- ggplot() +
  geom_raster(data = filter(mld.1deg.grid, !is.na(beta)),
              aes(x = lon, y = lat, fill = beta)) +
  geom_contour(data = filter(mld.1deg.grid, !is.na(qval)),
               aes(x = lon, y = lat, z = qval),
               breaks = 0.10, color = "black", linewidth = 0.4) +
  geom_polygon(data = mapWorld.clean,
               aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray30", linewidth = 0.25) +
  geom_path(data = al.box.df, aes(x = x, y = y),
            color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue4", mid = "white",
                       high = "darkorange4", midpoint = 0,
                       limits = c(-col.lim.1deg, col.lim.1deg),
                       name = "\u03b2 (AR(1) / z)") +
  # Revision: trim western boundary to 160E. Use coord_cartesian
  # (NOT coord_map) so world polygons don't get mangled and overdraw
  # the data area; matches Section 17 SST AR(1) treatment.
  coord_cartesian(xlim = c(160, 250), ylim = c(20, 66)) +
  scale_x_continuous(breaks = seq(160, 250, 20), labels = lon_label) +
  scale_y_continuous(breaks = seq(20, 60, 10),
                     labels = function(y) paste0(y, "\u00b0N")) +
  labs(title    = "1\u00b0-aggregated winter MLD AR(1) regressed on AL SLP SD (15-yr rolling)",
       subtitle = "GLS AR(1); contour = Benjamini-Hochberg FDR q \u2264 0.10",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/MLD_AR1_AL_SLP_SD_regression_1deg.png",
       plot = p.mld.1deg, width = 9, height = 6, dpi = 300)
message("Saved: Figures/MLD_AR1_AL_SLP_SD_regression_1deg.png")

# ============================================================
# SECTION 11: 2-DEG AGGREGATED MLD AR(1) regressed on AL SLP SD
# ============================================================
# Identical pipeline to Section 10 but with MLD averaged into
# 2 deg x 2 deg spatial blocks before anomaly/detrend/AR(1)/GLS.

mld.2deg.reg.cache <- "./Output/mld_ar1_al_slp_sd_regression_2deg.csv"

if (file.exists(mld.2deg.reg.cache)) {
  message("Loading cached 2-deg MLD AR(1) ~ AL SLP SD regression from: ",
          mld.2deg.reg.cache)
  mld.2deg.reg.df <- read.csv(mld.2deg.reg.cache)
} else {
  message("Loading raw MLD NC (streaming, aggregating to 2-deg) ...")
  nc.r   <- nc_open(raw.mld.file)
  lons.r <- ncvar_get(nc.r, "nav_lon")
  lats.r <- ncvar_get(nc.r, "nav_lat")
  time.r <- ncvar_get(nc.r, "time_counter")
  tun.r  <- ncatt_get(nc.r, "time_counter", "units")$value
  origin.r <- sub(".*since ", "", tun.r)
  dates.r  <- as.Date(as.POSIXct(time.r, origin = origin.r, tz = "UTC"))
  years.r  <- as.integer(format(dates.r, "%Y"))
  months.r <- as.integer(format(dates.r, "%m"))

  lons.r360 <- ifelse(as.vector(lons.r) < 0, as.vector(lons.r) + 360,
                      as.vector(lons.r))
  lats.v    <- as.vector(lats.r)

  in.domain <- lons.r360 >= 110 & lons.r360 <= 250 &
               lats.v    >= 20  & lats.v    <= 66
  # 2-deg bin centers: 111, 113, ..., 249 (lon) and 21, 23, ..., 65 (lat)
  bin.lon <- floor(lons.r360 / 2) * 2 + 1
  bin.lat <- floor(lats.v    / 2) * 2 + 1
  bin.key <- ifelse(in.domain, paste(bin.lon, bin.lat, sep = "_"), NA)

  uniq.keys <- sort(unique(stats::na.omit(bin.key)))
  bin.idx   <- match(bin.key, uniq.keys)
  nbin      <- length(uniq.keys)

  nx <- dim(lons.r)[1]; ny <- dim(lons.r)[2]; nt <- length(time.r)
  bin.coords <- do.call(rbind, strsplit(uniq.keys, "_"))
  bin.lon.v  <- as.numeric(bin.coords[, 1])
  bin.lat.v  <- as.numeric(bin.coords[, 2])

  chunk.size <- 120
  starts     <- seq(1, nt, by = chunk.size)
  mat.bin    <- matrix(NA_real_, nrow = nbin, ncol = nt)

  for (s in starts) {
    n.read <- min(chunk.size, nt - s + 1)
    cat(sprintf("  chunk %d-%d / %d\n", s, s + n.read - 1, nt))
    chunk <- ncvar_get(nc.r, "somxl030",
                       start = c(1, 1, s), count = c(-1, -1, n.read))
    chunk[chunk > 1e10] <- NA
    chunk.mat <- matrix(chunk, nrow = nx * ny, ncol = n.read)
    rm(chunk)
    valid <- !is.na(bin.idx)
    for (tt in seq_len(n.read)) {
      v  <- chunk.mat[, tt]
      ok <- valid & !is.na(v)
      if (!any(ok)) next
      sums   <- tapply(v[ok], bin.idx[ok], sum)
      counts <- tapply(v[ok], bin.idx[ok], length)
      keys.t <- as.integer(names(sums))
      mat.bin[keys.t, s + tt - 1] <- sums / counts
    }
    rm(chunk.mat); gc()
  }
  nc_close(nc.r)

  message("Computing climatology, anomalies, detrending, winter means per 2-deg cell...")
  clim.mat <- matrix(NA_real_, nrow = nbin, ncol = 12)
  for (m in 1:12) {
    ti <- which(months.r == m)
    clim.mat[, m] <- rowMeans(mat.bin[, ti, drop = FALSE], na.rm = TRUE)
  }
  anom.mat <- mat.bin - clim.mat[, months.r]
  rm(mat.bin); gc()

  tvec <- seq_len(nt)
  for (k in seq_len(nbin)) {
    y <- anom.mat[k, ]
    ok <- !is.na(y)
    if (sum(ok) < 24) { anom.mat[k, ] <- NA_real_; next }
    f  <- lm(y ~ tvec, subset = ok)
    anom.mat[k, ok] <- residuals(f)
  }

  winyrs.r <- ifelse(months.r %in% c(11, 12), years.r + 1L, years.r)
  in.win   <- months.r %in% c(11, 12, 1, 2, 3)
  uy       <- sort(unique(winyrs.r[in.win]))
  win.mld.2 <- matrix(NA_real_, nrow = nbin, ncol = length(uy))
  for (k in seq_along(uy)) {
    ti <- which(winyrs.r == uy[k] & in.win)
    if (length(ti) == 0) next
    win.mld.2[, k] <- rowMeans(anom.mat[, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(anom.mat); gc()

  safe_ar1 <- function(v) {
    if (any(is.na(v))) return(NA_real_)
    suppressWarnings(acf(v, lag.max = 1, plot = FALSE)$acf[2])
  }
  roll_ar1_cell <- function(x, width = 15) {
    rollapply(x, width = width, fill = NA, align = "right", FUN = safe_ar1)
  }

  al.df   <- al %>% select(year, AL.sd) %>% filter(!is.na(AL.sd))
  yr.use  <- intersect(uy, al.df$year)
  x.vec   <- al.df$AL.sd[match(yr.use, al.df$year)]
  col.idx <- match(yr.use, uy)

  good.cells <- which(rowSums(!is.na(win.mld.2)) >= max(25, 0.5 * length(uy)))
  message(sprintf("Computing 15-yr rolling AR(1) + GLS at %d two-degree bins...",
                  length(good.cells)))

  ar1.mat <- matrix(NA_real_, nrow = length(good.cells), ncol = length(yr.use))
  for (k in seq_along(good.cells)) {
    a <- roll_ar1_cell(win.mld.2[good.cells[k], ])
    ar1.mat[k, ] <- a[col.idx]
  }
  rm(win.mld.2); gc()

  beta.vec <- rep(NA_real_, length(good.cells))
  pval.vec <- rep(NA_real_, length(good.cells))
  df.reg   <- data.frame(y = NA_real_, x = x.vec)
  for (k in seq_along(good.cells)) {
    y <- ar1.mat[k, ]
    if (sum(!is.na(y)) < 10) next
    df.reg$y <- y
    fit <- tryCatch(
      gls(y ~ x, data = df.reg, correlation = corAR1(form = ~ 1),
          method = "ML", na.action = na.omit),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    s <- summary(fit)$tTable
    if (nrow(s) < 2) next
    beta.vec[k] <- s[2, 1]
    pval.vec[k] <- s[2, 4]
  }

  beta.full <- rep(NA_real_, nbin); pval.full <- rep(NA_real_, nbin)
  beta.full[good.cells] <- beta.vec
  pval.full[good.cells] <- pval.vec

  mld.2deg.reg.df <- data.frame(
    lon  = bin.lon.v,
    lat  = bin.lat.v,
    beta = beta.full,
    pval = pval.full
  ) %>% filter(!is.na(beta))

  write.csv(mld.2deg.reg.df, mld.2deg.reg.cache, row.names = FALSE)
  message("Saved: ", mld.2deg.reg.cache)
}

mld.2deg.reg.df <- mld.2deg.reg.df %>%
  filter(lon >= 110, lon <= 250, lat >= 20, lat <= 66)
mld.2deg.reg.df$qval <- p.adjust(mld.2deg.reg.df$pval, method = "BH")

n.sig.raw <- sum(mld.2deg.reg.df$pval <= 0.05, na.rm = TRUE)
n.sig.fdr <- sum(mld.2deg.reg.df$qval <= 0.10, na.rm = TRUE)
message(sprintf("2-deg regression: raw p<=0.05: %d | BH-FDR q<=0.10: %d (of %d)",
                n.sig.raw, n.sig.fdr, nrow(mld.2deg.reg.df)))

mld.2deg.grid <- mld.2deg.reg.df %>%
  complete(lon = seq(111, 249, by = 2),
           lat = seq(21,  65,  by = 2))

col.lim.2deg <- max(abs(mld.2deg.reg.df$beta), na.rm = TRUE)

p.mld.2deg <- ggplot() +
  geom_raster(data = filter(mld.2deg.grid, !is.na(beta)),
              aes(x = lon, y = lat, fill = beta)) +
  geom_contour(data = filter(mld.2deg.grid, !is.na(qval)),
               aes(x = lon, y = lat, z = qval),
               breaks = 0.10, color = "black", linewidth = 0.4) +
  geom_polygon(data = mapWorld.clean,
               aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray30", linewidth = 0.25) +
  geom_path(data = al.box.df, aes(x = x, y = y),
            color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue4", mid = "white",
                       high = "darkorange4", midpoint = 0,
                       limits = c(-col.lim.2deg, col.lim.2deg),
                       name = "\u03b2 (AR(1) / z)") +
  # Revision: trim western boundary to 160E. Use coord_cartesian
  # (NOT coord_map) so world polygons don't get mangled and overdraw
  # the data area; matches Section 17 SST AR(1) treatment.
  coord_cartesian(xlim = c(160, 250), ylim = c(20, 66)) +
  scale_x_continuous(breaks = seq(160, 250, 20), labels = lon_label) +
  scale_y_continuous(breaks = seq(20, 60, 10),
                     labels = function(y) paste0(y, "\u00b0N")) +
  labs(title    = "2\u00b0-aggregated winter MLD AR(1) regressed on AL SLP SD (15-yr rolling)",
       subtitle = "GLS AR(1); contour = Benjamini-Hochberg FDR q \u2264 0.10",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/MLD_AR1_AL_SLP_SD_regression_2deg.png",
       plot = p.mld.2deg, width = 9, height = 6, dpi = 300)
message("Saved: Figures/MLD_AR1_AL_SLP_SD_regression_2deg.png")

# ============================================================
# SECTION 12: Cellwise regression of MLD anomaly on AL SLP SD
# ============================================================
# Same pipeline as Section 9 but the response is the per-cell
# 15-yr right-aligned rolling MEAN of winter detrended MLD
# anomaly, NOT the rolling AR(1).

mld.anom.reg.cache <- "./Output/mld_anom_al_slp_sd_regression.csv"

if (file.exists(mld.anom.reg.cache)) {
  message("Loading cached MLD anomaly ~ AL SLP SD regression from: ",
          mld.anom.reg.cache)
  mld.anom.reg.df <- read.csv(mld.anom.reg.cache)
} else {
  message("Loading detrended winter-month MLD anomalies (full NP domain)...")
  nc.m    <- nc_open(mld.winmon.file)
  mld.w   <- ncvar_get(nc.m, "somxl030")
  lons.m  <- ncvar_get(nc.m, "nav_lon")
  lats.m  <- ncvar_get(nc.m, "nav_lat")
  time.m  <- ncvar_get(nc.m, "time_counter")
  tun.m   <- ncatt_get(nc.m, "time_counter", "units")$value
  nc_close(nc.m)

  mld.w[mld.w > 1e10] <- NA
  origin.m  <- sub(".*since ", "", tun.m)
  dates.m   <- as.Date(as.POSIXct(time.m, origin = origin.m, tz = "UTC"))
  years.m   <- as.integer(format(dates.m, "%Y"))
  months.m  <- as.integer(format(dates.m, "%m"))
  winyrs.m  <- ifelse(months.m %in% c(11, 12), years.m + 1L, years.m)
  lons.m360 <- ifelse(lons.m < 0, lons.m + 360, lons.m)

  nx <- dim(mld.w)[1]; ny <- dim(mld.w)[2]; nt <- dim(mld.w)[3]
  mat <- matrix(mld.w, nrow = nx * ny, ncol = nt)
  rm(mld.w); gc()

  uy <- sort(unique(winyrs.m))
  win.mld <- matrix(NA_real_, nrow = nx * ny, ncol = length(uy))
  for (k in seq_along(uy)) {
    ti <- which(winyrs.m == uy[k])
    win.mld[, k] <- rowMeans(mat[, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(mat); gc()

  safe_mean <- function(v) {
    if (any(is.na(v))) return(NA_real_)
    mean(v)
  }
  roll_mean_cell <- function(x, width = 15) {
    rollapply(x, width = width, fill = NA, align = "right", FUN = safe_mean)
  }

  al.df   <- al %>% select(year, AL.sd) %>% filter(!is.na(AL.sd))
  yr.use  <- intersect(uy, al.df$year)
  x.vec   <- al.df$AL.sd[match(yr.use, al.df$year)]
  col.idx <- match(yr.use, uy)
  message(sprintf("Regression years: %d-%d (%d years)",
                  min(yr.use), max(yr.use), length(yr.use)))

  good.cells <- which(rowSums(!is.na(win.mld)) >= max(25, 0.5 * length(uy)))
  message(sprintf("Computing cellwise 15-yr rolling MEAN at %d cells...",
                  length(good.cells)))

  rmean.mat <- matrix(NA_real_, nrow = length(good.cells), ncol = length(yr.use))
  for (k in seq_along(good.cells)) {
    a <- roll_mean_cell(win.mld[good.cells[k], ])
    rmean.mat[k, ] <- a[col.idx]
  }
  rm(win.mld); gc()

  message(sprintf("Fitting GLS AR(1) on %d-cell rolling MEAN series vs AL SLP SD...",
                  length(good.cells)))

  beta.vec <- rep(NA_real_, length(good.cells))
  pval.vec <- rep(NA_real_, length(good.cells))
  df.reg   <- data.frame(y = NA_real_, x = x.vec)

  for (k in seq_along(good.cells)) {
    y <- rmean.mat[k, ]
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

  beta.full <- rep(NA_real_, nx * ny)
  pval.full <- rep(NA_real_, nx * ny)
  beta.full[good.cells] <- beta.vec
  pval.full[good.cells] <- pval.vec

  mld.anom.reg.df <- data.frame(
    lon  = as.vector(lons.m360),
    lat  = as.vector(lats.m),
    beta = beta.full,
    pval = pval.full
  ) %>% filter(!is.na(beta))

  write.csv(mld.anom.reg.df, mld.anom.reg.cache, row.names = FALSE)
  message("Saved: ", mld.anom.reg.cache)
}

mld.anom.reg.df <- mld.anom.reg.df %>%
  filter(lon >= 160, lon <= 250, lat >= 20, lat <= 66)

mld.anom.reg.df$qval <- p.adjust(mld.anom.reg.df$pval, method = "BH")

n.sig.raw <- sum(mld.anom.reg.df$pval <= 0.05, na.rm = TRUE)
n.sig.fdr <- sum(mld.anom.reg.df$qval <= 0.05, na.rm = TRUE)
message(sprintf("MLD anom: raw p<=0.05: %d cells | BH-FDR q<=0.05: %d cells (of %d)",
                n.sig.raw, n.sig.fdr, nrow(mld.anom.reg.df)))

mld.anom.reg.grid <- mld.anom.reg.df %>%
  mutate(lon = round(lon / 0.5) * 0.5,
         lat = round(lat / 0.5) * 0.5) %>%
  group_by(lon, lat) %>%
  summarise(beta = mean(beta, na.rm = TRUE),
            qval = mean(qval, na.rm = TRUE),
            .groups = "drop") %>%
  complete(lon = seq(160, 250, by = 0.5),
           lat = seq(20,  66,  by = 0.5))

col.lim.anom <- max(abs(mld.anom.reg.df$beta), na.rm = TRUE)

p.mld.anom.reg <- ggplot() +
  geom_raster(data = filter(mld.anom.reg.grid, !is.na(beta)),
              aes(x = lon, y = lat, fill = beta)) +
  geom_contour(data = filter(mld.anom.reg.grid, !is.na(qval)),
               aes(x = lon, y = lat, z = qval),
               breaks = 0.05, color = "black", linewidth = 0.4) +
  geom_polygon(data = mapWorld.clean,
               aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray30", linewidth = 0.25) +
  geom_path(data = al.box.df, aes(x = x, y = y),
            color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue4", mid = "white",
                       high = "darkorange4", midpoint = 0,
                       limits = c(-col.lim.anom, col.lim.anom),
                       name = "\u03b2 (m / z)") +
  coord_cartesian(xlim = c(160, 250), ylim = c(20, 66)) +
  scale_x_continuous(breaks = seq(160, 250, 20), labels = lon_label) +
  scale_y_continuous(breaks = seq(20, 60, 10),
                     labels = function(y) paste0(y, "\u00b0N")) +
  labs(title    = "Winter MLD anomaly regressed on AL SLP SD (15-yr rolling windows)",
       subtitle = "GLS AR(1); contour = Benjamini-Hochberg FDR q \u2264 0.05 (Wilks 2016)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/MLD_anom_AL_SLP_SD_regression.png",
       plot = p.mld.anom.reg, width = 9, height = 6, dpi = 300)
message("Saved: Figures/MLD_anom_AL_SLP_SD_regression.png")

# ============================================================
# SECTION 13: 1-DEG AGGREGATED MLD anomaly regressed on AL SLP SD
# ============================================================
# Same workflow as Section 10 but the response is the per-cell
# 15-yr right-aligned rolling MEAN of winter detrended MLD anomaly.

mld.anom.1deg.cache <- "./Output/mld_anom_al_slp_sd_regression_1deg.csv"

if (file.exists(mld.anom.1deg.cache)) {
  message("Loading cached 1-deg MLD anom ~ AL SLP SD regression from: ",
          mld.anom.1deg.cache)
  mld.anom.1deg.df <- read.csv(mld.anom.1deg.cache)
} else {
  message("Loading raw MLD NC (streaming, aggregating to 1-deg) ...")
  nc.r   <- nc_open(raw.mld.file)
  lons.r <- ncvar_get(nc.r, "nav_lon")
  lats.r <- ncvar_get(nc.r, "nav_lat")
  time.r <- ncvar_get(nc.r, "time_counter")
  tun.r  <- ncatt_get(nc.r, "time_counter", "units")$value
  origin.r <- sub(".*since ", "", tun.r)
  dates.r  <- as.Date(as.POSIXct(time.r, origin = origin.r, tz = "UTC"))
  years.r  <- as.integer(format(dates.r, "%Y"))
  months.r <- as.integer(format(dates.r, "%m"))

  lons.r360 <- ifelse(as.vector(lons.r) < 0, as.vector(lons.r) + 360,
                      as.vector(lons.r))
  lats.v    <- as.vector(lats.r)

  in.domain <- lons.r360 >= 110 & lons.r360 <= 250 &
               lats.v    >= 20  & lats.v    <= 66
  bin.lon <- floor(lons.r360) + 0.5
  bin.lat <- floor(lats.v)    + 0.5
  bin.key <- ifelse(in.domain, paste(bin.lon, bin.lat, sep = "_"), NA)

  uniq.keys <- sort(unique(stats::na.omit(bin.key)))
  bin.idx   <- match(bin.key, uniq.keys)
  nbin      <- length(uniq.keys)

  nx <- dim(lons.r)[1]; ny <- dim(lons.r)[2]; nt <- length(time.r)
  bin.coords <- do.call(rbind, strsplit(uniq.keys, "_"))
  bin.lon.v  <- as.numeric(bin.coords[, 1])
  bin.lat.v  <- as.numeric(bin.coords[, 2])

  chunk.size <- 120
  starts     <- seq(1, nt, by = chunk.size)
  mat.bin    <- matrix(NA_real_, nrow = nbin, ncol = nt)

  for (s in starts) {
    n.read <- min(chunk.size, nt - s + 1)
    cat(sprintf("  chunk %d-%d / %d\n", s, s + n.read - 1, nt))
    chunk <- ncvar_get(nc.r, "somxl030",
                       start = c(1, 1, s), count = c(-1, -1, n.read))
    chunk[chunk > 1e10] <- NA
    chunk.mat <- matrix(chunk, nrow = nx * ny, ncol = n.read)
    rm(chunk)
    valid <- !is.na(bin.idx)
    for (tt in seq_len(n.read)) {
      v  <- chunk.mat[, tt]
      ok <- valid & !is.na(v)
      if (!any(ok)) next
      sums   <- tapply(v[ok], bin.idx[ok], sum)
      counts <- tapply(v[ok], bin.idx[ok], length)
      keys.t <- as.integer(names(sums))
      mat.bin[keys.t, s + tt - 1] <- sums / counts
    }
    rm(chunk.mat); gc()
  }
  nc_close(nc.r)

  message("Computing climatology, anomalies, detrending, winter means per 1-deg cell...")
  clim.mat <- matrix(NA_real_, nrow = nbin, ncol = 12)
  for (m in 1:12) {
    ti <- which(months.r == m)
    clim.mat[, m] <- rowMeans(mat.bin[, ti, drop = FALSE], na.rm = TRUE)
  }
  anom.mat <- mat.bin - clim.mat[, months.r]
  rm(mat.bin); gc()

  tvec <- seq_len(nt)
  for (k in seq_len(nbin)) {
    y <- anom.mat[k, ]
    ok <- !is.na(y)
    if (sum(ok) < 24) { anom.mat[k, ] <- NA_real_; next }
    f  <- lm(y ~ tvec, subset = ok)
    anom.mat[k, ok] <- residuals(f)
  }

  winyrs.r <- ifelse(months.r %in% c(11, 12), years.r + 1L, years.r)
  in.win   <- months.r %in% c(11, 12, 1, 2, 3)
  uy       <- sort(unique(winyrs.r[in.win]))
  win.mld.1 <- matrix(NA_real_, nrow = nbin, ncol = length(uy))
  for (k in seq_along(uy)) {
    ti <- which(winyrs.r == uy[k] & in.win)
    if (length(ti) == 0) next
    win.mld.1[, k] <- rowMeans(anom.mat[, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(anom.mat); gc()

  safe_mean <- function(v) {
    if (any(is.na(v))) return(NA_real_)
    mean(v)
  }
  roll_mean_cell <- function(x, width = 15) {
    rollapply(x, width = width, fill = NA, align = "right", FUN = safe_mean)
  }

  al.df   <- al %>% select(year, AL.sd) %>% filter(!is.na(AL.sd))
  yr.use  <- intersect(uy, al.df$year)
  x.vec   <- al.df$AL.sd[match(yr.use, al.df$year)]
  col.idx <- match(yr.use, uy)

  good.cells <- which(rowSums(!is.na(win.mld.1)) >= max(25, 0.5 * length(uy)))
  message(sprintf("Computing 15-yr rolling MEAN + GLS at %d one-degree bins...",
                  length(good.cells)))

  rmean.mat <- matrix(NA_real_, nrow = length(good.cells), ncol = length(yr.use))
  for (k in seq_along(good.cells)) {
    a <- roll_mean_cell(win.mld.1[good.cells[k], ])
    rmean.mat[k, ] <- a[col.idx]
  }
  rm(win.mld.1); gc()

  beta.vec <- rep(NA_real_, length(good.cells))
  pval.vec <- rep(NA_real_, length(good.cells))
  df.reg   <- data.frame(y = NA_real_, x = x.vec)
  for (k in seq_along(good.cells)) {
    y <- rmean.mat[k, ]
    if (sum(!is.na(y)) < 10) next
    df.reg$y <- y
    fit <- tryCatch(
      gls(y ~ x, data = df.reg, correlation = corAR1(form = ~ 1),
          method = "ML", na.action = na.omit),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    s <- summary(fit)$tTable
    if (nrow(s) < 2) next
    beta.vec[k] <- s[2, 1]
    pval.vec[k] <- s[2, 4]
  }

  beta.full <- rep(NA_real_, nbin); pval.full <- rep(NA_real_, nbin)
  beta.full[good.cells] <- beta.vec
  pval.full[good.cells] <- pval.vec

  mld.anom.1deg.df <- data.frame(
    lon  = bin.lon.v,
    lat  = bin.lat.v,
    beta = beta.full,
    pval = pval.full
  ) %>% filter(!is.na(beta))

  write.csv(mld.anom.1deg.df, mld.anom.1deg.cache, row.names = FALSE)
  message("Saved: ", mld.anom.1deg.cache)
}

mld.anom.1deg.df <- mld.anom.1deg.df %>%
  filter(lon >= 160, lon <= 250, lat >= 20, lat <= 66)
mld.anom.1deg.df$qval <- p.adjust(mld.anom.1deg.df$pval, method = "BH")

n.sig.raw <- sum(mld.anom.1deg.df$pval <= 0.05, na.rm = TRUE)
n.sig.fdr <- sum(mld.anom.1deg.df$qval <= 0.10, na.rm = TRUE)
message(sprintf("1-deg MLD anom: raw p<=0.05: %d | BH-FDR q<=0.10: %d (of %d)",
                n.sig.raw, n.sig.fdr, nrow(mld.anom.1deg.df)))

mld.anom.1deg.grid <- mld.anom.1deg.df %>%
  complete(lon = seq(160.5, 249.5, by = 1),
           lat = seq(20.5,  65.5,  by = 1))

col.lim.anom.1 <- max(abs(mld.anom.1deg.df$beta), na.rm = TRUE)

p.mld.anom.1deg <- ggplot() +
  geom_raster(data = filter(mld.anom.1deg.grid, !is.na(beta)),
              aes(x = lon, y = lat, fill = beta)) +
  geom_contour(data = filter(mld.anom.1deg.grid, !is.na(qval)),
               aes(x = lon, y = lat, z = qval),
               breaks = 0.10, color = "black", linewidth = 0.4) +
  geom_polygon(data = mapWorld.clean,
               aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray30", linewidth = 0.25) +
  geom_path(data = al.box.df, aes(x = x, y = y),
            color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue4", mid = "white",
                       high = "darkorange4", midpoint = 0,
                       limits = c(-col.lim.anom.1, col.lim.anom.1),
                       name = "\u03b2 (m / z)") +
  coord_cartesian(xlim = c(160, 250), ylim = c(20, 66)) +
  scale_x_continuous(breaks = seq(160, 250, 20), labels = lon_label) +
  scale_y_continuous(breaks = seq(20, 60, 10),
                     labels = function(y) paste0(y, "\u00b0N")) +
  labs(title    = "1\u00b0-aggregated winter MLD anomaly regressed on AL SLP SD (15-yr rolling)",
       subtitle = "GLS AR(1); contour = Benjamini-Hochberg FDR q \u2264 0.10",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/MLD_anom_AL_SLP_SD_regression_1deg.png",
       plot = p.mld.anom.1deg, width = 9, height = 6, dpi = 300)
message("Saved: Figures/MLD_anom_AL_SLP_SD_regression_1deg.png")

# ============================================================
# SECTION 14: 2-DEG AGGREGATED MLD anomaly regressed on AL SLP SD
# ============================================================
# Identical pipeline to Section 13 but with MLD averaged into
# 2 deg x 2 deg spatial blocks.

mld.anom.2deg.cache <- "./Output/mld_anom_al_slp_sd_regression_2deg.csv"

if (file.exists(mld.anom.2deg.cache)) {
  message("Loading cached 2-deg MLD anom ~ AL SLP SD regression from: ",
          mld.anom.2deg.cache)
  mld.anom.2deg.df <- read.csv(mld.anom.2deg.cache)
} else {
  message("Loading raw MLD NC (streaming, aggregating to 2-deg) ...")
  nc.r   <- nc_open(raw.mld.file)
  lons.r <- ncvar_get(nc.r, "nav_lon")
  lats.r <- ncvar_get(nc.r, "nav_lat")
  time.r <- ncvar_get(nc.r, "time_counter")
  tun.r  <- ncatt_get(nc.r, "time_counter", "units")$value
  origin.r <- sub(".*since ", "", tun.r)
  dates.r  <- as.Date(as.POSIXct(time.r, origin = origin.r, tz = "UTC"))
  years.r  <- as.integer(format(dates.r, "%Y"))
  months.r <- as.integer(format(dates.r, "%m"))

  lons.r360 <- ifelse(as.vector(lons.r) < 0, as.vector(lons.r) + 360,
                      as.vector(lons.r))
  lats.v    <- as.vector(lats.r)

  in.domain <- lons.r360 >= 110 & lons.r360 <= 250 &
               lats.v    >= 20  & lats.v    <= 66
  bin.lon <- floor(lons.r360 / 2) * 2 + 1
  bin.lat <- floor(lats.v    / 2) * 2 + 1
  bin.key <- ifelse(in.domain, paste(bin.lon, bin.lat, sep = "_"), NA)

  uniq.keys <- sort(unique(stats::na.omit(bin.key)))
  bin.idx   <- match(bin.key, uniq.keys)
  nbin      <- length(uniq.keys)

  nx <- dim(lons.r)[1]; ny <- dim(lons.r)[2]; nt <- length(time.r)
  bin.coords <- do.call(rbind, strsplit(uniq.keys, "_"))
  bin.lon.v  <- as.numeric(bin.coords[, 1])
  bin.lat.v  <- as.numeric(bin.coords[, 2])

  chunk.size <- 120
  starts     <- seq(1, nt, by = chunk.size)
  mat.bin    <- matrix(NA_real_, nrow = nbin, ncol = nt)

  for (s in starts) {
    n.read <- min(chunk.size, nt - s + 1)
    cat(sprintf("  chunk %d-%d / %d\n", s, s + n.read - 1, nt))
    chunk <- ncvar_get(nc.r, "somxl030",
                       start = c(1, 1, s), count = c(-1, -1, n.read))
    chunk[chunk > 1e10] <- NA
    chunk.mat <- matrix(chunk, nrow = nx * ny, ncol = n.read)
    rm(chunk)
    valid <- !is.na(bin.idx)
    for (tt in seq_len(n.read)) {
      v  <- chunk.mat[, tt]
      ok <- valid & !is.na(v)
      if (!any(ok)) next
      sums   <- tapply(v[ok], bin.idx[ok], sum)
      counts <- tapply(v[ok], bin.idx[ok], length)
      keys.t <- as.integer(names(sums))
      mat.bin[keys.t, s + tt - 1] <- sums / counts
    }
    rm(chunk.mat); gc()
  }
  nc_close(nc.r)

  message("Computing climatology, anomalies, detrending, winter means per 2-deg cell...")
  clim.mat <- matrix(NA_real_, nrow = nbin, ncol = 12)
  for (m in 1:12) {
    ti <- which(months.r == m)
    clim.mat[, m] <- rowMeans(mat.bin[, ti, drop = FALSE], na.rm = TRUE)
  }
  anom.mat <- mat.bin - clim.mat[, months.r]
  rm(mat.bin); gc()

  tvec <- seq_len(nt)
  for (k in seq_len(nbin)) {
    y <- anom.mat[k, ]
    ok <- !is.na(y)
    if (sum(ok) < 24) { anom.mat[k, ] <- NA_real_; next }
    f  <- lm(y ~ tvec, subset = ok)
    anom.mat[k, ok] <- residuals(f)
  }

  winyrs.r <- ifelse(months.r %in% c(11, 12), years.r + 1L, years.r)
  in.win   <- months.r %in% c(11, 12, 1, 2, 3)
  uy       <- sort(unique(winyrs.r[in.win]))
  win.mld.2 <- matrix(NA_real_, nrow = nbin, ncol = length(uy))
  for (k in seq_along(uy)) {
    ti <- which(winyrs.r == uy[k] & in.win)
    if (length(ti) == 0) next
    win.mld.2[, k] <- rowMeans(anom.mat[, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(anom.mat); gc()

  safe_mean <- function(v) {
    if (any(is.na(v))) return(NA_real_)
    mean(v)
  }
  roll_mean_cell <- function(x, width = 15) {
    rollapply(x, width = width, fill = NA, align = "right", FUN = safe_mean)
  }

  al.df   <- al %>% select(year, AL.sd) %>% filter(!is.na(AL.sd))
  yr.use  <- intersect(uy, al.df$year)
  x.vec   <- al.df$AL.sd[match(yr.use, al.df$year)]
  col.idx <- match(yr.use, uy)

  good.cells <- which(rowSums(!is.na(win.mld.2)) >= max(25, 0.5 * length(uy)))
  message(sprintf("Computing 15-yr rolling MEAN + GLS at %d two-degree bins...",
                  length(good.cells)))

  rmean.mat <- matrix(NA_real_, nrow = length(good.cells), ncol = length(yr.use))
  for (k in seq_along(good.cells)) {
    a <- roll_mean_cell(win.mld.2[good.cells[k], ])
    rmean.mat[k, ] <- a[col.idx]
  }
  rm(win.mld.2); gc()

  beta.vec <- rep(NA_real_, length(good.cells))
  pval.vec <- rep(NA_real_, length(good.cells))
  df.reg   <- data.frame(y = NA_real_, x = x.vec)
  for (k in seq_along(good.cells)) {
    y <- rmean.mat[k, ]
    if (sum(!is.na(y)) < 10) next
    df.reg$y <- y
    fit <- tryCatch(
      gls(y ~ x, data = df.reg, correlation = corAR1(form = ~ 1),
          method = "ML", na.action = na.omit),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    s <- summary(fit)$tTable
    if (nrow(s) < 2) next
    beta.vec[k] <- s[2, 1]
    pval.vec[k] <- s[2, 4]
  }

  beta.full <- rep(NA_real_, nbin); pval.full <- rep(NA_real_, nbin)
  beta.full[good.cells] <- beta.vec
  pval.full[good.cells] <- pval.vec

  mld.anom.2deg.df <- data.frame(
    lon  = bin.lon.v,
    lat  = bin.lat.v,
    beta = beta.full,
    pval = pval.full
  ) %>% filter(!is.na(beta))

  write.csv(mld.anom.2deg.df, mld.anom.2deg.cache, row.names = FALSE)
  message("Saved: ", mld.anom.2deg.cache)
}

mld.anom.2deg.df <- mld.anom.2deg.df %>%
  filter(lon >= 160, lon <= 250, lat >= 20, lat <= 66)
mld.anom.2deg.df$qval <- p.adjust(mld.anom.2deg.df$pval, method = "BH")

n.sig.raw <- sum(mld.anom.2deg.df$pval <= 0.05, na.rm = TRUE)
n.sig.fdr <- sum(mld.anom.2deg.df$qval <= 0.10, na.rm = TRUE)
message(sprintf("2-deg MLD anom: raw p<=0.05: %d | BH-FDR q<=0.10: %d (of %d)",
                n.sig.raw, n.sig.fdr, nrow(mld.anom.2deg.df)))

mld.anom.2deg.grid <- mld.anom.2deg.df %>%
  complete(lon = seq(161, 249, by = 2),
           lat = seq(21,  65,  by = 2))

col.lim.anom.2 <- max(abs(mld.anom.2deg.df$beta), na.rm = TRUE)

p.mld.anom.2deg <- ggplot() +
  geom_raster(data = filter(mld.anom.2deg.grid, !is.na(beta)),
              aes(x = lon, y = lat, fill = beta)) +
  geom_contour(data = filter(mld.anom.2deg.grid, !is.na(qval)),
               aes(x = lon, y = lat, z = qval),
               breaks = 0.10, color = "black", linewidth = 0.4) +
  geom_polygon(data = mapWorld.clean,
               aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray30", linewidth = 0.25) +
  geom_path(data = al.box.df, aes(x = x, y = y),
            color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue4", mid = "white",
                       high = "darkorange4", midpoint = 0,
                       limits = c(-col.lim.anom.2, col.lim.anom.2),
                       name = "\u03b2 (m / z)") +
  coord_cartesian(xlim = c(160, 250), ylim = c(20, 66)) +
  scale_x_continuous(breaks = seq(160, 250, 20), labels = lon_label) +
  scale_y_continuous(breaks = seq(20, 60, 10),
                     labels = function(y) paste0(y, "\u00b0N")) +
  labs(title    = "2\u00b0-aggregated winter MLD anomaly regressed on AL SLP SD (15-yr rolling)",
       subtitle = "GLS AR(1); contour = Benjamini-Hochberg FDR q \u2264 0.10",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/MLD_anom_AL_SLP_SD_regression_2deg.png",
       plot = p.mld.anom.2deg, width = 9, height = 6, dpi = 300)
message("Saved: Figures/MLD_anom_AL_SLP_SD_regression_2deg.png")

# ============================================================
# SECTION 15: Cellwise CORRELATION of MLD AR(1) with AL SLP SD
# ============================================================
# Same inputs as Section 9 (native-resolution winter detrended MLD
# anomaly; AL SLP SD 15-yr rolling z-scored) but stores Pearson r and
# per-cell p-values from the Modified Chelton (Pyper & Peterman 1998)
# effective-sample-size adjustment, then BH-FDR control.

mld.ar1.cor.cache <- "./Output/mld_ar1_al_slp_sd_correlation_chelton.csv"

if (file.exists(mld.ar1.cor.cache)) {
  message("Loading cached MLD AR(1) ~ AL SLP SD correlation from: ",
          mld.ar1.cor.cache)
  mld.ar1.cor.df <- read.csv(mld.ar1.cor.cache)
} else {
  message("Loading detrended winter-month MLD anomalies for correlation map...")
  nc.m    <- nc_open(mld.winmon.file)
  mld.w   <- ncvar_get(nc.m, "somxl030")
  lons.m  <- ncvar_get(nc.m, "nav_lon")
  lats.m  <- ncvar_get(nc.m, "nav_lat")
  time.m  <- ncvar_get(nc.m, "time_counter")
  tun.m   <- ncatt_get(nc.m, "time_counter", "units")$value
  nc_close(nc.m)

  mld.w[mld.w > 1e10] <- NA
  origin.m  <- sub(".*since ", "", tun.m)
  dates.m   <- as.Date(as.POSIXct(time.m, origin = origin.m, tz = "UTC"))
  years.m   <- as.integer(format(dates.m, "%Y"))
  months.m  <- as.integer(format(dates.m, "%m"))
  winyrs.m  <- ifelse(months.m %in% c(11, 12), years.m + 1L, years.m)
  lons.m360 <- ifelse(lons.m < 0, lons.m + 360, lons.m)

  nx <- dim(mld.w)[1]; ny <- dim(mld.w)[2]; nt <- dim(mld.w)[3]
  mat <- matrix(mld.w, nrow = nx * ny, ncol = nt)
  rm(mld.w); gc()

  uy <- sort(unique(winyrs.m))
  win.mld <- matrix(NA_real_, nrow = nx * ny, ncol = length(uy))
  for (k in seq_along(uy)) {
    ti <- which(winyrs.m == uy[k])
    win.mld[, k] <- rowMeans(mat[, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(mat); gc()

  safe_ar1 <- function(v) {
    if (any(is.na(v))) return(NA_real_)
    suppressWarnings(acf(v, lag.max = 1, plot = FALSE)$acf[2])
  }
  roll_ar1_cell <- function(x, width = 15) {
    rollapply(x, width = width, fill = NA, align = "right", FUN = safe_ar1)
  }

  al.df   <- al %>% select(year, AL.sd) %>% filter(!is.na(AL.sd))
  yr.use  <- intersect(uy, al.df$year)
  x.vec   <- al.df$AL.sd[match(yr.use, al.df$year)]
  col.idx <- match(yr.use, uy)

  # Two-sided Modified Chelton test: uses existing chelton_test() helper
  # (Section 5) for n.eff and t-stat, then converts upper-tail one-sided
  # pval to two-sided so negative correlations are scored symmetrically.
  chelton_pval_two_sided <- function(x, y) {
    ct <- tryCatch(chelton_test(x, y), error = function(e) NULL)
    if (is.null(ct) || !is.finite(ct$t.obs) || !is.finite(ct$df)) {
      return(c(r = NA_real_, p = NA_real_))
    }
    p2 <- 2 * pt(-abs(ct$t.obs), df = ct$df)
    c(r = ct$r.obs, p = p2)
  }

  good.cells <- which(rowSums(!is.na(win.mld)) >= max(25, 0.5 * length(uy)))
  message(sprintf("Computing cellwise r (Modified Chelton) at %d cells...",
                  length(good.cells)))

  r.vec   <- rep(NA_real_, length(good.cells))
  p.vec   <- rep(NA_real_, length(good.cells))
  for (k in seq_along(good.cells)) {
    a <- roll_ar1_cell(win.mld[good.cells[k], ])
    y <- a[col.idx]
    ok <- !is.na(y) & !is.na(x.vec)
    if (sum(ok) < 10) next
    rp <- chelton_pval_two_sided(x.vec[ok], y[ok])
    r.vec[k] <- rp["r"]
    p.vec[k] <- rp["p"]
  }
  rm(win.mld); gc()

  r.full <- rep(NA_real_, nx * ny); p.full <- rep(NA_real_, nx * ny)
  r.full[good.cells] <- r.vec
  p.full[good.cells] <- p.vec

  mld.ar1.cor.df <- data.frame(
    lon  = as.vector(lons.m360),
    lat  = as.vector(lats.m),
    r    = r.full,
    pval = p.full
  ) %>% filter(!is.na(r))

  write.csv(mld.ar1.cor.df, mld.ar1.cor.cache, row.names = FALSE)
  message("Saved: ", mld.ar1.cor.cache)
}

mld.ar1.cor.df <- mld.ar1.cor.df %>%
  filter(lon >= 110, lon <= 250, lat >= 20, lat <= 66)
mld.ar1.cor.df$qval <- p.adjust(mld.ar1.cor.df$pval, method = "BH")

n.sig.raw <- sum(mld.ar1.cor.df$pval <= 0.05, na.rm = TRUE)
n.sig.fdr <- sum(mld.ar1.cor.df$qval <= 0.05, na.rm = TRUE)
message(sprintf("Correlation (Modified Chelton): raw p<=0.05: %d | BH-FDR q<=0.05: %d (of %d)",
                n.sig.raw, n.sig.fdr, nrow(mld.ar1.cor.df)))

mld.ar1.cor.grid <- mld.ar1.cor.df %>%
  mutate(lon = round(lon / 0.5) * 0.5,
         lat = round(lat / 0.5) * 0.5) %>%
  group_by(lon, lat) %>%
  summarise(r = mean(r, na.rm = TRUE),
            qval = mean(qval, na.rm = TRUE),
            .groups = "drop") %>%
  complete(lon = seq(110, 250, by = 0.5),
           lat = seq(20,  66,  by = 0.5))

p.mld.ar1.cor <- ggplot() +
  geom_raster(data = filter(mld.ar1.cor.grid, !is.na(r)),
              aes(x = lon, y = lat, fill = r)) +
  geom_contour(data = filter(mld.ar1.cor.grid, !is.na(qval)),
               aes(x = lon, y = lat, z = qval),
               breaks = 0.05, color = "black", linewidth = 0.4) +
  geom_polygon(data = mapWorld.clean,
               aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray30", linewidth = 0.25) +
  geom_path(data = al.box.df, aes(x = x, y = y),
            color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue4", mid = "white",
                       high = "darkorange4", midpoint = 0,
                       limits = c(-1, 1), name = "r") +
  coord_map(projection = "rectangular", parameters = 55,
            xlim = c(110, 250), ylim = c(20, 66)) +
  scale_x_continuous(breaks = seq(120, 250, 30), labels = lon_label) +
  scale_y_continuous(breaks = seq(20, 60, 10),
                     labels = function(y) paste0(y, "\u00b0N")) +
  labs(title    = "Cellwise correlation: winter MLD AR(1) vs AL SLP SD (15-yr rolling)",
       subtitle = "Pearson r; Modified Chelton p-values with BH-FDR contour at q \u2264 0.05",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/MLD_AR1_AL_SLP_SD_correlation.png",
       plot = p.mld.ar1.cor, width = 9, height = 6, dpi = 300)
message("Saved: Figures/MLD_AR1_AL_SLP_SD_correlation.png")

# ============================================================
# SECTION 16: ERA-DIFFERENCE MAPS — SST AR(1), MLD AR(1), MLD mean
# ============================================================
# Three-panel cellwise map of era 2 minus era 1 for:
#   1. SST AR(1)   (winter Nov-Mar detrended anomalies)
#   2. MLD AR(1)   (winter Nov-Mar detrended anomalies)
#   3. MLD mean    (era mean of winter detrended MLD anomalies, m)
# Era 1: 1989-2006 (18 winters); Era 2: 2007-2024 (18 winters).
# Ported from ORAS5_MLD_summary.R.

era1.yrs <- 1989:2006
era2.yrs <- 2007:2024
stopifnot(length(era1.yrs) == 18, length(era2.yrs) == 18)

sst.src         <- "./Data/era5_sst_NP_1950_2026.nc"
sst.winmon.file <- "./Data/era5_sst_NP_detr_anom_winmon.nc"
cdo             <- "cdo"

if (!file.exists(sst.winmon.file)) {
  message("CDO: building ERA5 SST detrended winter-month anomalies...")
  system(paste(
    cdo, "selmon,11,12,1,2,3 -detrend -ymonsub", sst.src,
    paste0("-ymonmean ", sst.src),
    sst.winmon.file
  ))
}

compute_era_ar1 <- function(win.mat, all.yrs, era.yrs, min.n = 10) {
  idx <- which(all.yrs %in% era.yrs)
  apply(win.mat[, idx, drop = FALSE], 1, function(v) {
    v <- v[!is.na(v)]
    if (length(v) < min.n) return(NA_real_)
    tryCatch(acf(v, lag.max = 1, plot = FALSE)$acf[2],
             error = function(e) NA_real_)
  })
}
compute_era_mean <- function(win.mat, all.yrs, era.yrs, min.n = 10) {
  idx <- which(all.yrs %in% era.yrs)
  apply(win.mat[, idx, drop = FALSE], 1, function(v) {
    v <- v[!is.na(v)]
    if (length(v) < min.n) return(NA_real_)
    mean(v)
  })
}

era.diff.cache <- "./Output/era_diff_sst_ar1_mld_ar1_mld_mean.csv"

if (file.exists(era.diff.cache)) {
  message("Loading cached era-difference maps: ", era.diff.cache)
  era.diff <- read.csv(era.diff.cache)
} else {
  # --- SST ---
  message("Loading ERA5 SST detrended winter anomalies...")
  nc.s      <- nc_open(sst.winmon.file)
  sst.w     <- ncvar_get(nc.s, "sst")
  lons.s    <- ncvar_get(nc.s, "longitude")
  lats.s    <- ncvar_get(nc.s, "latitude")
  all.names <- c(names(nc.s$var), names(nc.s$dim))
  time.var  <- intersect(c("valid_time", "time"), all.names)[1]
  time.s    <- ncvar_get(nc.s, time.var)
  tun.s     <- ncatt_get(nc.s, time.var, "units")$value
  nc_close(nc.s)

  origin.s   <- sub(".*since ", "", tun.s)
  dates.s    <- as.Date(as.POSIXct(time.s, origin = origin.s, tz = "UTC"))
  years.s    <- as.integer(format(dates.s, "%Y"))
  months.s   <- as.integer(format(dates.s, "%m"))
  winyears.s <- ifelse(months.s %in% c(11, 12), years.s + 1L, years.s)
  lons.s360  <- ifelse(lons.s < 0, lons.s + 360, lons.s)

  lon.idx.s <- which(lons.s360 >= 110 & lons.s360 <= 250)
  lat.idx.s <- which(lats.s   >= 20  & lats.s   <= 66)
  sst.w     <- sst.w[lon.idx.s, lat.idx.s, ]
  lons.s360 <- lons.s360[lon.idx.s]
  lats.s    <- lats.s[lat.idx.s]

  stopifnot(all(era1.yrs %in% winyears.s), all(era2.yrs %in% winyears.s))

  nlon.s <- length(lons.s360); nlat.s <- length(lats.s)
  sst.mat <- matrix(sst.w, nrow = nlon.s * nlat.s, ncol = dim(sst.w)[3])
  rm(sst.w); gc()

  all.wyrs.s <- sort(unique(winyears.s))
  sst.win    <- matrix(NA_real_, nrow = nlon.s * nlat.s,
                       ncol = length(all.wyrs.s))
  for (k in seq_along(all.wyrs.s)) {
    ti <- which(winyears.s == all.wyrs.s[k])
    sst.win[, k] <- rowMeans(sst.mat[, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(sst.mat); gc()

  sst.ar1.diff <- compute_era_ar1(sst.win, all.wyrs.s, era2.yrs) -
                  compute_era_ar1(sst.win, all.wyrs.s, era1.yrs)
  rm(sst.win); gc()

  sst.df <- expand.grid(lon = lons.s360, lat = lats.s) %>%
    mutate(diff = sst.ar1.diff, variable = "SST_AR1") %>%
    filter(!is.na(diff))

  # --- MLD (AR(1) and mean) ---
  message("Loading ORAS5 MLD detrended winter anomalies...")
  nc.m   <- nc_open(mld.winmon.file)
  mld.w  <- ncvar_get(nc.m, "somxl030")
  lons.m <- ncvar_get(nc.m, "nav_lon")
  lats.m <- ncvar_get(nc.m, "nav_lat")
  time.m <- ncvar_get(nc.m, "time_counter")
  tun.m  <- ncatt_get(nc.m, "time_counter", "units")$value
  nc_close(nc.m)

  mld.w[mld.w > 1e10] <- NA
  origin.m   <- sub(".*since ", "", tun.m)
  dates.m    <- as.Date(as.POSIXct(time.m, origin = origin.m, tz = "UTC"))
  years.m    <- as.integer(format(dates.m, "%Y"))
  months.m   <- as.integer(format(dates.m, "%m"))
  winyears.m <- ifelse(months.m %in% c(11, 12), years.m + 1L, years.m)
  lons.m360  <- ifelse(lons.m < 0, lons.m + 360, lons.m)

  nx.m <- dim(mld.w)[1]; ny.m <- dim(mld.w)[2]; nt.m <- dim(mld.w)[3]
  np.mask.m <- as.vector(lons.m360 >= 110 & lons.m360 <= 250 &
                         lats.m    >= 20  & lats.m    <= 66)
  mld.mat   <- matrix(mld.w, nrow = nx.m * ny.m, ncol = nt.m)
  good.m    <- which(np.mask.m & rowSums(!is.na(mld.mat)) > nt.m * 0.5)
  rm(mld.w); gc()

  all.wyrs.m <- sort(unique(winyears.m))
  mld.win    <- matrix(NA_real_, nrow = length(good.m),
                       ncol = length(all.wyrs.m))
  for (k in seq_along(all.wyrs.m)) {
    ti <- which(winyears.m == all.wyrs.m[k])
    mld.win[, k] <- rowMeans(mld.mat[good.m, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(mld.mat); gc()

  mld.ar1.diff  <- compute_era_ar1(mld.win,  all.wyrs.m, era2.yrs) -
                   compute_era_ar1(mld.win,  all.wyrs.m, era1.yrs)
  mld.mean.diff <- compute_era_mean(mld.win, all.wyrs.m, era2.yrs) -
                   compute_era_mean(mld.win, all.wyrs.m, era1.yrs)
  rm(mld.win); gc()

  ar1.full  <- rep(NA_real_, nx.m * ny.m)
  mean.full <- rep(NA_real_, nx.m * ny.m)
  ar1.full[good.m]  <- mld.ar1.diff
  mean.full[good.m] <- mld.mean.diff

  mld.ar1.df <- data.frame(
    lon      = as.vector(lons.m360),
    lat      = as.vector(lats.m),
    diff     = ar1.full,
    variable = "MLD_AR1"
  ) %>% filter(!is.na(diff), lon >= 110, lon <= 250, lat >= 20, lat <= 66)

  mld.mean.df <- data.frame(
    lon      = as.vector(lons.m360),
    lat      = as.vector(lats.m),
    diff     = mean.full,
    variable = "MLD_MEAN"
  ) %>% filter(!is.na(diff), lon >= 110, lon <= 250, lat >= 20, lat <= 66)

  era.diff <- bind_rows(sst.df, mld.ar1.df, mld.mean.df)
  write.csv(era.diff, era.diff.cache, row.names = FALSE)
  message("Saved: ", era.diff.cache)
}

# Plot domain: N. Pacific 160-250E, 20-66N (western edge trimmed)
era.plot.sst <- era.diff %>% filter(variable == "SST_AR1", lon >= 160)
era.plot.mar <- era.diff %>% filter(variable == "MLD_AR1", lon >= 160)
era.plot.mme <- era.diff %>% filter(variable == "MLD_MEAN", lon >= 160)

goa.line.df <- data.frame(lon = c(goa.x, goa.x[1]),
                          lat = c(goa.y, goa.y[1]),
                          grp = "GOA")
ebs.line.df <- data.frame(lon = c(ebs.x, ebs.x[1]),
                          lat = c(ebs.y, ebs.y[1]),
                          grp = "EBS")
message("GOA poly range: lon ", paste(range(goa.line.df$lon), collapse = "-"),
        "  lat ", paste(range(goa.line.df$lat), collapse = "-"))
message("EBS poly range: lon ", paste(range(ebs.line.df$lon), collapse = "-"),
        "  lat ", paste(range(ebs.line.df$lat), collapse = "-"))

col.sst <- max(abs(era.plot.sst$diff), na.rm = TRUE)
col.mar <- max(abs(era.plot.mar$diff), na.rm = TRUE)
col.mme <- max(abs(era.plot.mme$diff), na.rm = TRUE)

make_era_panel <- function(df, col.lim, title.label, legend.label) {
  ggplot() +
    geom_point(data = df, aes(x = lon, y = lat, color = diff),
               shape = 15, size = 0.3) +
    geom_polygon(data = mapWorld.clean,
                 aes(x = long, y = lat, group = group),
                 fill = "gray85", color = "gray30", linewidth = 0.25) +
    scale_color_gradient2(low = "steelblue4", mid = "white",
                          high = "darkorange4", midpoint = 0,
                          limits = c(-col.lim, col.lim),
                          name = legend.label) +
    geom_path(data = goa.line.df,
              aes(x = lon, y = lat),
              color = "black", linewidth = 1.0,
              inherit.aes = FALSE) +
    geom_path(data = ebs.line.df,
              aes(x = lon, y = lat),
              color = "black", linewidth = 1.0,
              inherit.aes = FALSE) +
    coord_cartesian(xlim = c(160, 250), ylim = c(20, 66)) +
    scale_x_continuous(breaks = seq(160, 250, 20), labels = lon_label) +
    scale_y_continuous(breaks = seq(20, 60, 10),
                       labels = function(y) paste0(y, "\u00b0N")) +
    labs(title = title.label, x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title      = element_text(hjust = 0.5, size = 10),
          legend.position = "right",
          panel.grid.minor = element_blank())
}

p.era.sst <- make_era_panel(era.plot.sst, col.sst,
                            "\u0394 SST AR(1)",  "\u0394AR(1)")
p.era.mar <- make_era_panel(era.plot.mar, col.mar,
                            "\u0394 MLD AR(1)",  "\u0394AR(1)")
p.era.mme <- make_era_panel(era.plot.mme, col.mme,
                            "\u0394 MLD mean",   "\u0394MLD (m)")

p.era3 <- (p.era.sst | p.era.mar | p.era.mme) +
  plot_annotation(
    title    = "Era differences (era 2 \u2212 era 1)",
    subtitle = paste0("Era 1: 1989\u20132006 (n = 18)  |  ",
                      "Era 2: 2007\u20132024 (n = 18)  |  ",
                      "winter (Nov\u2013Mar) annual means"),
    theme = theme(plot.title    = element_text(hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5, size = 9))
  )

ggsave("./Figures/ERA_diff_SST_AR1_MLD_AR1_MLD_mean.png",
       plot = p.era3, width = 16, height = 5, dpi = 150)
message("Saved: Figures/ERA_diff_SST_AR1_MLD_AR1_MLD_mean.png")

# ============================================================
# SECTION 17: Cellwise regression of SST AR(1) on AL SLP SD
# ============================================================
# Same workflow as Section 9 but with SST AR(1) as the response.
# x: al$AL.sd (15-yr rolling SD of winter AL SLP anomaly, z-scored)
# y: per-cell 15-yr right-aligned rolling AR(1) of winter-mean
#    detrended SST anomaly (from sst.winmon.file).
# No significance overlay.

sst.ar1.reg.cache <- "./Output/sst_ar1_al_slp_sd_regression.csv"

if (file.exists(sst.ar1.reg.cache)) {
  message("Loading cached SST AR(1) ~ AL SLP SD regression: ",
          sst.ar1.reg.cache)
  sst.ar1.reg.df <- read.csv(sst.ar1.reg.cache)
} else {
  message("Loading ERA5 detrended winter-month SST anomalies...")
  nc.s      <- nc_open(sst.winmon.file)
  sst.w     <- ncvar_get(nc.s, "sst")
  lons.s    <- ncvar_get(nc.s, "longitude")
  lats.s    <- ncvar_get(nc.s, "latitude")
  all.names <- c(names(nc.s$var), names(nc.s$dim))
  time.var  <- intersect(c("valid_time", "time"), all.names)[1]
  time.s    <- ncvar_get(nc.s, time.var)
  tun.s     <- ncatt_get(nc.s, time.var, "units")$value
  nc_close(nc.s)

  origin.s   <- sub(".*since ", "", tun.s)
  dates.s    <- as.Date(as.POSIXct(time.s, origin = origin.s, tz = "UTC"))
  years.s    <- as.integer(format(dates.s, "%Y"))
  months.s   <- as.integer(format(dates.s, "%m"))
  winyears.s <- ifelse(months.s %in% c(11, 12), years.s + 1L, years.s)
  lons.s360  <- ifelse(lons.s < 0, lons.s + 360, lons.s)

  lon.idx.s <- which(lons.s360 >= 160 & lons.s360 <= 250)
  lat.idx.s <- which(lats.s   >= 20  & lats.s   <= 66)
  sst.w     <- sst.w[lon.idx.s, lat.idx.s, ]
  lons.s360 <- lons.s360[lon.idx.s]
  lats.s    <- lats.s[lat.idx.s]

  nlon.s <- length(lons.s360); nlat.s <- length(lats.s)
  sst.mat <- matrix(sst.w, nrow = nlon.s * nlat.s, ncol = dim(sst.w)[3])
  rm(sst.w); gc()

  all.wyrs.s <- sort(unique(winyears.s))
  sst.win    <- matrix(NA_real_, nrow = nlon.s * nlat.s,
                       ncol = length(all.wyrs.s))
  for (k in seq_along(all.wyrs.s)) {
    ti <- which(winyears.s == all.wyrs.s[k])
    sst.win[, k] <- rowMeans(sst.mat[, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(sst.mat); gc()

  safe_ar1 <- function(v) {
    if (any(is.na(v))) return(NA_real_)
    suppressWarnings(acf(v, lag.max = 1, plot = FALSE)$acf[2])
  }
  roll_ar1_cell <- function(x, width = 15) {
    rollapply(x, width = width, fill = NA, align = "right", FUN = safe_ar1)
  }

  al.df   <- al %>% select(year, AL.sd) %>% filter(!is.na(AL.sd))
  yr.use  <- intersect(all.wyrs.s, al.df$year)
  x.vec   <- al.df$AL.sd[match(yr.use, al.df$year)]
  col.idx <- match(yr.use, all.wyrs.s)
  message(sprintf("SST regression years: %d-%d (%d years)",
                  min(yr.use), max(yr.use), length(yr.use)))

  good.cells <- which(rowSums(!is.na(sst.win)) >= max(25, 0.5 * length(all.wyrs.s)))
  message(sprintf("Computing cellwise 15-yr rolling AR(1) at %d cells...",
                  length(good.cells)))

  ar1.mat <- matrix(NA_real_, nrow = length(good.cells), ncol = length(yr.use))
  for (k in seq_along(good.cells)) {
    a <- roll_ar1_cell(sst.win[good.cells[k], ])
    ar1.mat[k, ] <- a[col.idx]
  }
  rm(sst.win); gc()

  message(sprintf("Fitting GLS AR(1) at %d cells...", length(good.cells)))
  beta.vec <- rep(NA_real_, length(good.cells))
  pval.vec <- rep(NA_real_, length(good.cells))
  df.reg   <- data.frame(y = NA_real_, x = x.vec)
  for (k in seq_along(good.cells)) {
    y <- ar1.mat[k, ]
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

  beta.full <- rep(NA_real_, nlon.s * nlat.s)
  pval.full <- rep(NA_real_, nlon.s * nlat.s)
  beta.full[good.cells] <- beta.vec
  pval.full[good.cells] <- pval.vec

  lonlat <- expand.grid(lon = lons.s360, lat = lats.s)
  sst.ar1.reg.df <- data.frame(lon  = lonlat$lon,
                               lat  = lonlat$lat,
                               beta = beta.full,
                               pval = pval.full) %>%
    filter(!is.na(beta))

  write.csv(sst.ar1.reg.df, sst.ar1.reg.cache, row.names = FALSE)
  message("Saved: ", sst.ar1.reg.cache)
}

sst.ar1.reg.df <- sst.ar1.reg.df %>% filter(lon >= 160, lon <= 250)
sst.ar1.reg.df$qval <- p.adjust(sst.ar1.reg.df$pval, method = "BH")

n.sig.raw.sst <- sum(sst.ar1.reg.df$pval <= 0.05, na.rm = TRUE)
n.sig.fdr.sst <- sum(sst.ar1.reg.df$qval <= 0.05, na.rm = TRUE)
message(sprintf("SST AR(1) ~ AL SLP SD: raw p<=0.05: %d | BH-FDR q<=0.05: %d (of %d)",
                n.sig.raw.sst, n.sig.fdr.sst, nrow(sst.ar1.reg.df)))

sst.ar1.reg.grid <- sst.ar1.reg.df %>%
  mutate(lon = round(lon / 0.5) * 0.5,
         lat = round(lat / 0.5) * 0.5) %>%
  group_by(lon, lat) %>%
  summarise(beta = mean(beta, na.rm = TRUE),
            qval = mean(qval, na.rm = TRUE),
            .groups = "drop") %>%
  complete(lon = seq(160, 250, by = 0.5),
           lat = seq(20,  66,  by = 0.5))

col.lim.sst <- max(abs(sst.ar1.reg.df$beta), na.rm = TRUE)

p.sst.ar1.reg <- ggplot() +
  geom_raster(data = filter(sst.ar1.reg.grid, !is.na(beta)),
              aes(x = lon, y = lat, fill = beta)) +
  geom_contour(data = filter(sst.ar1.reg.grid, !is.na(qval)),
               aes(x = lon, y = lat, z = qval),
               breaks = 0.05, color = "black", linewidth = 0.4) +
  geom_polygon(data = mapWorld.clean,
               aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray30", linewidth = 0.25) +
  geom_path(data = al.box.df, aes(x = x, y = y),
            color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue4", mid = "white",
                       high = "darkorange4", midpoint = 0,
                       limits = c(-col.lim.sst, col.lim.sst),
                       name = "\u03b2 (AR(1) / z)") +
  coord_cartesian(xlim = c(160, 250), ylim = c(20, 66)) +
  scale_x_continuous(breaks = seq(160, 250, 20), labels = lon_label) +
  scale_y_continuous(breaks = seq(20, 60, 10),
                     labels = function(y) paste0(y, "\u00b0N")) +
  labs(title    = "Winter SST AR(1) regressed on AL SLP SD (15-yr rolling windows)",
       subtitle = "GLS AR(1); contour = Benjamini-Hochberg FDR q \u2264 0.05 (Wilks 2016)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/SST_AR1_AL_SLP_SD_regression.png",
       plot = p.sst.ar1.reg, width = 9, height = 6, dpi = 300)
message("Saved: Figures/SST_AR1_AL_SLP_SD_regression.png")

# ============================================================
# SECTION 18: SEASONAL LAG SENSITIVITY (Newman et al. 2016, J. Climate)
# ============================================================
# Same dual-axis plots as Section 5 (15-yr AL SLP SD vs SST AR(1) for
# GOA and EBS), but with a seasonal lag between AL forcing and ocean
# response: AL SLP SD uses Nov-Dec-Jan means; SST AR(1) uses Feb-Mar-Apr
# means. Year = January year for both, so they refer to the same winter.

# --- AL SLP SD with N-D-J season ---
al.ndj <- al.monthly %>%
  filter(month %in% c(11, 12, 1)) %>%
  mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
  group_by(win.year) %>%
  summarise(SLP = mean(SLP, na.rm = TRUE), .groups = "drop") %>%
  rename(year = win.year) %>%
  arrange(year) %>%
  mutate(AL.sd = roll_sd(SLP),
         AL.sd = as.numeric(scale(AL.sd)))

# --- SST AR(1) with F-M-A season (year = calendar year = January year) ---
fma.anom <- monthly.anom %>%
  filter(month %in% c(2, 3, 4)) %>%
  group_by(year) %>%
  summarise(GOA = mean(GOA.anom, na.rm = TRUE),
            EBS = mean(EBS.anom, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(year)

fma.detr <- fma.anom %>%
  mutate(GOA = detrend(GOA),
         EBS = detrend(EBS))

fma.roll <- fma.detr %>%
  mutate(GOA.ar1 = roll_ar1(GOA),
         EBS.ar1 = roll_ar1(EBS))

al.sst.lag <- fma.roll %>%
  select(year, GOA.ar1, EBS.ar1) %>%
  left_join(al.ndj %>% select(year, AL.sd), by = "year") %>%
  filter(!is.na(AL.sd))

al.sd.range.lag   <- range(al.sst.lag$AL.sd,   na.rm = TRUE)
goa.ar1.range.lag <- range(al.sst.lag$GOA.ar1, na.rm = TRUE)
ebs.ar1.range.lag <- range(al.sst.lag$EBS.ar1, na.rm = TRUE)

dat.goa.lag <- al.sst.lag %>% filter(!is.na(GOA.ar1), !is.na(AL.sd))
dat.ebs.lag <- al.sst.lag %>% filter(!is.na(EBS.ar1), !is.na(AL.sd))

chel.goa.lag <- chelton_test(dat.goa.lag$GOA.ar1, dat.goa.lag$AL.sd)
chel.ebs.lag <- chelton_test(dat.ebs.lag$EBS.ar1, dat.ebs.lag$AL.sd)

message(sprintf("GOA NDJ-FMA Chelton: r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
                chel.goa.lag$r.obs, chel.goa.lag$n.eff, chel.goa.lag$t.obs, chel.goa.lag$pval))
message(sprintf("EBS NDJ-FMA Chelton: r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
                chel.ebs.lag$r.obs, chel.ebs.lag$n.eff, chel.ebs.lag$t.obs, chel.ebs.lag$pval))

make_dual_axis_plot_lag <- function(dat, ar1.col, ar1.color, region.label, ar1.range, chel) {
  k <- diff(ar1.range) / diff(al.sd.range.lag)
  b <- ar1.range[1] - k * al.sd.range.lag[1]

  dat <- dat %>% filter(!is.na(.data[[ar1.col]]), !is.na(AL.sd))

  lbl <- sprintf("r = %.2f, N* = %.1f\nt = %.2f, p = %.3f",
                 chel$r.obs, chel$n.eff, chel$t.obs, chel$pval)

  ggplot(dat, aes(x = year)) +
    geom_line(aes(y = .data[[ar1.col]], color = "SST AR(1) [FMA]"), linewidth = 0.7) +
    geom_point(aes(y = .data[[ar1.col]], color = "SST AR(1) [FMA]"), size = 1.5) +
    geom_line(aes(y = k * AL.sd + b, color = "AL SLP SD [NDJ]"),
              linewidth = 0.7, linetype = "dashed") +
    geom_point(aes(y = k * AL.sd + b, color = "AL SLP SD [NDJ]"), size = 1.5, shape = 1) +
    annotate("text", x = Inf, y = -Inf, label = lbl,
             hjust = 1.1, vjust = -0.4, size = 3.5, lineheight = 1.3) +
    scale_y_continuous(
      name     = "SST AR(1)  [FMA]",
      sec.axis = sec_axis(~ (. - b) / k, name = "AL SLP SD (z-scored)  [NDJ]")
    ) +
    scale_color_manual(values = c("SST AR(1) [FMA]" = ar1.color,
                                  "AL SLP SD [NDJ]" = "gray30")) +
    labs(title = region.label, x = "Year", color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position    = "bottom",
          axis.title.y.left  = element_text(color = ar1.color),
          axis.text.y.left   = element_text(color = ar1.color),
          axis.title.y.right = element_text(color = "gray30"),
          axis.text.y.right  = element_text(color = "gray30"))
}

p.goa.dual.lag <- make_dual_axis_plot_lag(al.sst.lag, "GOA.ar1", "darkorange4",
                                           "GOA", goa.ar1.range.lag, chel.goa.lag)
p.ebs.dual.lag <- make_dual_axis_plot_lag(al.sst.lag, "EBS.ar1", "steelblue4",
                                           "EBS", ebs.ar1.range.lag, chel.ebs.lag)

p.dual.lag <- p.goa.dual.lag / p.ebs.dual.lag +
  plot_annotation(
    title = "15-yr Rolling AL SLP SD (NDJ) vs. SST AR(1) (FMA) \u2014 Newman et al. 2016 lag",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 12))
  )

ggsave("./Figures/AL_SD_NDJ_SST_AR1_FMA_15yr_dual_axis.png",
       plot = p.dual.lag, width = 7, height = 7, dpi = 300)
message("Saved: Figures/AL_SD_NDJ_SST_AR1_FMA_15yr_dual_axis.png")

# ------------------------------------------------------------
# Section 18 INDEPENDENT CHECKS
# ------------------------------------------------------------
# Verify: (1) month filtering & year-alignment; (2) recompute one
# NDJ mean from monthly table; (3) recompute one rolling SD value
# by hand; (4) compare NDJ vs NDJFM AL series; (5) overlay plot.

message("\n========== Section 18: Independent checks ==========")

# --- Check 1: month counts per year for NDJ and FMA ---
ndj.counts <- al.monthly %>%
  filter(month %in% c(11, 12, 1)) %>%
  mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
  count(win.year, name = "n.months")

fma.counts <- monthly.anom %>%
  filter(month %in% c(2, 3, 4)) %>%
  count(year, name = "n.months")

message(sprintf("NDJ month counts: years with 3 months = %d; <3 months = %d (edge years: %s)",
                sum(ndj.counts$n.months == 3),
                sum(ndj.counts$n.months <  3),
                paste(ndj.counts$win.year[ndj.counts$n.months < 3], collapse = ", ")))
message(sprintf("FMA month counts: years with 3 months = %d; <3 months = %d (edge years: %s)",
                sum(fma.counts$n.months == 3),
                sum(fma.counts$n.months <  3),
                paste(fma.counts$year[fma.counts$n.months < 3], collapse = ", ")))

# --- Check 2: hand-recompute one NDJ mean (target year = 2000 if available) ---
target.yr <- 2000
ndj.hand <- mean(al.monthly$SLP[
  (al.monthly$year == target.yr - 1L & al.monthly$month %in% c(11, 12)) |
  (al.monthly$year == target.yr      & al.monthly$month == 1)
], na.rm = TRUE)
ndj.code <- al.ndj$SLP[al.ndj$year == target.yr]
message(sprintf("NDJ %d hand=%.6f | code=%.6f | diff=%.2e",
                target.yr, ndj.hand, ndj.code, ndj.hand - ndj.code))

# --- Check 3: hand-recompute one rolling-SD value (right-aligned 15-yr) ---
# For year y, roll_sd uses years (y-14) through y of the SLP series (pre-scale).
y.test <- 2010
slp.window <- al.ndj %>% filter(year >= y.test - 14, year <= y.test) %>% pull(SLP)
sd.hand <- if (length(slp.window) == 15) sd(slp.window) else NA_real_
sd.code.raw <- roll_sd(al.ndj$SLP)[al.ndj$year == y.test]
message(sprintf("Rolling SD %d window length=%d | hand=%.6f | code=%.6f | diff=%.2e",
                y.test, length(slp.window), sd.hand, sd.code.raw, sd.hand - sd.code.raw))

# --- Check 4: compare NDJ (3-mo) annual SLP vs NDJFM (5-mo) annual SLP ---
slp.compare <- al %>%
  select(year, SLP.ndjfm = SLP) %>%
  left_join(al.ndj %>% select(year, SLP.ndj = SLP), by = "year") %>%
  filter(!is.na(SLP.ndj), !is.na(SLP.ndjfm))

cor.slp <- cor(slp.compare$SLP.ndj, slp.compare$SLP.ndjfm)
message(sprintf("Correlation of annual mean SLP, NDJ vs NDJFM: r = %.3f (n=%d) -- expect high (~0.85-0.95)",
                cor.slp, nrow(slp.compare)))

sd.compare <- al %>%
  select(year, sd.ndjfm = AL.sd) %>%
  left_join(al.ndj %>% select(year, sd.ndj = AL.sd), by = "year") %>%
  filter(!is.na(sd.ndj), !is.na(sd.ndjfm))

cor.sd <- cor(sd.compare$sd.ndj, sd.compare$sd.ndjfm)
message(sprintf("Correlation of 15-yr rolling AL SD (z-scored), NDJ vs NDJFM: r = %.3f (n=%d)",
                cor.sd, nrow(sd.compare)))

# --- Check 4b: same correlation on RAW (pre-scale) rolling SDs, in case
#     scale() amplified noise, plus side-by-side preview every 5 yrs ---
raw.sd.compare <- data.frame(
  year         = al$year,
  raw.sd.ndjfm = roll_sd(al$SLP)
) %>%
  left_join(data.frame(year = al.ndj$year,
                       raw.sd.ndj = roll_sd(al.ndj$SLP)),
            by = "year") %>%
  filter(!is.na(raw.sd.ndjfm), !is.na(raw.sd.ndj))

cor.sd.raw <- cor(raw.sd.compare$raw.sd.ndj, raw.sd.compare$raw.sd.ndjfm)
message(sprintf("Correlation of 15-yr rolling AL SD on RAW (pre-scale) SLP: r = %.3f (n=%d)",
                cor.sd.raw, nrow(raw.sd.compare)))

# Lag check: does shifting NDJ by +/-1 yr improve agreement?
y.both <- raw.sd.compare %>% arrange(year)
cor.lag.m1 <- cor(y.both$raw.sd.ndj[-1], y.both$raw.sd.ndjfm[-nrow(y.both)])
cor.lag.p1 <- cor(y.both$raw.sd.ndj[-nrow(y.both)], y.both$raw.sd.ndjfm[-1])
message(sprintf("  NDJ shifted -1 yr: r = %.3f | NDJ shifted +1 yr: r = %.3f (large diff => alignment bug)",
                cor.lag.m1, cor.lag.p1))

message("Rolling SD preview (every 5 yr) -- year | raw.NDJFM | raw.NDJ:")
preview <- raw.sd.compare[seq(1, nrow(raw.sd.compare), by = 5), ]
for (i in seq_len(nrow(preview))) {
  message(sprintf("  %4d  %7.4f   %7.4f",
                  preview$year[i], preview$raw.sd.ndjfm[i], preview$raw.sd.ndj[i]))
}

# --- Check 5: post-scale AL.sd should have mean ~0, sd ~1 (NaNs ignored) ---
sd.vals <- al.ndj$AL.sd[!is.na(al.ndj$AL.sd)]
message(sprintf("NDJ AL.sd post-scale: mean = %.4f, sd = %.4f (target: 0, 1)",
                mean(sd.vals), sd(sd.vals)))

# --- Diagnostic overlay plot: NDJ vs NDJFM AL SD ---
sd.long <- bind_rows(
  al     %>% select(year, AL.sd) %>% mutate(season = "NDJFM (5-mo)"),
  al.ndj %>% select(year, AL.sd) %>% mutate(season = "NDJ (3-mo)")
) %>% filter(!is.na(AL.sd))

p.sd.overlay <- ggplot(sd.long, aes(x = year, y = AL.sd, color = season)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.3) +
  scale_color_manual(values = c("NDJFM (5-mo)" = "gray30",
                                "NDJ (3-mo)"   = "firebrick")) +
  labs(title = "AL SLP SD: NDJ vs NDJFM (15-yr rolling, z-scored)",
       subtitle = sprintf("r = %.3f over overlapping years", cor.sd),
       x = "Year", y = "AL SLP SD (z-scored)", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("./Figures/AL_SD_NDJ_vs_NDJFM_check.png",
       plot = p.sd.overlay, width = 7, height = 4, dpi = 300)
message("Saved: Figures/AL_SD_NDJ_vs_NDJFM_check.png")

# --- Check 6: drop partial-coverage years from NDJ before rolling SD ---
# Edge year 1950 has only Jan; recomputing rolling SD on full-3-month
# years only should remove any 1964-window inflation.
ndj.full.years <- ndj.counts$win.year[ndj.counts$n.months == 3]
al.ndj.full <- al.ndj %>%
  filter(year %in% ndj.full.years) %>%
  arrange(year) %>%
  mutate(AL.sd.full = roll_sd(SLP))

cmp.full <- al %>%
  select(year, sd.ndjfm = AL.sd) %>%
  left_join(al.ndj.full %>%
              mutate(sd.ndj.full = as.numeric(scale(AL.sd.full))) %>%
              select(year, sd.ndj.full),
            by = "year") %>%
  filter(!is.na(sd.ndjfm), !is.na(sd.ndj.full))

cor.sd.full <- cor(cmp.full$sd.ndj.full, cmp.full$sd.ndjfm)
message(sprintf("NDJ rolling SD (full-coverage years only) vs NDJFM: r = %.3f (n=%d)",
                cor.sd.full, nrow(cmp.full)))
ndj.1964 <- al.ndj.full$AL.sd.full[al.ndj.full$year == 1964]
message(sprintf("  1964 NDJ rolling SD (full-coverage): %.4f (was 7.74 with partial 1950)",
                ifelse(length(ndj.1964) == 0, NA_real_, ndj.1964)))

# --- Check 7: FMA (Feb-Mar-Apr) AL SLP SD vs NDJFM ---
# If NDJFM's decadal signal is dominated by late winter, FMA should track it.
al.fma <- al.monthly %>%
  filter(month %in% c(2, 3, 4)) %>%
  group_by(year) %>%
  summarise(SLP = mean(SLP, na.rm = TRUE), .groups = "drop") %>%
  arrange(year) %>%
  mutate(AL.sd = roll_sd(SLP),
         AL.sd = as.numeric(scale(AL.sd)))

cmp.fma <- al %>%
  select(year, sd.ndjfm = AL.sd) %>%
  left_join(al.fma %>% select(year, sd.fma = AL.sd), by = "year") %>%
  filter(!is.na(sd.ndjfm), !is.na(sd.fma))

cor.sd.fma <- cor(cmp.fma$sd.fma, cmp.fma$sd.ndjfm)
message(sprintf("FMA AL SLP SD vs NDJFM: r = %.3f (n=%d) -- if high, NDJFM signal is driven by late winter",
                cor.sd.fma, nrow(cmp.fma)))

cmp.ndj.fma <- al.ndj %>%
  select(year, sd.ndj = AL.sd) %>%
  left_join(al.fma %>% select(year, sd.fma = AL.sd), by = "year") %>%
  filter(!is.na(sd.ndj), !is.na(sd.fma))

cor.ndj.fma <- cor(cmp.ndj.fma$sd.ndj, cmp.ndj.fma$sd.fma)
message(sprintf("NDJ AL SLP SD vs FMA: r = %.3f (n=%d) -- if low, early/late winter volatility decouples",
                cor.ndj.fma, nrow(cmp.ndj.fma)))

# Three-season overlay plot
sd.three <- bind_rows(
  al     %>% select(year, AL.sd) %>% mutate(season = "NDJFM (5-mo)"),
  al.ndj %>% select(year, AL.sd) %>% mutate(season = "NDJ (3-mo)"),
  al.fma %>% select(year, AL.sd) %>% mutate(season = "FMA (3-mo)")
) %>% filter(!is.na(AL.sd))

p.sd.three <- ggplot(sd.three, aes(x = year, y = AL.sd, color = season)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.3) +
  scale_color_manual(values = c("NDJFM (5-mo)" = "gray30",
                                "NDJ (3-mo)"   = "firebrick",
                                "FMA (3-mo)"   = "steelblue4")) +
  labs(title = "AL SLP SD by season (15-yr rolling, z-scored)",
       subtitle = sprintf("NDJFM~NDJ r=%.2f | NDJFM~FMA r=%.2f | NDJ~FMA r=%.2f",
                          cor.sd, cor.sd.fma, cor.ndj.fma),
       x = "Year", y = "AL SLP SD (z-scored)", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("./Figures/AL_SD_three_seasons_check.png",
       plot = p.sd.three, width = 8, height = 4, dpi = 300)
message("Saved: Figures/AL_SD_three_seasons_check.png")
message("===================================================\n")

# ============================================================
# SECTION 19: MLD AR(1) regression vs. PDO (SST EOF1) pattern
# ============================================================
# Compares the cellwise MLD AR(1) ~ AL SLP SD regression map
# (native resolution, from Section 9) against the leading EOF
# of winter (Nov-Mar) detrended monthly ERA5 SST anomalies over
# the same time domain - NO rolling windows for the EOF.
#
# EOF1 of winter monthly NDJFM detrended SST anomalies is the
# canonical PDO-like leading mode of N. Pacific variability.
# Computed on cos(lat)-weighted covariance via truncated SVD
# (irlba), sign flipped to match conventional PDO (negative
# loading in central N. Pacific, positive along NA coast).
#
# Spatial coherence is evaluated on the native MLD-regression grid
# (nearest-neighbor join of EOF cells, no spatial aggregation) via
# THREE robust methods:
#   (1) Pearson r of paired cellwise (beta, EOF1) values, with
#       p-value from spatial-block bootstrap (1000 reps, 5x5-deg
#       blocks) to handle spatial autocorrelation. This is the
#       standard practical approach when raw cells are not
#       independent (Wilks 2011, ch. 14).
#   (2) Tucker congruence coefficient (uncentered cosine; phi).
#       phi >= 0.95 = essentially identical patterns; 0.85-0.94
#       = fair similarity; < 0.85 = patterns differ. Robust to
#       sign and scale; standard in factor-pattern comparison
#       (Lorenzo-Seva & ten Berge 2006, Methodology).
#   (3) Sign-agreement fraction: proportion of common cells
#       where beta and EOF1 share sign. Tested against the null
#       expectation of 0.5 with a binomial test using the
#       spatially effective sample size (Bretherton et al. 1999,
#       J. Climate, eq. for spatial DOF from EOF eigenvalues).

library(irlba)

# --- (a) Load Section 9 MLD AR(1) regression map (native res) ---
mld.ar1.reg.cache <- "./Output/mld_ar1_al_slp_sd_regression.csv"
stopifnot(file.exists(mld.ar1.reg.cache))
mld.beta <- read.csv(mld.ar1.reg.cache) %>%
  filter(lon >= 180, lon <= 250, lat >= 30, lat <= 66, !is.na(beta))

# Time domain that underlies Section 9 (winter years used as the
# regression series; we re-derive from the cached MLD winter file
# so the EOF spans the same calendar window).
nc.m9   <- nc_open(mld.winmon.file)
time.m9 <- ncvar_get(nc.m9, "time_counter")
tun.m9  <- ncatt_get(nc.m9, "time_counter", "units")$value
nc_close(nc.m9)
origin.m9 <- sub(".*since ", "", tun.m9)
dates.m9  <- as.Date(as.POSIXct(time.m9, origin = origin.m9, tz = "UTC"))
mld.years <- as.integer(format(dates.m9, "%Y"))
mld.mons  <- as.integer(format(dates.m9, "%m"))
mld.winyr <- ifelse(mld.mons %in% c(11, 12), mld.years + 1L, mld.years)
sec19.win.range <- range(mld.winyr)
message(sprintf("Section 19: EOF time domain = winter years %d-%d (matches Sec 9 source)",
                sec19.win.range[1], sec19.win.range[2]))

# --- (b) Load winter monthly NDJFM detrended SST anomalies ---
# sst.winmon.file already contains -selmon,11,12,1,2,3 -detrend -ymonsub
# (built in Section 16). No rolling windows.
nc.s   <- nc_open(sst.winmon.file)
sst.w  <- ncvar_get(nc.s, "sst")
lons.s <- ncvar_get(nc.s, "longitude")
lats.s <- ncvar_get(nc.s, "latitude")
all.names.s <- c(names(nc.s$var), names(nc.s$dim))
time.var.s  <- intersect(c("valid_time", "time"), all.names.s)[1]
time.s <- ncvar_get(nc.s, time.var.s)
tun.s  <- ncatt_get(nc.s, time.var.s, "units")$value
nc_close(nc.s)

origin.s <- sub(".*since ", "", tun.s)
dates.s  <- as.Date(as.POSIXct(time.s, origin = origin.s, tz = "UTC"))
years.s  <- as.integer(format(dates.s, "%Y"))
mons.s   <- as.integer(format(dates.s, "%m"))
winyrs.s <- ifelse(mons.s %in% c(11, 12), years.s + 1L, years.s)
lons.s360 <- ifelse(lons.s < 0, lons.s + 360, lons.s)

# Restrict to plotting domain (matches MLD AR(1) map: 180-250E, 30-66N)
lon.idx <- which(lons.s360 >= 180 & lons.s360 <= 250)
lat.idx <- which(lats.s   >= 30  & lats.s   <= 66)
sst.w   <- sst.w[lon.idx, lat.idx, ]
lons.s360 <- lons.s360[lon.idx]
lats.s    <- lats.s[lat.idx]

# Restrict time to winters that overlap MLD window (full time domain
# of Section 9 source, no rolling-window restriction).
t.keep <- which(winyrs.s >= sec19.win.range[1] & winyrs.s <= sec19.win.range[2])
sst.w  <- sst.w[, , t.keep]
nlon <- length(lons.s360); nlat <- length(lats.s); nt <- dim(sst.w)[3]
message(sprintf("Section 19: %d winter months across %d cells for SST EOF",
                nt, nlon * nlat))

# --- (c) Cellwise EOF1 via cos(lat)-weighted truncated SVD ---
sst.mat <- matrix(sst.w, nrow = nlon * nlat, ncol = nt)
rm(sst.w); gc()

# Drop cells with any NA across the time series (land/ice mask)
ok.cells <- which(rowSums(is.na(sst.mat)) == 0)
sst.ok   <- sst.mat[ok.cells, , drop = FALSE]

grid.sst <- expand.grid(lon = lons.s360, lat = lats.s)
w.sst    <- sqrt(cos(grid.sst$lat[ok.cells] * pi / 180))
sv.sst   <- irlba(sweep(sst.ok, 1, w.sst, "*"), nv = 3)

# EOF1 spatial loadings (un-weighted), variance fraction
eof1.sst <- sv.sst$u[, 1] / w.sst
var.frac <- (sv.sst$d^2) / sum(apply(sst.ok, 1, var, na.rm = TRUE) *
                               (nt - 1))   # rough variance share
message(sprintf("Section 19: EOF1 explains ~%.1f%% of total winter SST variance",
                100 * var.frac[1]))

# Sign convention: PDO has negative loading in the central N. Pacific
# (around 180E-200E, 35-50N) and positive along the North American coast.
# Flip if mean loading in central box is positive.
cen.mask <- grid.sst$lon[ok.cells] >= 180 & grid.sst$lon[ok.cells] <= 200 &
            grid.sst$lat[ok.cells] >=  35 & grid.sst$lat[ok.cells] <=  50
if (mean(eof1.sst[cen.mask], na.rm = TRUE) > 0) eof1.sst <- -eof1.sst

eof.full <- rep(NA_real_, nlon * nlat)
eof.full[ok.cells] <- eof1.sst
eof.df <- grid.sst %>%
  mutate(EOF1 = eof.full) %>%
  filter(!is.na(EOF1))

# --- (d) Plot EOF1 map in same format as Section 9 regression ---
col.lim.eof <- max(abs(eof.df$EOF1), na.rm = TRUE)

p.sst.eof1 <- ggplot() +
  geom_raster(data = eof.df, aes(x = lon, y = lat, fill = EOF1)) +
  geom_polygon(data = mapWorld.clean,
               aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray30", linewidth = 0.25) +
  geom_path(data = al.box.df, aes(x = x, y = y),
            color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue4", mid = "white",
                       high = "darkorange4", midpoint = 0,
                       limits = c(-col.lim.eof, col.lim.eof),
                       name = "EOF1 loading") +
  coord_cartesian(xlim = c(180, 250), ylim = c(30, 66)) +
  scale_x_continuous(breaks = seq(180, 250, 20), labels = lon_label) +
  scale_y_continuous(breaks = seq(30, 60, 10),
                     labels = function(y) paste0(y, "\u00b0N")) +
  labs(title    = sprintf("ERA5 SST EOF1 (PDO-like): winter NDJFM, %d-%d",
                          sec19.win.range[1], sec19.win.range[2]),
       subtitle = sprintf("Leading mode of detrended monthly anomalies; explains ~%.1f%% of total variance",
                          100 * var.frac[1]),
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/SST_EOF1_PDO_pattern.png",
       plot = p.sst.eof1, width = 9, height = 6, dpi = 300)
message("Saved: Figures/SST_EOF1_PDO_pattern.png")

# --- (e) Match patterns on the original spatial grid (no aggregation) -
# Pair each native MLD-regression cell with its nearest SST-EOF cell
# (the two reanalyses live on slightly different native grids of
# similar resolution, so we use a nearest-neighbor join with a small
# tolerance rather than averaging into coarse blocks).
mld.native <- mld.beta %>%
  filter(is.finite(beta)) %>%
  select(lon, lat, beta)
eof.native <- eof.df %>%
  filter(is.finite(EOF1)) %>%
  select(lon, lat, EOF1)

# For each MLD cell, find nearest EOF cell.
eof.lons <- sort(unique(eof.native$lon))
eof.lats <- sort(unique(eof.native$lat))
nn.lon <- eof.lons[ pmax(1, pmin(length(eof.lons),
                                 findInterval(mld.native$lon, eof.lons,
                                              all.inside = TRUE))) ]
nn.lat <- eof.lats[ pmax(1, pmin(length(eof.lats),
                                 findInterval(mld.native$lat, eof.lats,
                                              all.inside = TRUE))) ]
# findInterval gives the floor; pick whichever neighbor is closer
nn.lon.up  <- eof.lons[ pmin(length(eof.lons),
                             findInterval(mld.native$lon, eof.lons,
                                          all.inside = TRUE) + 1) ]
nn.lat.up  <- eof.lats[ pmin(length(eof.lats),
                             findInterval(mld.native$lat, eof.lats,
                                          all.inside = TRUE) + 1) ]
nn.lon <- ifelse(abs(mld.native$lon - nn.lon.up) <
                 abs(mld.native$lon - nn.lon),  nn.lon.up, nn.lon)
nn.lat <- ifelse(abs(mld.native$lat - nn.lat.up) <
                 abs(mld.native$lat - nn.lat),  nn.lat.up, nn.lat)
match.tol <- 0.30   # degrees; native grids are ~0.25 deg
keep <- abs(mld.native$lon - nn.lon) <= match.tol &
        abs(mld.native$lat - nn.lat) <= match.tol
paired <- mld.native[keep, ] %>%
  mutate(eof.lon = nn.lon[keep], eof.lat = nn.lat[keep]) %>%
  inner_join(eof.native,
             by = c("eof.lon" = "lon", "eof.lat" = "lat")) %>%
  filter(is.finite(beta), is.finite(EOF1))
n.cells <- nrow(paired)
message(sprintf("Section 19: %d paired native-grid cells for coherence stats",
                n.cells))

# (1) Pearson r with spatial-block bootstrap p-value
r.obs <- cor(paired$beta, paired$EOF1)

set.seed(42)
n.boot <- 1000
block.size <- 5     # 5-deg lon/lat blocks (~500 km mid-lats)
paired$lon.blk <- block.size * floor(paired$lon / block.size)
paired$lat.blk <- block.size * floor(paired$lat / block.size)
paired$blk     <- paste(paired$lon.blk, paired$lat.blk, sep = "_")
blocks         <- unique(paired$blk)
n.blk          <- length(blocks)

r.boot <- replicate(n.boot, {
  shuf <- sample(blocks)
  lookup <- setNames(shuf, blocks)
  paired$blk2 <- lookup[paired$blk]
  shuffled <- paired %>%
    group_by(blk2) %>%
    mutate(EOF1.s = EOF1[sample(seq_len(n()))]) %>%
    ungroup()
  cor(paired$beta, shuffled$EOF1.s)
})
p.boot <- mean(abs(r.boot) >= abs(r.obs))
message(sprintf("(1) Pearson r(beta, EOF1) = %+.3f, n = %d, block-bootstrap p = %.3f",
                r.obs, n.cells, p.boot))

# (2) Tucker congruence coefficient (uncentered cosine)
tucker.phi <- sum(paired$beta * paired$EOF1) /
              sqrt(sum(paired$beta^2) * sum(paired$EOF1^2))
phi.label <- if (abs(tucker.phi) >= 0.95) "essentially identical"  else
             if (abs(tucker.phi) >= 0.85) "fair similarity"        else
                                          "patterns differ"
message(sprintf("(2) Tucker congruence phi = %+.3f (%s)",
                tucker.phi, phi.label))

# (3) Sign-agreement fraction; binomial test with spatial DOF
sign.agree <- mean(sign(paired$beta) == sign(paired$EOF1))
# Spatial DOF from EOF eigenvalue spread (Bretherton et al. 1999 eq. 31):
# N_eff ~= n.cells * (sum(lambda)^2 / sum(lambda^2)), capped at n.cells.
lambda.full <- (sv.sst$d^2)
n.eff <- min(n.cells,
             round(n.cells * (sum(lambda.full)^2 / sum(lambda.full^2))))
k.agree <- round(sign.agree * n.eff)
p.sign  <- binom.test(k.agree, n.eff, p = 0.5, alternative = "two.sided")$p.value
message(sprintf("(3) Sign agreement = %.1f%% (%d/%d cells), N_eff = %d, p = %.4f",
                100 * sign.agree, sum(sign(paired$beta) == sign(paired$EOF1)),
                n.cells, n.eff, p.sign))

# Save coherence metrics to CSV
coh.summary <- data.frame(
  metric    = c("pearson_r", "block_bootstrap_p", "tucker_phi",
                "sign_agreement_frac", "binomial_p_Neff"),
  value     = c(r.obs, p.boot, tucker.phi, sign.agree, p.sign),
  notes     = c(sprintf("n=%d cells (native grid, nearest-neighbor join)", n.cells),
                sprintf("%d block-shuffles, %d-deg blocks", n.boot, block.size),
                phi.label,
                sprintf("%d agree of %d cells", sum(sign(paired$beta) == sign(paired$EOF1)), n.cells),
                sprintf("Bretherton 1999 N_eff = %d", n.eff))
)
write.csv(coh.summary, "./Output/mld_ar1_vs_pdo_coherence.csv", row.names = FALSE)
message("Saved: Output/mld_ar1_vs_pdo_coherence.csv")

# --- (f) Side-by-side composite figure ---
mld.ar1.for.plot <- read.csv(mld.ar1.reg.cache) %>%
  filter(lon >= 180, lon <= 250, lat >= 30, lat <= 66) %>%
  mutate(lon = round(lon / 0.5) * 0.5,
         lat = round(lat / 0.5) * 0.5) %>%
  group_by(lon, lat) %>%
  summarise(beta = mean(beta, na.rm = TRUE), .groups = "drop") %>%
  complete(lon = seq(180, 250, by = 0.5),
           lat = seq(30,  66,  by = 0.5))
col.lim.ar1.r <- max(abs(mld.ar1.for.plot$beta), na.rm = TRUE)

p.mld.ar1.bare <- ggplot() +
  geom_raster(data = filter(mld.ar1.for.plot, !is.na(beta)),
              aes(x = lon, y = lat, fill = beta)) +
  geom_polygon(data = mapWorld.clean,
               aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray30", linewidth = 0.25) +
  geom_path(data = al.box.df, aes(x = x, y = y),
            color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_fill_gradient2(low = "steelblue4", mid = "white",
                       high = "darkorange4", midpoint = 0,
                       limits = c(-col.lim.ar1.r, col.lim.ar1.r),
                       name = "\u03b2 (AR(1)/z)") +
  coord_cartesian(xlim = c(180, 250), ylim = c(30, 66)) +
  scale_x_continuous(breaks = seq(180, 250, 20), labels = lon_label) +
  scale_y_continuous(breaks = seq(30, 60, 10),
                     labels = function(y) paste0(y, "\u00b0N")) +
  labs(title = "MLD AR(1) ~ AL SLP SD (Section 9)", x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3))

p.eof.bare <- p.sst.eof1 +
  labs(title = "ERA5 SST EOF1 (PDO-like)", subtitle = NULL)

p.compare <- p.mld.ar1.bare + p.eof.bare +
  plot_annotation(
    title    = "MLD AR(1) regression on AL SLP SD vs PDO spatial pattern",
    subtitle = sprintf("Pearson r = %+.2f (block-bootstrap p = %.3f) | Tucker \u03c6 = %+.2f | sign agreement = %.0f%% (p = %.3f)",
                       r.obs, p.boot, tucker.phi, 100 * sign.agree, p.sign)
  )

ggsave("./Figures/MLD_AR1_vs_PDO_pattern_comparison.png",
       plot = p.compare, width = 14, height = 6, dpi = 300)
message("Saved: Figures/MLD_AR1_vs_PDO_pattern_comparison.png")
message("===================================================\n")

# ============================================================
# SECTION 20: AR(1) effect on extreme event risk (simulation)
# ============================================================
# Motivation: how does decadal-scale change in winter SST AR(1)
# alter the probability of extreme heatwaves and cold spells?
#
# Workflow:
#   (a) Fit a linear trend to GOA and EBS winter (Nov-Mar) mean SST
#       anomaly (winter.anom from Section 2). Plot with the linear
#       trend overlaid; save trend + SD per system.
#   (b) Simulate 1000 winter time series per system under three
#       AR(1) conditions: min observed (Sec 3), 0, max observed.
#       Marginal SD of simulated series = observed SD; innovation
#       SD = obs.sd * sqrt(1 - rho^2). Mean trajectory = fitted
#       linear trend.
#   (c) 1950-1979 climatology (mean + SD) of the observed winter
#       anomaly. Heatwave threshold = clim.mean + 2*clim.sd;
#       cold-spell threshold = clim.mean - 2*clim.sd.
#   (d) For each year, compute % of 1000 sims above heatwave
#       threshold and below cold-spell threshold.
#   (e) Four-panel plot: 2 systems x {heatwave %, cold spell %};
#       color-coded by AR(1) condition.

# --- (a) Linear trends + SDs (saved for downstream use) -------

trend.tbl <- bind_rows(
  data.frame(region = "GOA", winter.anom %>% select(year, anom = GOA)),
  data.frame(region = "EBS", winter.anom %>% select(year, anom = EBS))
) %>%
  filter(!is.na(anom))

trend.fits <- trend.tbl %>%
  group_by(region) %>%
  summarise(intercept = coef(lm(anom ~ year))[1],
            slope     = coef(lm(anom ~ year))[2],
            obs.sd    = sd(anom),
            .groups   = "drop")

write.csv(trend.fits, "./Output/winter_SST_anom_trend_sd.csv", row.names = FALSE)
message("Saved: Output/winter_SST_anom_trend_sd.csv")
message("Linear trend & SD per system:")
print(trend.fits)

p.win.lin <- trend.tbl %>%
  mutate(region = factor(region, levels = c("GOA", "EBS"))) %>%
  ggplot(aes(x = year, y = anom, color = region)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              color = "black", linewidth = 0.8) +
  scale_color_manual(values = region.colors, labels = region.labels) +
  facet_wrap(~ region, ncol = 1, scales = "free_y",
             labeller = labeller(region = region.labels)) +
  scale_x_continuous(breaks = seq(1950, end.yr, 10)) +
  labs(title    = paste0("Winter (Nov\u2013Mar) SST Anomaly with Linear Trend, 1950\u2013", end.yr),
       subtitle = "Black line = OLS linear fit; trend + SD saved for AR(1) simulations",
       x = "Year", y = "SST Anomaly (\u00b0C)") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "none",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/ERA5_SST_winter_anom_linear.png",
       plot = p.win.lin, width = 7, height = 6, dpi = 300)
message("Saved: Figures/ERA5_SST_winter_anom_linear.png")

# --- (b) AR(1) extremes per system (min, ~0, max from Section 3 rolling)
# Marginal SD is held fixed at the full time series SD (trend.fits$obs.sd)
# so the conditions differ ONLY in rho. This isolates the AR(1) effect:
# AR(1) does not change the marginal per-year exceedance probability,
# but it does change the probability of CONSECUTIVE exceedance years.

ar1.range <- winter.roll %>%
  summarise(GOA.min = min(GOA.ar1, na.rm = TRUE),
            GOA.max = max(GOA.ar1, na.rm = TRUE),
            EBS.min = min(EBS.ar1, na.rm = TRUE),
            EBS.max = max(EBS.ar1, na.rm = TRUE))
message("Observed 15-yr rolling AR(1) extremes:")
print(ar1.range)

ar1.conditions <- list(
  GOA = c(min = ar1.range$GOA.min, max = ar1.range$GOA.max),
  EBS = c(min = ar1.range$EBS.min, max = ar1.range$EBS.max)
)

# --- (c) 1950-1979 climatology (mean + SD) per system --------

clim.tbl <- trend.tbl %>%
  filter(year >= 1950, year <= 1979) %>%
  group_by(region) %>%
  summarise(clim.mean = mean(anom),
            clim.sd   = sd(anom),
            .groups   = "drop")
message("1950-1979 climatology per system:")
print(clim.tbl)

# --- (d) Simulations -----------------------------------------

set.seed(42)
n.sim   <- 10000
sim.years <- sort(unique(trend.tbl$year))
n.t     <- length(sim.years)

simulate_ar1 <- function(rho, n, sd.marginal) {
  # AR(1) with marginal SD = sd.marginal: innovation sd = sd.marg * sqrt(1 - rho^2)
  if (abs(rho) >= 1) stop("|rho| must be < 1")
  sd.inn <- sd.marginal * sqrt(1 - rho^2)
  x <- numeric(n)
  x[1] <- rnorm(1, 0, sd.marginal)
  for (i in 2:n) x[i] <- rho * x[i - 1] + rnorm(1, 0, sd.inn)
  x
}

results <- list()
for (sys in c("GOA", "EBS")) {
  fit.s  <- trend.fits %>% filter(region == sys)
  cl.s   <- clim.tbl   %>% filter(region == sys)
  trend.line <- fit.s$intercept + fit.s$slope * sim.years
  hi.thr <- cl.s$clim.mean + 2 * cl.s$clim.sd
  lo.thr <- cl.s$clim.mean - 2 * cl.s$clim.sd
  sd.fix <- fit.s$obs.sd

  for (cond in c("min", "max")) {
    rho <- ar1.conditions[[sys]][[cond]]
    # n.t x n.sim matrix of simulated anomalies (trend + AR(1) noise)
    sim.mat <- replicate(n.sim,
      trend.line + simulate_ar1(rho = rho, n = n.t,
                                sd.marginal = sd.fix))
    above <- sim.mat > hi.thr   # n.t x n.sim logical
    below <- sim.mat < lo.thr
    # Probability that year y is the SECOND of >= 2 consecutive
    # exceedance years: indicator that years (y-1, y) both exceed.
    pair.hi <- above[-1, ] & above[-n.t, ]   # (n.t-1) x n.sim
    pair.lo <- below[-1, ] & below[-n.t, ]
    n.hi <- c(NA, rowSums(pair.hi))
    n.lo <- c(NA, rowSums(pair.lo))
    pct.hi <- 100 * n.hi / n.sim
    pct.lo <- 100 * n.lo / n.sim
    # Wilson 95% CI on each binomial proportion (n.sim Bernoulli trials)
    wilson <- function(k, n, conf = 0.95) {
      z <- qnorm(1 - (1 - conf) / 2)
      p <- k / n
      denom <- 1 + z^2 / n
      ctr   <- (p + z^2 / (2 * n)) / denom
      hw    <- (z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / denom
      list(lo = 100 * (ctr - hw), hi = 100 * (ctr + hw))
    }
    ci.hi <- wilson(n.hi, n.sim)
    ci.lo <- wilson(n.lo, n.sim)
    results[[paste(sys, cond, sep = "_")]] <- data.frame(
      region    = sys,
      condition = cond,
      rho       = as.numeric(rho),
      sd.fixed  = sd.fix,
      year      = sim.years,
      pct.hot      = pct.hi,
      pct.hot.lo   = ci.hi$lo,
      pct.hot.hi   = ci.hi$hi,
      pct.cold     = pct.lo,
      pct.cold.lo  = ci.lo$lo,
      pct.cold.hi  = ci.lo$hi
    )
  }
}
sim.out <- bind_rows(results)
write.csv(sim.out, "./Output/AR1_extreme_event_simulation.csv", row.names = FALSE)
message("Saved: Output/AR1_extreme_event_simulation.csv")

# --- (e) Four-panel plot: 2 systems x {heatwave %, cold spell %}

# Colorblind palette (Wong 2011) — pick min = blue (#0072B2), max = vermillion (#D55E00)
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cond.labels <- c(min = "AR(1) = min observed",
                 max = "AR(1) = max observed")
cond.colors <- c(min = cb[6], max = cb[7])

extreme.levels <- c("% \u2265 2 consecutive heatwave years (> +2 SD)",
                    "% \u2265 2 consecutive cold-spell years (< -2 SD)")

sim.long <- bind_rows(
  sim.out %>%
    transmute(region, condition, year,
              extreme = extreme.levels[1],
              pct = pct.hot, pct.lo = pct.hot.lo, pct.hi = pct.hot.hi),
  sim.out %>%
    transmute(region, condition, year,
              extreme = extreme.levels[2],
              pct = pct.cold, pct.lo = pct.cold.lo, pct.hi = pct.cold.hi)
) %>%
  mutate(extreme   = factor(extreme, levels = extreme.levels),
         region    = factor(region, levels = c("GOA", "EBS")),
         condition = factor(condition, levels = c("min", "max")))

p.extreme <- ggplot(sim.long,
                    aes(x = year, y = pct,
                        color = condition, fill = condition)) +
  geom_ribbon(aes(ymin = pct.lo, ymax = pct.hi),
              alpha = 0.25, color = NA) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = cond.colors, labels = cond.labels,
                     name = NULL) +
  scale_fill_manual(values = cond.colors, labels = cond.labels,
                    name = NULL) +
  facet_grid(extreme ~ region, scales = "free_y",
             labeller = labeller(region = region.labels)) +
  scale_x_continuous(breaks = seq(1950, end.yr, 10)) +
  labs(title    = sprintf("Effect of AR(1) on consecutive-extreme-year risk (%s simulations)",
                          format(n.sim, big.mark = ",")),
       subtitle = "Mean = OLS trend; marginal SD = full-series SD (held fixed across AR(1) conditions); thresholds = 1950-1979 mean +/- 2 SD; ribbons = Wilson 95% CI",
       x = "Year",
       y = "% of simulations with year y and y-1 both exceeding threshold") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        plot.subtitle    = element_text(size = 9))

ggsave("./Figures/AR1_extreme_event_simulation.png",
       plot = p.extreme, width = 10, height = 7, dpi = 300)
message("Saved: Figures/AR1_extreme_event_simulation.png")
message("===================================================\n")
