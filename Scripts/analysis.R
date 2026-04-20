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

area.colors    <- c("GOA" = "steelblue4",
                    "EBS" = "darkorange4")

area.linetypes <- c("GOA" = "solid",
                    "EBS" = "solid")

area.linewidths <- c("GOA" = 0.9,
                     "EBS" = 0.9)

area.labels <- c(
  "GOA" = "Gulf of Alaska (GOA)",
  "EBS" = "Eastern Bering Sea (EBS)"
)

region.fills <- c("GOA" = "steelblue", "EBS" = "darkorange")

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
region.colors <- c("GOA" = "steelblue4", "EBS" = "darkorange4", "NP" = "forestgreen")
region.labels <- c("GOA" = "Gulf of Alaska", "EBS" = "Eastern Bering Sea",
                   "NP"  = "North Pacific")

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
  geom_line(color = "steelblue4") +
  geom_point(color = "steelblue4", size = 1.4) +
  labs(title = "Gulf of Alaska", x = NULL, y = "AR(1)") +
  theme_bw(base_size = 11)

p.ebs.ar1 <- ggplot(winter.roll, aes(x = year, y = EBS.ar1)) +
  geom_line(color = "darkorange4") +
  geom_point(color = "darkorange4", size = 1.4) +
  labs(title = "Eastern Bering Sea", x = NULL, y = "AR(1)") +
  theme_bw(base_size = 11)

p.goa.sd <- ggplot(winter.roll, aes(x = year, y = GOA.sd)) +
  geom_line(color = "steelblue4") +
  geom_point(color = "steelblue4", size = 1.4) +
  labs(x = "Year", y = "SD (\u00b0C)") +
  theme_bw(base_size = 11)

p.ebs.sd <- ggplot(winter.roll, aes(x = year, y = EBS.sd)) +
  geom_line(color = "darkorange4") +
  geom_point(color = "darkorange4", size = 1.4) +
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
# SECTION 4: SLP EOF1 LOADINGS MAP — MERCATOR PROJECTION
# ============================================================
# Reproduces the EOF1 analysis from ERA5_SLP_download.R and plots
# the spatial loadings on a Mercator projection, with the Aleutian
# Low box (191-208E, 44-55N) overlaid. Sign is flipped if needed
# so the AL box region holds positive loadings.

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

# Flip sign so AL box (191-208E, 44-55N) has positive mean loading.
# grid.slp$lon is in 0-360 space (ERA5 returns 120-250E for this domain).
al.box.mask <- grid.slp$lon >= 191 & grid.slp$lon <= 208 &
               grid.slp$lat >=  44 & grid.slp$lat <=  55
al.box.mean <- mean(eof1.slp[al.box.mask], na.rm = TRUE)
if (!is.na(al.box.mean) && al.box.mean < 0) eof1.slp <- -eof1.slp

# Convert lons to 0-360 for Pacific-centered plot
eof.df.slp <- grid.slp %>%
  mutate(lon  = ifelse(lon < 0, lon + 360, lon),
         EOF1 = eof1.slp)

# AL box in 0-360 space
al.box.df <- data.frame(
  x = c(191, 208, 208, 191, 191),
  y = c(44,  44,  55,  55,  44)
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
# Loads the pre-computed AL winter SLP anomaly, calculates 15-yr
# rolling SD (z-scored), joins with SST AR(1) from Section 3,
# and plots dual-axis time series for GOA and EBS.

al <- read.csv("./Output/AL_winter_SLP_anomaly.csv") %>%
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

message(sprintf("GOA: r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
                chel.goa$r.obs, chel.goa$n.eff, chel.goa$t.obs, chel.goa$pval))
message(sprintf("EBS: r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
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

p.goa.dual <- make_dual_axis_plot(al.sst, "GOA.ar1", "steelblue4", "GOA", goa.ar1.range, chel.goa)
p.ebs.dual <- make_dual_axis_plot(al.sst, "EBS.ar1", "darkorange4", "EBS", ebs.ar1.range, chel.ebs)

p.dual <- p.goa.dual / p.ebs.dual +
  plot_annotation(
    title = "15-yr Rolling AL SLP SD vs. SST AR(1) \u2014 Winter (Nov\u2013Mar)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 12))
  )

ggsave("./Figures/AL_SD_SST_AR1_15yr_dual_axis.png",
       plot = p.dual, width = 7, height = 7, dpi = 300)
message("Saved: Figures/AL_SD_SST_AR1_15yr_dual_axis.png")

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

chel.goa <- chelton_test(dat.goa$GOA.ar1, dat.goa$AL.sd)
chel.ebs <- chelton_test(dat.ebs$EBS.ar1, dat.ebs$AL.sd)

message(sprintf("GOA Chelton: r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
                chel.goa$r.obs, chel.goa$n.eff, chel.goa$t.obs, chel.goa$pval))
message(sprintf("EBS Chelton: r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
                chel.ebs$r.obs, chel.ebs$n.eff, chel.ebs$t.obs, chel.ebs$pval))

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

