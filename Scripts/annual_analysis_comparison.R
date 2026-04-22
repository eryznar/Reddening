# PURPOSE: Reddening project — annual-means counterpart to analysis.R.
# Mirrors analysis.R but replaces every winter (Nov-Mar, January-year)
# aggregation with calendar-year annual means, so results can be compared
# to the winter pipeline. All figure/output filenames are suffixed
# '_annual' to keep the two sets of outputs side-by-side.
#
# Omitted vs. analysis.R:
#   * Reference map (Section 1 of analysis.R) — not time-averaged,
#     identical regardless of winter/annual; see NP_reference_map.png.
#   * SLP EOF1 loadings map — not time-averaged at the plotting step
#     (fit to all months), identical to ERA5_SLP_EOF1_loadings.png.
#     The EOF computation is kept because Section 5 needs al.box.mask
#     and detrend.slp.
#
# Sections:
#   2. Annual SST anomaly time series
#   3. 15-yr rolling AR(1) and SD — annual
#   4. SLP EOF1 computation (no map output)
#   5. Annual AL SLP SD vs SST AR(1) — 15-yr dual-axis plot

library(ggplot2)
library(dplyr)
library(maps)
library(ncdf4)
library(zoo)
library(tidyr)
library(sf)
library(patchwork)
library(irlba)

# ============================================================
# SECTION 2: ERA5 SST ANNUAL ANOMALY TIME SERIES
# ============================================================
# Monthly anomalies (1950-1979 climatology) → calendar-year means for
# GOA, EBS, North Pacific. Only years with all 12 months present are
# retained; partial years (e.g. 2026 Jan-Mar only) are dropped.

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

# Study-area polygons (same vertices as analysis.R)
goa.x <- c(191, 191, 203, 205, 208, 223, 234)
goa.y <- c(50,  53,  57.5, 59,  61,  61,  50)
ebs.x <- c(183, 183, 203, 203, 191)
ebs.y <- c(53,  65,  65,  57.5, 53)

to180 <- function(x) ifelse(x > 180, x - 360, x)

make_sf_poly <- function(lon360, lat) {
  coords <- rbind(cbind(to180(lon360), lat), c(to180(lon360[1]), lat[1]))
  st_sf(geometry = st_sfc(st_polygon(list(coords)), crs = 4326))
}

goa.sf <- make_sf_poly(goa.x, goa.y)
ebs.sf <- make_sf_poly(ebs.x, ebs.y)

grid    <- expand.grid(lon = lons, lat = lats)
grid.sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)

in.goa <- lengths(st_within(grid.sf, goa.sf)) > 0
in.ebs <- lengths(st_within(grid.sf, ebs.sf)) > 0
in.np  <- !is.na(as.vector(sst[,,1]))

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

# Monthly climatology (1950-1979)
clim.sst <- monthly.sst %>%
  filter(year >= 1950, year <= 1979) %>%
  group_by(month) %>%
  summarise(GOA.clim = mean(GOA, na.rm = TRUE),
            EBS.clim = mean(EBS, na.rm = TRUE),
            NP.clim  = mean(NP,  na.rm = TRUE),
            .groups = "drop")

monthly.anom <- monthly.sst %>%
  left_join(clim.sst, by = "month") %>%
  mutate(GOA.anom = GOA - GOA.clim,
         EBS.anom = EBS - EBS.clim,
         NP.anom  = NP  - NP.clim) %>%
  select(year, month, GOA.anom, EBS.anom, NP.anom)

# --- Annual means: calendar year, all 12 months present ---
annual.anom <- monthly.anom %>%
  group_by(year) %>%
  summarise(GOA = mean(GOA.anom, na.rm = TRUE),
            EBS = mean(EBS.anom, na.rm = TRUE),
            NP  = mean(NP.anom,  na.rm = TRUE),
            n.mo = n(),
            .groups = "drop") %>%
  filter(n.mo == 12) %>%
  select(-n.mo) %>%
  arrange(year)

region.colors <- c("GOA" = "darkorange4", "EBS" = "steelblue4", "NP" = "forestgreen")
region.labels <- c("GOA" = "Gulf of Alaska", "EBS" = "Eastern Bering Sea",
                   "NP"  = "North Pacific")

end.yr <- max(annual.anom$year, na.rm = TRUE)

p.ann.anom <- annual.anom %>%
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
  labs(title   = paste0("Annual SST Anomaly, 1950\u2013", end.yr),
       x = "Year", y = "SST Anomaly (\u00b0C)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"))

ggsave("./Figures/ERA5_SST_annual_anom_ts.png",
       plot = p.ann.anom, width = 7, height = 8, dpi = 300)
message("Saved: Figures/ERA5_SST_annual_anom_ts.png")

# ============================================================
# SECTION 3: 15-YEAR ROLLING AR(1) AND SD — annual, GOA and EBS
# ============================================================

detrend <- function(x) residuals(lm(x ~ seq_along(x)))

roll_ar1 <- function(x, width = 15) {
  rollapply(x, width = width, fill = NA, align = "right",
            FUN = function(v) acf(v, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2])
}
roll_sd <- function(x, width = 15) {
  rollapply(x, width = width, fill = NA, align = "right", FUN = sd)
}

annual.detr <- annual.anom %>%
  mutate(GOA = detrend(GOA),
         EBS = detrend(EBS),
         NP  = detrend(NP))

annual.roll <- annual.detr %>%
  mutate(GOA.ar1 = roll_ar1(GOA),
         EBS.ar1 = roll_ar1(EBS),
         GOA.sd  = roll_sd(GOA),
         EBS.sd  = roll_sd(EBS))

p.goa.ar1 <- ggplot(annual.roll, aes(x = year, y = GOA.ar1)) +
  geom_line(color = "darkorange4") +
  geom_point(color = "darkorange4", size = 1.4) +
  labs(title = "Gulf of Alaska", x = NULL, y = "AR(1)") +
  theme_bw(base_size = 11)

p.ebs.ar1 <- ggplot(annual.roll, aes(x = year, y = EBS.ar1)) +
  geom_line(color = "steelblue4") +
  geom_point(color = "steelblue4", size = 1.4) +
  labs(title = "Eastern Bering Sea", x = NULL, y = "AR(1)") +
  theme_bw(base_size = 11)

p.goa.sd <- ggplot(annual.roll, aes(x = year, y = GOA.sd)) +
  geom_line(color = "darkorange4") +
  geom_point(color = "darkorange4", size = 1.4) +
  labs(x = "Year", y = "SD (\u00b0C)") +
  theme_bw(base_size = 11)

p.ebs.sd <- ggplot(annual.roll, aes(x = year, y = EBS.sd)) +
  geom_line(color = "steelblue4") +
  geom_point(color = "steelblue4", size = 1.4) +
  labs(x = "Year", y = "SD (\u00b0C)") +
  theme_bw(base_size = 11)

p.roll <- (p.goa.ar1 | p.ebs.ar1) / (p.goa.sd | p.ebs.sd) +
  plot_annotation(
    title   = "15-yr Rolling Annual SST AR(1) and SD (right-aligned)",
    theme   = theme(plot.title = element_text(size = 12, hjust = 0.5))
  )

ggsave("./Figures/ERA5_SST_AR1_SD_15yr_rolling_annual.png",
       plot = p.roll, width = 10, height = 6.5, dpi = 300)
message("Saved: Figures/ERA5_SST_AR1_SD_15yr_rolling_annual.png")

# ============================================================
# SECTION 4: SLP EOF1 COMPUTATION (no map output)
# ============================================================
# Identical to analysis.R Section 4 but no plot is saved — loadings map
# is the same under winter/annual workflows. Result: detrend.slp (cell ×
# time matrix of detrended monthly anomalies) and al.box.mask, both
# needed for Section 5.

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

# Monthly anomalies (1950-1979 climatology)
slp.anom <- array(NA_real_, dim = dim(slp))
base.slp <- years.slp >= 1950 & years.slp <= 1979

for (m in 1:12) {
  t.base <- which(base.slp & months.slp == m)
  t.all  <- which(months.slp == m)
  clim.m <- apply(slp[, , t.base, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  for (t in t.all) slp.anom[, , t] <- slp[, , t] - clim.m
}

# Detrend each cell
X.slp    <- cbind(1, seq_len(n.time.slp))
IH.slp   <- diag(n.time.slp) - X.slp %*% solve(t(X.slp) %*% X.slp) %*% t(X.slp)
anom.mat.slp  <- matrix(slp.anom,
                        nrow = length(lons.slp) * length(lats.slp),
                        ncol = n.time.slp)
detrend.slp <- anom.mat.slp %*% t(IH.slp)

# Area-weighted truncated SVD for EOF1
grid.slp <- expand.grid(lon = lons.slp, lat = lats.slp)
w.slp    <- sqrt(cos(grid.slp$lat * pi / 180))
sv.slp   <- irlba(sweep(detrend.slp, 1, w.slp, "*"), nv = 1)
eof1.slp <- sv.slp$u[, 1] / w.slp

# AL box (Litzow et al. 2020 PNAS: 192.5-207.5E, 45-55N)
al.box.mask <- grid.slp$lon >= 192.5 & grid.slp$lon <= 207.5 &
               grid.slp$lat >=  45   & grid.slp$lat <=  55

# (No sign flip or plot — detrend.slp is what Section 5 needs, and SD
#  is sign-invariant.)

# ============================================================
# SECTION 5: ANNUAL AL SLP SD vs SST AR(1) — 15-yr dual-axis plot
# ============================================================
# Calendar-year AL SLP anomaly (mean over detrended monthly anomalies
# inside al.box.mask, all 12 months). 15-yr rolling SD (z-scored),
# joined with annual SST AR(1) from Section 3.

al.monthly <- data.frame(
  year  = years.slp,
  month = months.slp,
  SLP   = apply(detrend.slp[al.box.mask, ], 2, mean, na.rm = TRUE)
)

al <- al.monthly %>%
  group_by(year) %>%
  summarise(SLP = mean(SLP, na.rm = TRUE), n.mo = n(), .groups = "drop") %>%
  filter(n.mo == 12) %>%
  select(-n.mo) %>%
  arrange(year) %>%
  mutate(AL.sd = roll_sd(SLP),
         AL.sd = as.numeric(scale(AL.sd)))

al.sst <- annual.roll %>%
  select(year, GOA.ar1, EBS.ar1) %>%
  left_join(al %>% select(year, AL.sd), by = "year") %>%
  filter(!is.na(AL.sd))

al.sd.range   <- range(al.sst$AL.sd,   na.rm = TRUE)
goa.ar1.range <- range(al.sst$GOA.ar1, na.rm = TRUE)
ebs.ar1.range <- range(al.sst$EBS.ar1, na.rm = TRUE)

# --- Modified Chelton method (Pyper & Peterman 1998) ---
chelton_test <- function(x, y) {
  n       <- length(x)
  r.obs   <- cor(x, y, use = "complete.obs")
  max.lag <- floor(n / 5)

  acf.x <- acf(x, lag.max = max.lag, plot = FALSE, na.action = na.pass)$acf[-1]
  acf.y <- acf(y, lag.max = max.lag, plot = FALSE, na.action = na.pass)$acf[-1]

  j     <- seq_len(max.lag)
  n.eff <- n / (1 + 2 * sum((1 - j / n) * acf.x * acf.y))
  n.eff <- max(3, min(n, n.eff))

  t.obs <- r.obs * sqrt(n.eff - 2) / sqrt(1 - r.obs^2)
  df    <- n.eff - 2
  pval  <- pt(t.obs, df = df, lower.tail = FALSE)

  list(r.obs = r.obs, t.obs = t.obs, n.eff = n.eff, df = df, pval = pval)
}

dat.goa  <- al.sst %>% filter(!is.na(GOA.ar1), !is.na(AL.sd))
dat.ebs  <- al.sst %>% filter(!is.na(EBS.ar1), !is.na(AL.sd))

chel.goa <- chelton_test(dat.goa$GOA.ar1, dat.goa$AL.sd)
chel.ebs <- chelton_test(dat.ebs$EBS.ar1, dat.ebs$AL.sd)

message(sprintf("GOA Chelton (annual): r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
                chel.goa$r.obs, chel.goa$n.eff, chel.goa$t.obs, chel.goa$pval))
message(sprintf("EBS Chelton (annual): r = %.2f, N* = %.1f, t = %.2f, p = %.4f",
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
    title = "15-yr Rolling AL SLP SD vs. SST AR(1) \u2014 Annual",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 12))
  )

ggsave("./Figures/AL_SD_SST_AR1_15yr_dual_axis_annual.png",
       plot = p.dual, width = 7, height = 7, dpi = 300)
message("Saved: Figures/AL_SD_SST_AR1_15yr_dual_axis_annual.png")

make_chelton_panel <- function(t.obs, df, n.eff, pval, region.label) {
  t.range <- seq(-4, max(4, t.obs + 0.5), length.out = 500)
  lbl <- sprintf("p = %.4f\nN* = %.1f", pval, n.eff)
  ggplot(data.frame(t = t.range, d = dt(t.range, df = df)), aes(x = t, y = d)) +
    geom_area(fill = "gray80", color = "gray50", linewidth = 0.3) +
    geom_vline(xintercept = t.obs, linetype = "dashed", linewidth = 0.9) +
    annotate("text", x = Inf, y = Inf, label = lbl,
             hjust = 1.1, vjust = 1.5, size = 3.5, lineheight = 1.3) +
    labs(title = paste(region.label, "\u2014 Chelton method (annual)"),
         x = sprintf("t  (df = %.1f)", df), y = "Density") +
    theme_bw(base_size = 11)
}

p.goa.chel <- make_chelton_panel(chel.goa$t.obs, chel.goa$df, chel.goa$n.eff, chel.goa$pval, "GOA")
p.ebs.chel <- make_chelton_panel(chel.ebs$t.obs, chel.ebs$df, chel.ebs$n.eff, chel.ebs$pval, "EBS")

p.chelton <- p.goa.chel / p.ebs.chel +
  plot_annotation(
    title = "Chelton significance test: AL SLP SD vs. SST AR(1) \u2014 Annual",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 12))
  )

ggsave("./Figures/AL_SD_SST_AR1_chelton_dist_annual.png",
       plot = p.chelton, width = 5, height = 7, dpi = 300)
message("Saved: Figures/AL_SD_SST_AR1_chelton_dist_annual.png")
