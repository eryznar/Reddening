# PURPOSE: Compare observed and modeled AR1 values for SST in the GOA and EBS
# Author: Mike Litzow

# LOAD LIBS/FUNCTIONS ----------------------------------
source("./Scripts/load.libs.functions.R")

# MAP OF GOA AND EBS POLYGONS ----------------------------------

# Pacific-centered world map (longitudes wrapped to 20-380 range)
mapWorld <- map_data("world", wrap = c(20, 380))

# Define GOA polygon in 0-360 space (for masking, convert to -180/180 below)
goa.x <- c(198, 198, 203, 205, 208, 225, 231, 201)
goa.y <- c(54, 55.5, 57.5, 59, 61, 61, 54, 54)

goa.poly <- data.frame(lon = goa.x, lat = goa.y, region = "GOA")

# Define EBS polygon in 0-360 space
ebs.x <- c(183, 183, 203, 203, 191)
ebs.y <- c(53, 65, 65, 57.5, 53)

ebs.poly <- data.frame(lon = ebs.x, lat = ebs.y, region = "EBS")

# Convert to -180/180 for sf masking
goa.x <- ifelse(goa.x > 180, goa.x - 360, goa.x)
ebs.x <- ifelse(ebs.x > 180, ebs.x - 360, ebs.x)

# Combine GOA and EBS polygons (0-360 lons for Pacific-centered plot)
polys <- bind_rows(goa.poly, ebs.poly)

# North Pacific ERA5 domain (20-66N, 110-250E) — single rectangle in 0-360 space
np.rect <- data.frame(xmin = 110, xmax = 250, ymin = 20, ymax = 66)

# Plot
map.plot <- ggplot() +
  geom_polygon(data = mapWorld,
               aes(x = long, y = lat, group = group),
               fill = "gray80", color = "gray50", linewidth = 0.2) +
  geom_rect(data = np.rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "forestgreen", color = "forestgreen", alpha = 0.1, linewidth = 0.7) +
  geom_polygon(data = polys,
               aes(x = lon, y = lat, group = region, fill = region, color = region),
               alpha = 0.4, linewidth = 1) +
  scale_fill_manual(values = c("GOA" = "steelblue", "EBS" = "darkorange")) +
  scale_color_manual(values = c("GOA" = "steelblue4", "EBS" = "darkorange4")) +
  coord_cartesian(xlim = c(100, 260), ylim = c(15, 70)) +
  labs(title = "North Pacific ERA5 Domain with GOA and EBS Study Regions",
       x = "Longitude", y = "Latitude",
       fill = "Region", color = "Region") +
  theme_bw() +
  theme(legend.position = "bottom")

print(map.plot)

# DOWNLOAD ERA5 SST ----------------------------------

library(ecmwfr)

# Set CDS API key (run once, then comment out — stores credentials in macOS Keychain)
# Find your UID and API key at: https://cds.climate.copernicus.eu -> your profile
# wf_set_key(user = "YOUR_CDS_UID", key = "YOUR_CDS_API_KEY", service = "cds")

# Full North Pacific domain: 20-66N, 110-250E (area: N, W, S, E in -180/180)
request <- list(
  dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
  product_type    = "monthly_averaged_reanalysis",
  variable        = "sea_surface_temperature",
  year            = as.character(1950:2025),
  month           = sprintf("%02d", 1:12),
  time            = "00:00",
  area            = c(66, 110, 20, -110),    # N, W, S, E (250E = -110W)
  data_format     = "netcdf",
  download_format = "unarchived",
  target          = "era5_sst_NP_1950_2025.nc"
)

# Slow down polling to avoid CDS rate limit
options(ecmwfr.sleep = 120)   # wait 120s between status checks

wf_request(
  request  = request,
  transfer = TRUE,
  path     = "./Data"
)

# CALCULATE OBSERVED AR1 AND SD FROM ERA5 ----------------------------------

# Load ERA5 SST netCDF (North Pacific domain)
nc <- nc_open("./Data/era5_sst_NP_1950_2025.nc")

sst  <- ncvar_get(nc, "sst")        # [longitude x latitude x valid_time], units: K
lons <- ncvar_get(nc, "longitude")
lats <- ncvar_get(nc, "latitude")
time <- ncvar_get(nc, "valid_time") # seconds since 1970-01-01
nc_close(nc)

# Convert time to dates
dates  <- as.Date(as.POSIXct(time, origin = "1970-01-01", tz = "UTC"))
months <- as.integer(format(dates, "%m"))
years  <- as.integer(format(dates, "%Y"))

# Convert SST from Kelvin to Celsius
sst <- sst - 273.15

# Build polygon sf objects for masking (close rings by repeating first coordinate)
goa.coords <- rbind(cbind(goa.x, goa.y), c(goa.x[1], goa.y[1]))
ebs.coords <- rbind(cbind(ebs.x, ebs.y), c(ebs.x[1], ebs.y[1]))

goa.sf <- st_sf(geometry = st_sfc(st_polygon(list(goa.coords)), crs = 4326))
ebs.sf <- st_sf(geometry = st_sfc(st_polygon(list(ebs.coords)), crs = 4326))

# Build grid of all lon/lat combinations matching array order
grid    <- expand.grid(lon = lons, lat = lats)
grid.sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)

# Identify which cells fall inside each polygon / domain
in.goa <- lengths(st_within(grid.sf, goa.sf)) > 0
in.ebs <- lengths(st_within(grid.sf, ebs.sf)) > 0
in.np  <- !is.na(as.vector(sst[,,1]))   # all ocean cells

# Latitude-based weights (cos(lat) proportional to cell area)
weights <- cos(grid$lat * pi / 180)

# --- Step 1: area-wide weighted mean SST for each domain ---
weighted_region_mean <- function(mask) {
  w <- weights[mask]
  apply(sst, 3, function(slice) {
    vals <- as.vector(slice)[mask]
    if (all(is.na(vals))) return(NA_real_)
    weighted.mean(vals, w, na.rm = TRUE)
  })
}

goa.sst <- weighted_region_mean(in.goa)
ebs.sst <- weighted_region_mean(in.ebs)
np.sst  <- weighted_region_mean(in.np)

# Assemble monthly data frame
monthly <- data.frame(
  year  = years,
  month = months,
  GOA   = goa.sst,
  EBS   = ebs.sst,
  NP    = np.sst
)

# --- Step 2: monthly climatology for each domain (1950-1979) ---
clim <- monthly %>%
  filter(year >= 1950, year <= 1979) %>%
  group_by(month) %>%
  summarise(GOA.clim = mean(GOA, na.rm = TRUE),
            EBS.clim = mean(EBS, na.rm = TRUE),
            NP.clim  = mean(NP,  na.rm = TRUE),
            .groups = "drop")

# --- Step 3: anomaly time series (°C) for full period ---
monthly <- monthly %>%
  left_join(clim, by = "month") %>%
  mutate(GOA.anom = GOA - GOA.clim,
         EBS.anom = EBS - EBS.clim,
         NP.anom  = NP  - NP.clim)

# Save non-detrended anomalies for plotting
monthly.anom <- monthly %>%
  mutate(date = as.Date(paste(year, month, "15", sep = "-"))) %>%
  select(date, year, month, GOA.anom, EBS.anom, NP.anom)

# --- Step 4: detrend each domain's anomaly time series ---
detrend <- function(x) residuals(lm(x ~ seq_along(x)))

monthly <- monthly %>%
  mutate(GOA.anom = detrend(GOA.anom),
         EBS.anom = detrend(EBS.anom),
         NP.anom  = detrend(NP.anom))

# --- Step 5: restrict to winter (Nov-Mar), year = January year ---
monthly.win <- monthly %>%
  filter(month %in% c(11, 12, 1, 2, 3)) %>%
  mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year))

# Winter mean anomaly per domain per year
winter <- monthly.win %>%
  group_by(win.year) %>%
  summarise(GOA = mean(GOA.anom, na.rm = TRUE),
            EBS = mean(EBS.anom, na.rm = TRUE),
            NP  = mean(NP.anom,  na.rm = TRUE),
            .groups = "drop") %>%
  rename(year = win.year) %>%
  arrange(year)

# --- Step 6: 15-year rolling AR(1) and SD (right-aligned) ---
roll_ar1 <- function(x, width = 15) {
  rollapply(x, width = width, fill = NA, align = "right",
            FUN = function(v) acf(v, lag.max = 1, plot = FALSE)$acf[2])
}
roll_sd <- function(x, width = 15) {
  rollapply(x, width = width, fill = NA, align = "right", FUN = sd)
}

winter <- winter %>%
  mutate(GOA.ar1 = roll_ar1(GOA),
         EBS.ar1 = roll_ar1(EBS),
         NP.ar1  = roll_ar1(NP),
         GOA.sd  = roll_sd(GOA),
         EBS.sd  = roll_sd(EBS),
         NP.sd   = roll_sd(NP))

# --- Step 7: plot ---
region.colors <- c("GOA" = "steelblue4", "EBS" = "darkorange4", "NP" = "forestgreen")

ar1.plot <- winter %>%
  select(year, GOA.ar1, EBS.ar1, NP.ar1) %>%
  pivot_longer(-year, names_to = "region", values_to = "ar1") %>%
  mutate(region = sub(".ar1", "", region, fixed = TRUE)) %>%
  ggplot(aes(x = year, y = ar1, color = region)) +
  geom_line() +
  geom_point(size = 1.5) +
  scale_color_manual(values = region.colors) +
  labs(title = "15-year Rolling AR(1) — Winter SST Anomaly",
       x = "Year", y = "AR(1)", color = "Region") +
  theme_bw() +
  theme(legend.position = "bottom")

sd.plot <- winter %>%
  select(year, GOA.sd, EBS.sd, NP.sd) %>%
  pivot_longer(-year, names_to = "region", values_to = "sd") %>%
  mutate(region = sub(".sd", "", region, fixed = TRUE)) %>%
  ggplot(aes(x = year, y = sd, color = region)) +
  geom_line() +
  geom_point(size = 1.5) +
  scale_color_manual(values = region.colors) +
  labs(title = "15-year Rolling SD — Winter SST Anomaly",
       x = "Year", y = "SD (°C)", color = "Region") +
  theme_bw() +
  theme(legend.position = "bottom")

print(ar1.plot / sd.plot)

# AUTOCORRELATION BY ERA ----------------------------------

# Use full monthly detrended anomaly time series (all months)
era.breaks <- list(
  "1950-1988" = c(1950, 1988),
  "1989-2000" = c(1989, 2000),
  "2001-2025" = c(2001, 2025)
)
era.colors <- c("1950-1988" = "steelblue4", "1989-2000" = "darkorange3", "2001-2025" = "firebrick3")

max.lag <- 60

# Compute ACF for one time series and return as data frame
acf_df <- function(x, max.lag) {
  a <- acf(x, lag.max = max.lag, plot = FALSE, na.action = na.pass)
  data.frame(lag = as.integer(a$lag), acf = as.numeric(a$acf))
}

# Build ACF data frame for all domains and eras
acf.out <- data.frame()

for (era in names(era.breaks)) {
  yr.range <- era.breaks[[era]]
  dat <- monthly %>%
    filter(year >= yr.range[1], year <= yr.range[2]) %>%
    arrange(year, month)

  for (dom in c("GOA", "EBS", "NP")) {
    col <- paste0(dom, ".anom")
    a   <- acf_df(dat[[col]], max.lag)
    acf.out <- rbind(acf.out,
                     data.frame(era = era, domain = dom, lag = a$lag, acf = a$acf))
  }
}

# Plot — one panel per domain
acf.plot <- ggplot(acf.out, aes(x = lag, y = acf, color = era)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "gray50") +
  geom_line() +
  geom_point(size = 1.2) +
  facet_wrap(~ domain, ncol = 1) +
  scale_color_manual(values = era.colors) +
  scale_x_continuous(breaks = seq(0, max.lag, by = 12)) +
  labs(title = "Autocorrelation of Monthly SST Anomaly by Era",
       x = "Lag (months)", y = "Autocorrelation", color = "Era") +
  theme_bw() +
  theme(legend.position = "bottom")

print(acf.plot)

# MONTHLY SST ANOMALY TIME SERIES ----------------------------------

ts.plot <- monthly.anom %>%
  pivot_longer(cols = c(GOA.anom, EBS.anom, NP.anom),
               names_to = "region", values_to = "anom") %>%
  mutate(region = recode(region,
                         "GOA.anom" = "GOA",
                         "EBS.anom" = "EBS",
                         "NP.anom"  = "NP")) %>%
  ggplot(aes(x = date, y = anom, color = region)) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = c("GOA" = "steelblue4", "EBS" = "darkorange4", "NP" = "forestgreen")) +
  facet_wrap(~ region, ncol = 1, scale = "free_y") +
  labs(title = "Monthly SST Anomaly (1950–2025)",
       x = "Year", y = "SST Anomaly (°C)", color = "Region") +
  theme_bw() +
  theme(legend.position = "none")

print(ts.plot)

# ANNUAL WINTER SST ANOMALY TIME SERIES (non-detrended) ----------------------------------

winter.anom <- monthly.anom %>%
  filter(month %in% c(11, 12, 1, 2, 3)) %>%
  mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
  group_by(win.year) %>%
  summarise(GOA = mean(GOA.anom, na.rm = TRUE),
            EBS = mean(EBS.anom, na.rm = TRUE),
            NP  = mean(NP.anom,  na.rm = TRUE),
            .groups = "drop") %>%
  rename(year = win.year)

winter.anom.plot <- winter.anom %>%
  pivot_longer(-year, names_to = "region", values_to = "anom") %>%
  ggplot(aes(x = year, y = anom, color = region)) +
  geom_line() +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("GOA" = "steelblue4", "EBS" = "darkorange4", "NP" = "forestgreen")) +
  facet_wrap(~ region, ncol = 1, scales = "free_y") +
  labs(title = "Annual Winter SST Anomaly (Nov–Mar, 1950–2025)",
       x = "Year", y = "SST Anomaly (°C)") +
  theme_bw() +
  theme(legend.position = "none")

print(winter.anom.plot)

# AL SLP SD vs. SST AR1 ----------------------------------

# Load Aleutian Low winter SLP anomaly
al <- read.csv("./Output/AL_winter_SLP_anomaly.csv")

# 15-year rolling SD of AL SLP (right-aligned), z-scored to center on 0
al <- al %>%
  arrange(year) %>%
  mutate(AL.sd = roll_sd(SLP),
         AL.sd = as.numeric(scale(AL.sd)))

# Join AL SD with SST rolling AR1 from winter data frame
al.sst <- winter %>%
  select(year, GOA.ar1, EBS.ar1) %>%
  left_join(al %>% select(year, AL.sd), by = "year") %>%
  filter(!is.na(AL.sd))

# Fit GLS with AR(1) residuals; extract slope p-value and marginal R²
fit_gls_stats <- function(dat, y.col) {
  dat <- dat %>% filter(!is.na(.data[[y.col]]), !is.na(AL.sd))
  fit <- gls(as.formula(paste(y.col, "~ AL.sd")), data = dat,
             correlation = corAR1(form = ~ 1), method = "ML")
  p   <- summary(fit)$tTable["AL.sd", "p-value"]
  # marginal R²: cor² between fitted and observed
  r2  <- cor(dat[[y.col]], fitted(fit))^2
  list(p = p, r2 = r2)
}

goa.stats <- fit_gls_stats(al.sst, "GOA.ar1")
ebs.stats <- fit_gls_stats(al.sst, "EBS.ar1")

fmt_label <- function(s) {
  sprintf("R² = %.2f,  p = %.3f", s$r2, s$p)
}

# Two-panel scatter plot
goa.scatter <- ggplot(al.sst, aes(x = AL.sd, y = GOA.ar1)) +
  geom_point(color = "steelblue4", size = 2) +
  geom_smooth(method = "lm", color = "steelblue4", fill = "steelblue", alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
           label = fmt_label(goa.stats), size = 3.5) +
  labs(title = "GOA",
       x = "AL SLP SD (z-scored, 15-yr window)", y = "SST AR(1)") +
  theme_bw()

ebs.scatter <- ggplot(al.sst, aes(x = AL.sd, y = EBS.ar1)) +
  geom_point(color = "darkorange4", size = 2) +
  geom_smooth(method = "lm", color = "darkorange4", fill = "darkorange", alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
           label = fmt_label(ebs.stats), size = 3.5) +
  labs(title = "EBS",
       x = "AL SLP SD (z-scored, 15-yr window)", y = "SST AR(1)") +
  theme_bw()

print(goa.scatter / ebs.scatter)

# TIME SERIES: SST AR1 vs. AL SLP SD — multiple window lengths ----------------------------------

windows <- c(10, 15, 20, 25)

# Build long-format data frame with AR1 and AL SD for each window
multi.win <- lapply(windows, function(w) {

  # SST AR1 for this window
  sst.w <- winter %>%
    arrange(year) %>%
    mutate(GOA.ar1 = rollapply(GOA, width = w, fill = NA, align = "right",
                               FUN = function(v) acf(v, lag.max = 1, plot = FALSE)$acf[2]),
           EBS.ar1 = rollapply(EBS, width = w, fill = NA, align = "right",
                               FUN = function(v) acf(v, lag.max = 1, plot = FALSE)$acf[2]))

  # AL SLP SD for this window, z-scored to center on 0
  al.w <- al %>%
    arrange(year) %>%
    mutate(AL.sd = rollapply(SLP, width = w, fill = NA, align = "right", FUN = sd),
           AL.sd = as.numeric(scale(AL.sd)))

  sst.w %>%
    select(year, GOA.ar1, EBS.ar1) %>%
    left_join(al.w %>% select(year, AL.sd), by = "year") %>%
    mutate(window = paste0(w, "-yr"))

}) %>% bind_rows()

multi.win$window <- factor(multi.win$window,
                           levels = paste0(windows, "-yr"))

# Long format: one row per year/window/variable
ts.long <- multi.win %>%
  pivot_longer(cols = c(GOA.ar1, EBS.ar1, AL.sd),
               names_to = "variable", values_to = "value") %>%
  mutate(
    region = case_when(
      variable == "GOA.ar1" ~ "GOA",
      variable == "EBS.ar1" ~ "EBS",
      variable == "AL.sd"   ~ "GOA"   # AL SD plotted in both rows via duplication below
    ),
    series = ifelse(variable == "AL.sd", "AL SLP SD", "SST AR(1)")
  )

# Duplicate AL SD rows for EBS panel
al.ebs <- ts.long %>%
  filter(variable == "AL.sd") %>%
  mutate(region = "EBS")

ts.long <- bind_rows(ts.long %>% filter(variable != "AL.sd" | region == "GOA"),
                     al.ebs)

ts.multi.plot <- ggplot(ts.long %>% filter(!is.na(value)),
                        aes(x = year, y = value,
                            color = series, linetype = series, shape = series)) +
  geom_line() +
  geom_point(size = 1.2) +
  scale_color_manual(values   = c("SST AR(1)" = "steelblue4", "AL SLP SD" = "gray30")) +
  scale_linetype_manual(values = c("SST AR(1)" = "solid",     "AL SLP SD" = "dashed")) +
  scale_shape_manual(values    = c("SST AR(1)" = 16,          "AL SLP SD" = 1)) +
  facet_grid(region ~ window, scales = "free_y") +
  labs(title = "SST AR(1) and AL SLP SD — multiple rolling windows",
       x = "Year", y = NULL,
       color = NULL, linetype = NULL, shape = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

print(ts.multi.plot)

# SCATTER PLOTS — all window x region combinations ----------------------------------

# Build long-format scatter data from multi.win
scatter.long <- multi.win %>%
  pivot_longer(cols = c(GOA.ar1, EBS.ar1),
               names_to = "region", values_to = "ar1") %>%
  mutate(region = sub(".ar1", "", region, fixed = TRUE)) %>%
  filter(!is.na(ar1), !is.na(AL.sd))

# GLS stats for each window x region combination
gls.stats <- scatter.long %>%
  group_by(window, region) %>%
  group_modify(~ {
    fit <- tryCatch(
      gls(ar1 ~ AL.sd, data = .x,
          correlation = corAR1(form = ~ 1), method = "ML"),
      error = function(e) NULL
    )
    if (is.null(fit)) return(data.frame(label = "p = NA,  R² = NA"))
    p  <- summary(fit)$tTable["AL.sd", "p-value"]
    r2 <- cor(.x$ar1, fitted(fit))^2
    data.frame(label = sprintf("R² = %.2f,  p = %.3f", r2, p))
  }) %>%
  ungroup()

region.colors <- c("GOA" = "steelblue4", "EBS" = "darkorange4")

scatter.multi.plot <- ggplot(scatter.long, aes(x = AL.sd, y = ar1, color = region)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", aes(fill = region), alpha = 0.2) +
  geom_text(data = gls.stats, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
            inherit.aes = FALSE, size = 3) +
  scale_color_manual(values = region.colors) +
  scale_fill_manual(values  = region.colors) +
  facet_grid(region ~ window, scales = "free_y") +
  labs(title = "AL SLP SD vs. SST AR(1) — multiple rolling windows",
       x = "AL SLP SD (z-scored)", y = "SST AR(1)") +
  theme_bw() +
  theme(legend.position = "none")

print(scatter.multi.plot)


