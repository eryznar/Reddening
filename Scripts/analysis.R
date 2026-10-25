# PURPOSE: Reddening project — main analysis script
# Author: Mike Litzow
#
# Sections:
#   1. Reference map — North Pacific study areas

library(ggplot2)
library(dplyr)
library(maps)
library(ncdf4)
library(tidyr)
library(patchwork)

# ============================================================
# SHARED: MANUSCRIPT FIGURE STYLE (GRL) + MAP PROJECTION HELPERS
# ============================================================
# Style targets Geophysical Research Letters:
#   - Column widths: 95 mm (single) / 190 mm (full). Use grl.width.*.
#   - Sans-serif (Helvetica/Arial), 8 pt base, 6 pt hard minimum.
#   - No plot titles — GRL figures carry captions, not in-figure titles.
#   - Hairline strokes (<= 0.4 pt) and no heavy panel gridlines.
# Projection helpers implement a spherical north polar stereographic
# transform directly (no sf/mapproj dependency) so that map panels can be
# drawn in a projection appropriate for high-latitude North Pacific fields.

grl.width.single <- 95  / 25.4   # inches
grl.width.full   <- 190 / 25.4   # inches

theme_grl <- function(base_size = 8, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      text             = element_text(colour = "black"),
      axis.text        = element_text(size = base_size - 1, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.25, colour = "black"),
      legend.title     = element_text(size = base_size),
      legend.text      = element_text(size = base_size - 1),
      legend.key.width = unit(3.5, "mm"),
      legend.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.2),
      panel.border     = element_rect(colour = "black", linewidth = 0.4, fill = NA),
      plot.title       = element_text(size = base_size + 1, hjust = 0),
      plot.subtitle    = element_text(size = base_size - 1, hjust = 0),
      plot.margin      = margin(2, 2, 2, 2, "mm"),
      strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.3),
      strip.text       = element_text(size = base_size))
}

# Map variant: graticule and labels are drawn as data, so panel axes,
# gridlines and border are all switched off.
theme_grl_map <- function(base_size = 8, base_family = "sans") {
  theme_grl(base_size, base_family) +
    theme(axis.text         = element_blank(),
          axis.ticks        = element_blank(),
          axis.title        = element_blank(),
          panel.grid.major  = element_blank(),
          panel.grid.minor  = element_blank(),
          panel.border      = element_blank())
}

# --- North polar stereographic (spherical) --------------------------------
# lon/lat in degrees (lon in 0-360 space); lon0 = central meridian, drawn
# pointing straight down from the pole. Returns x/y in Earth radii, pole at
# the origin. East of lon0 falls to the right, matching a view from above
# the North Pole.
ps_project <- function(lon, lat, lon0 = 205) {
  lam <- (lon - lon0) * pi / 180
  r   <- 2 * tan(pi / 4 - (lat * pi / 180) / 2)
  data.frame(x = r * sin(lam), y = -r * cos(lam))
}

# Insert intermediate vertices so that straight lon/lat edges render as the
# curves they actually are once projected.
densify_lonlat <- function(lon, lat, max.step = 0.5) {
  n <- length(lon)
  if (n < 2) return(list(lon = lon, lat = lat))
  k    <- pmax(1, ceiling(pmax(abs(diff(lon)), abs(diff(lat))) / max.step))
  idx  <- rep(seq_len(n - 1), k)
  frac <- unlist(lapply(k, function(m) seq(0, 1, length.out = m + 1)[-(m + 1)]))
  list(lon = c(lon[idx] + frac * (lon[idx + 1] - lon[idx]), lon[n]),
       lat = c(lat[idx] + frac * (lat[idx + 1] - lat[idx]), lat[n]))
}

# Sutherland-Hodgman clip of a polygon against a lon/lat rectangle. Keeps
# coastlines from spilling outside the mapped domain once projected.
clip_lonlat <- function(lon, lat, lon.rng, lat.rng) {
  clip1 <- function(px, py, keep, cut) {
    n <- length(px)
    if (n < 3) return(list(x = numeric(0), y = numeric(0)))
    j  <- c(n, seq_len(n - 1))
    ci <- keep(px, py); cj <- ci[j]
    it <- cut(px[j], py[j], px, py)
    xs <- rbind(it$x, px); ys <- rbind(it$y, py)
    sel <- rbind(ci != cj, ci)
    list(x = xs[sel], y = ys[sel])
  }
  cut_x <- function(xc) function(x1, y1, x2, y2) {
    t <- (xc - x1) / (x2 - x1); list(x = rep(xc, length(t)), y = y1 + t * (y2 - y1))
  }
  cut_y <- function(yc) function(x1, y1, x2, y2) {
    t <- (yc - y1) / (y2 - y1); list(x = x1 + t * (x2 - x1), y = rep(yc, length(t)))
  }
  p <- list(x = lon, y = lat)
  for (e in list(
    list(k = function(x, y) x >= lon.rng[1], c = cut_x(lon.rng[1])),
    list(k = function(x, y) x <= lon.rng[2], c = cut_x(lon.rng[2])),
    list(k = function(x, y) y >= lat.rng[1], c = cut_y(lat.rng[1])),
    list(k = function(x, y) y <= lat.rng[2], c = cut_y(lat.rng[2])))) {
    p <- clip1(p$x, p$y, e$k, e$c)
    if (!length(p$x)) break
  }
  p
}

# Coastline polygons clipped to a lon/lat window, densified and projected.
# The ring is closed (first vertex repeated) BEFORE densifying: geom_polygon
# closes a ring with a straight segment in projected space, and along the
# northern clip edge that chord bows outside the curved boundary, spilling
# land past the map outline. Densifying the closing edge removes it.
ps_coast <- function(map.df, lon.rng, lat.rng, lon0 = 205) {
  parts <- lapply(split(map.df, map.df$group), function(g) {
    cp <- clip_lonlat(g$long, g$lat, lon.rng, lat.rng)
    if (length(cp$x) < 3) return(NULL)
    d <- densify_lonlat(c(cp$x, cp$x[1]), c(cp$y, cp$y[1]))
    cbind(ps_project(d$lon, d$lat, lon0), group = g$group[1])
  })
  bind_rows(parts[!vapply(parts, is.null, logical(1))])
}

# Densify a lon/lat path, then project it.
ps_path <- function(lon, lat, lon0 = 205, n = 200) {
  i  <- seq_along(lon)
  ii <- seq(1, length(lon), length.out = n)
  ps_project(approx(i, lon, ii)$y, approx(i, lat, ii)$y, lon0)
}

# Expand cell centres into projected quadrilaterals. `half` is half the cell
# size in degrees; `value` is carried through for the fill aesthetic. Cells on
# the domain edge are centred on the boundary, so their corners are clamped to
# lon.rng/lat.rng — without this the field overhangs the map outline by half a
# cell. Quads are axis-aligned in lon/lat, so clamping corners clips exactly.
ps_cells <- function(lon, lat, value, half = 0.25, lon0 = 205,
                     lon.rng = NULL, lat.rng = NULL) {
  n     <- length(lon)
  dlon  <- c(-half,  half, half, -half)
  dlat  <- c(-half, -half, half,  half)
  lonc  <- rep(lon, each = 4) + rep(dlon, n)
  latc  <- rep(lat, each = 4) + rep(dlat, n)
  if (!is.null(lon.rng)) lonc <- pmin(pmax(lonc, lon.rng[1]), lon.rng[2])
  if (!is.null(lat.rng)) latc <- pmin(pmax(latc, lat.rng[1]), lat.rng[2])
  xy <- ps_project(lonc, latc, lon0)
  data.frame(x = xy$x, y = xy$y,
             id = rep(seq_len(n), each = 4),
             value = rep(value, each = 4))
}

# Text with a white halo, so labels stay legible over a filled field without
# an opaque label box hiding data. Returns a list of layers: `n` offset copies
# in the halo colour, then the text itself. `r` is the halo radius in data
# units (projected Earth radii). A `hjust` column in `data` aligns each label;
# `...` passes further fixed parameters (vjust, lineheight) to geom_text().
geom_text_halo <- function(data, size = 7 / .pt, colour = "black",
                           halo = "white", r = 0.0028, n = 12,
                           fontface = "bold", ...) {
  map <- if ("hjust" %in% names(data))
    aes(x = x, y = y, label = lab, hjust = hjust) else aes(x = x, y = y, label = lab)
  ang <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
  bg  <- lapply(ang, function(a) {
    d <- data; d$x <- d$x + r * cos(a); d$y <- d$y + r * sin(a)
    geom_text(data = d, map, colour = halo, size = size, fontface = fontface, ...)
  })
  c(bg, list(geom_text(data = data, map, colour = colour, size = size,
                       fontface = fontface, ...)))
}

# Clean 0-360 world map used by every projected map panel.
mapWorld.clean <- map_data("world2")

# Complete manuscript map panel: a cellwise regression field drawn on the
# north polar stereographic projection, with the FDR contour, clipped
# coastlines, graticule, EBS/GOA polygons and the AL box. Every regression
# map in the manuscript goes through this so the panels are guaranteed to
# share a projection, domain and layout. `grid` is a regular lon/lat data
# frame carrying `beta` and optionally `qval`. Study polygons and the AL box
# are read from the global environment, as elsewhere in the script.
# `land.fill = NA` leaves land unfilled so an atmospheric field shows
# through it. `fill.scale` replaces the default diverging scale.
# `grat.size` and `line.scale` enlarge the lat/lon labels and every stroke
# for maps drawn bigger than a journal panel (the slides of Section 27).
ps_regression_map <- function(grid, fill.name, col.lim,
                              lon.rng = c(160, 250), lat.rng = c(20, 66),
                              lon0 = 205, cell.half = 0.25, q.level = 0.05,
                              sig.colour = "#1B7837",
                              label.size = 7 / .pt,     # EBS/GOA/AL Index labels
                              grat.size  = 7 / .pt,     # lat/lon labels
                              line.scale = 1,           # multiplies every stroke
                              land.fill = "grey85",
                              fill.scale = NULL) {
  # Cells are drawn 20% oversize so neighbours overlap: exactly abutting
  # polygons leave hairline anti-aliasing seams in both PNG and PDF output.
  # Colour boundaries move by a tenth of a cell, which is not visible.
  cells <- grid %>% filter(!is.na(beta)) %>%
    { ps_cells(.$lon, .$lat, .$beta, half = cell.half * 1.2, lon0 = lon0,
               lon.rng = lon.rng, lat.rng = lat.rng) }

  # Significance is contoured, in a colour that separates it from the black
  # study polygons and the dashed AL box — stippling at this figure scale is
  # too fine to register. The contour is traced on the lon/lat grid and only
  # then projected, so reprojection cannot alter the statistics shown.
  qpath <- NULL
  if ("qval" %in% names(grid) && any(!is.na(grid$qval))) {
    q.lon <- sort(unique(grid$lon)); q.lat <- sort(unique(grid$lat))
    q.mat <- matrix(NA_real_, nrow = length(q.lon), ncol = length(q.lat))
    q.mat[cbind(match(grid$lon, q.lon), match(grid$lat, q.lat))] <- grid$qval
    q.lines <- contourLines(x = q.lon, y = q.lat, z = q.mat, levels = q.level)
    if (length(q.lines))
      qpath <- bind_rows(lapply(seq_along(q.lines), function(i)
        cbind(ps_project(q.lines[[i]]$x, q.lines[[i]]$y, lon0), seg = i)))
  }

  coast <- ps_coast(mapWorld.clean, lon.rng, lat.rng, lon0)
  areas <- bind_rows(
    cbind(ps_path(c(goa.poly$lon, goa.poly$lon[1]),
                  c(goa.poly$lat, goa.poly$lat[1]), lon0, n = 500), area = "GOA"),
    cbind(ps_path(c(ebs.poly$lon, ebs.poly$lon[1]),
                  c(ebs.poly$lat, ebs.poly$lat[1]), lon0, n = 500), area = "EBS"))
  # Anchored toward each polygon's top-right, at the most top-right position
  # whose whole text extent still clears the mainland coast. Extent is wider
  # in degrees than it looks: meridians converge poleward, so "EBS" spans
  # about +/-4 deg of longitude at 63N. EBS cannot sit at its true top-right
  # corner, which is mainland Alaska. "AL Index" spans ~10 deg at 47N, so its
  # anchor is the leftmost point keeping the whole label inside the AL box.
  area.lab <- bind_rows(
    cbind(ps_project(217.0, 58.5, lon0), lab = "GOA"),
    cbind(ps_project(190.5, 63.5, lon0), lab = "EBS"),
    cbind(ps_project(198.0, 46.8, lon0), lab = "AL Index"))
  al.box <- ps_path(al.box.df$x, al.box.df$y, lon0, n = 400)
  domain <- ps_path(lon.rng[c(1, 2, 2, 1, 1)], lat.rng[c(1, 1, 2, 2, 1)],
                    lon0, n = 800)
  grat   <- ps_graticule(lon.rng, lat.rng, lon.by = 20, lat.by = 10, lon0 = lon0)

  xlim <- range(domain$x) + c(-0.13, 0.03)   # room for lat labels
  ylim <- range(domain$y) + c(-0.10, 0.02)   # room for lon labels

  ggplot() +
    geom_polygon(data = cells, aes(x = x, y = y, group = id, fill = value)) +
    { if (!is.null(qpath))
        geom_path(data = qpath, aes(x = x, y = y, group = seg),
                  colour = sig.colour, linewidth = 0.45 * line.scale) } +
    geom_polygon(data = coast, aes(x = x, y = y, group = group),
                 fill = land.fill, colour = "grey30", linewidth = 0.2 * line.scale) +
    geom_path(data = grat$lines, aes(x = x, y = y, group = grat),
              colour = "grey45", linewidth = 0.15 * line.scale, linetype = "dotted") +
    geom_path(data = al.box, aes(x = x, y = y),
              colour = "black", linewidth = 0.25 * line.scale, linetype = "dashed") +
    geom_path(data = areas, aes(x = x, y = y, group = area),
              colour = "black", linewidth = 0.4 * line.scale) +
    geom_path(data = domain, aes(x = x, y = y),
              colour = "black", linewidth = 0.4 * line.scale) +
    geom_text(data = grat$labels, aes(x = x, y = y, label = lab),
              size = grat.size, colour = "black") +
    geom_text_halo(area.lab, size = label.size) +
    { if (is.null(fill.scale))
        scale_fill_gradient2(low = "steelblue4", mid = "white",
                             high = "darkorange4", midpoint = 0,
                             limits = c(-col.lim, col.lim), name = fill.name,
                             oob = scales::squish)
      else fill.scale } +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    guides(fill = guide_colourbar(barwidth  = unit(55, "mm"),
                                  barheight = unit(3, "mm"),
                                  title.position = "left",
                                  title.vjust = 1,
                                  frame.colour = "black",
                                  frame.linewidth = 0.3,
                                  ticks.colour = "black")) +
    theme_grl_map() +
    theme(legend.position  = "bottom",
          legend.direction = "horizontal",
          legend.margin    = margin(1, 0, 0, 0, "mm"),
          legend.box.spacing = unit(1, "mm"))
}

# Graticule as projected paths. Labels are pushed a fixed distance outside
# the domain edge in projected units — a constant offset in degrees would
# crowd the high-latitude labels, where the meridians converge.
ps_graticule <- function(lon.rng, lat.rng, lon.by = 20, lat.by = 10,
                         lon0 = 205, off = 0.055) {
  mers <- seq(ceiling(lon.rng[1] / lon.by) * lon.by, lon.rng[2], by = lon.by)
  pars <- seq(ceiling(lat.rng[1] / lat.by) * lat.by, lat.rng[2], by = lat.by)
  lines <- bind_rows(
    bind_rows(lapply(mers, function(L)
      cbind(ps_path(c(L, L), lat.rng, lon0), grat = paste0("m", L)))),
    bind_rows(lapply(pars, function(P)
      cbind(ps_path(lon.rng, c(P, P), lon0), grat = paste0("p", P)))))

  # Unit vectors at each anchor: radial (outward, away from the pole) and
  # tangential (eastward). Meridian labels step radially outward past the
  # southern arc; parallel labels step west of the domain's western edge.
  lam.m <- (mers - lon0) * pi / 180
  lam.p <- (lon.rng[1] - lon0) * pi / 180
  m.xy  <- ps_project(mers, lat.rng[1], lon0)
  p.xy  <- ps_project(lon.rng[1], pars, lon0)

  labs <- bind_rows(
    data.frame(x = m.xy$x + off * sin(lam.m),
               y = m.xy$y - off * cos(lam.m),
               lab = lon_label(mers), edge = "lon"),
    data.frame(x = p.xy$x - off * cos(lam.p),
               y = p.xy$y - off * sin(lam.p),
               lab = paste0(pars, "°N"), edge = "lat"))
  list(lines = lines, labels = labs)
}

# Normalise longitude to 0-360, whichever convention a file supplies.
to360 <- function(x) ifelse(x < 0, x + 360, x)

# Generic right-aligned rolling window: the base-R equivalent of
# zoo::rollapply(x, width, fill = NA, align = "right", FUN = FUN).
roll_right <- function(x, width, FUN, ...) {
  n <- length(x); out <- rep(NA_real_, n)
  if (n < width) return(out)
  for (i in width:n) out[i] <- FUN(x[(i - width + 1):i], ...)
  out
}

# Leading `nv` singular triplets of a tall matrix, replacing irlba::irlba.
# Every use here has far fewer columns (times) than rows (cells), so the
# Gram matrix t(A) %*% A is small and its eigen-decomposition gives the
# thin SVD directly. Returns d/u/v like irlba, so sign conventions are
# still fixed downstream where the script flips EOFs explicitly.
partial_svd <- function(A, nv) {
  e <- eigen(crossprod(A), symmetric = TRUE)
  k <- seq_len(nv)
  d <- sqrt(pmax(e$values[k], 0))
  v <- e$vectors[, k, drop = FALSE]
  list(d = d, u = A %*% v %*% diag(1 / d, nv, nv), v = v)
}

# coord_map(projection = "rectangular", parameters = lat0) without mapproj:
# that projection is just x = lon * cos(lat0), y = lat, i.e. a fixed aspect.
coord_rect <- function(xlim, ylim, lat0 = 55) {
  coord_fixed(ratio = 1 / cos(lat0 * pi / 180), xlim = xlim, ylim = ylim)
}

# --- Small numeric helpers (base R, so no zoo/sf dependency) --------------
# Vectorised ray-casting point-in-polygon: mask many grid cells at once.
inpoly_vec <- function(x, y, px, py) {
  n <- length(px); inside <- rep(FALSE, length(x)); j <- n
  for (i in seq_len(n)) {
    cond <- ((py[i] > y) != (py[j] > y)) &
            (x < (px[j] - px[i]) * (y - py[i]) / (py[j] - py[i]) + px[i])
    inside <- xor(inside, cond %in% TRUE)
    j <- i
  }
  inside
}

# Linear detrend, NA-safe.
detrend_lin <- function(x) {
  i <- seq_along(x); ok <- !is.na(x); out <- rep(NA_real_, length(x))
  out[ok] <- residuals(lm(x[ok] ~ i[ok]))
  out
}

# Right-aligned rolling lag-1 autocorrelation, matching Section 3's
# rollapply(..., align = "right") without requiring zoo.
roll_ar1_right <- function(x, width = 15) {
  roll_right(x, width,
             function(v) acf(v, lag.max = 1, plot = FALSE,
                             na.action = na.pass)$acf[2])
}

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
  coord_rect(xlim = c(120, 250), ylim = c(20, 66)) +
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

# Grid and polygons are both held in 0-360 space. A previous version
# converted the polygons to -180/180 while the grid kept the file's own
# longitudes -- which the 2026 ERA5 file supplies as 110-250 -- so polygons
# and grid never overlapped and both masks came back empty. to360() now
# normalises the grid, so either file convention works.
grid <- expand.grid(lon = to360(as.vector(lons)), lat = as.vector(lats))

in.goa <- inpoly_vec(grid$lon, grid$lat, goa.x, goa.y)
in.ebs <- inpoly_vec(grid$lon, grid$lat, ebs.x, ebs.y)
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
  roll_right(x, width, FUN = function(v) acf(v, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2])
}
roll_sd <- function(x, width = 15) {
  roll_right(x, width, FUN = sd)
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
  roll_right(x, width, FUN = mean)
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
sv.slp   <- partial_svd(sweep(detrend.slp, 1, w.slp, "*"), nv = 1)
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
  coord_rect(xlim = c(120, 250), ylim = c(20, 66)) +
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
  sv.e <- partial_svd(Aw, nv = 5)

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
    coord_rect(xlim = mld.eof.domain$lon, ylim = mld.eof.domain$lat) +
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

# Regular 0.5-deg grid for the projected field + FDR contouring (ORAS5 is
# irregular). Restricted to the Section 17 map domain so the two manuscript
# maps share identical lat/lon limits. Note mld.reg.df itself keeps its full
# 110-250E extent: Section 6C still plots that field from 150E.
# FDR is recomputed over the plotted cells alone: correcting across the full
# 110-250E field would penalise this map for tests it never shows. (Section
# 6C still uses the wider-domain qval carried on mld.reg.df.)
mld.reg.plot <- mld.reg.df %>% filter(lon >= 160, lon <= 250)
mld.reg.plot$qval <- p.adjust(mld.reg.plot$pval, method = "BH")
message(sprintf("Plotted domain only: BH-FDR q<=0.05: %d of %d cells",
                sum(mld.reg.plot$qval <= 0.05, na.rm = TRUE), nrow(mld.reg.plot)))

mld.reg.grid <- mld.reg.plot %>%
  mutate(lon = round(lon / 0.5) * 0.5,
         lat = round(lat / 0.5) * 0.5) %>%
  group_by(lon, lat) %>%
  summarise(beta = mean(beta, na.rm = TRUE),
            qval = mean(qval, na.rm = TRUE),
            .groups = "drop") %>%
  complete(lon = seq(160, 250, by = 0.5),
           lat = seq(20,  66,  by = 0.5))

# Colour limit clipped at the 98th percentile of |beta| (values beyond are
# squished to the end colours). max|beta| is set by a thin tail and leaves
# the bulk of the field washed out.
col.lim.mld <- as.numeric(quantile(abs(mld.reg.plot$beta), 0.98, na.rm = TRUE))

p.mld.reg <- ps_regression_map(
  mld.reg.grid,
  fill.name = expression(beta~"(m / hPa)"),
  col.lim   = col.lim.mld)

ggsave("./Figures/MLD_AL_SLP_regression.png",
       plot = p.mld.reg, width = grl.width.full, height = 4.8,
       dpi = 400, bg = "white")
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
  coord_rect(xlim = c(110, 250), ylim = c(20, 66)) +
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
    coord_rect(xlim = c(150, 250), ylim = c(20, 66)) +
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
# MLD box NetCDFs, masks to the GOA/EBS polygons (goa.poly / ebs.poly
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

  region_ts_mld <- function(d, poly.lon, poly.lat) {
    lons.v <- as.vector(d$lons); lats.v <- as.vector(d$lats)
    in.poly <- inpoly_vec(to360(lons.v), lats.v, poly.lon, poly.lat)
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

  goa.mld.ts <- region_ts_mld(goa.d, goa.poly$lon, goa.poly$lat)
  ebs.mld.ts <- region_ts_mld(ebs.d, ebs.poly$lon, ebs.poly$lat)

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
  roll_right(x, width, FUN = mean, na.rm = TRUE)
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
# SECTION 8B: CHECKS — AL & MLD ROLLING SD vs. AR(1)
# ============================================================
# Rules out an alternate reading of the SD-driven results: that the AL
# SD signal is really the AL forcing itself reddening (SD confounded with
# AR(1)), and likewise that MLD SD tracks MLD AR(1). Each panel shows the
# 15-yr right-aligned rolling SD (dashed/open, left axis) and AR(1) (solid,
# right axis) on dual y-axes. Produces two figures:
#   (1) AL_SD_vs_AR1_check.png  — standalone Aleutian Low panel (al$SLP
#       from Section 5; SD in hPa, not the z-scored series).
#   (2) MLD_SD_vs_AR1_check.png — two panels, GOA (top) and EBS (bottom),
#       from mld.win (Section 7).
# The per-series Pearson r(SD, AR(1)) is annotated: |r| small => SD and
# AR(1) carry distinct information (alternate signal ruled out).

# --- Rolling series (reuse roll_ar1 / roll_sd from Section 3) ---
al.chk <- al %>% arrange(year) %>%
  mutate(AL.sd = roll_sd(SLP), AL.ar1 = roll_ar1(SLP))

mld.chk <- mld.win %>% arrange(year) %>%
  mutate(GOA.sd = roll_sd(GOA), GOA.ar1 = roll_ar1(GOA),
         EBS.sd = roll_sd(EBS), EBS.ar1 = roll_ar1(EBS))

# --- Dual-axis SD-vs-AR(1) panel builder (single series) ---
sd_ar1_panel <- function(dat, sd.col, ar1.col, ttl, sd.name, ar1.name, ar1.color) {
  sd.rng  <- range(dat[[sd.col]],  na.rm = TRUE)
  ar1.rng <- range(dat[[ar1.col]], na.rm = TRUE)
  k <- diff(sd.rng) / diff(ar1.rng)
  b <- sd.rng[1] - k * ar1.rng[1]
  rlab <- sprintf("r(SD, AR1) = %.2f",
                  cor(dat[[sd.col]], dat[[ar1.col]], use = "complete.obs"))
  ggplot(dat, aes(x = year)) +
    geom_line(aes(y = .data[[sd.col]], color = "SD"), linetype = "dashed", linewidth = 0.7) +
    geom_point(aes(y = .data[[sd.col]], color = "SD"), shape = 1, size = 1.5) +
    geom_line(aes(y = k * .data[[ar1.col]] + b, color = "AR(1)"), linewidth = 0.7) +
    geom_point(aes(y = k * .data[[ar1.col]] + b, color = "AR(1)"), size = 1.5) +
    annotate("text", x = -Inf, y = Inf, label = rlab, hjust = -0.1, vjust = 1.4, size = 3.5) +
    scale_y_continuous(name = sd.name, sec.axis = sec_axis(~ (. - b) / k, name = ar1.name)) +
    scale_color_manual(values = c("SD" = "gray30", "AR(1)" = ar1.color)) +
    labs(title = ttl, x = "Year", color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
}

# --- Figure 1: Aleutian Low (standalone) ---
p.al.chk <- sd_ar1_panel(al.chk, "AL.sd", "AL.ar1",
                         "Aleutian Low: 15-yr rolling SD vs AR(1)",
                         "AL SLP SD (hPa)", "AL SLP AR(1)", "firebrick")
ggsave("./Figures/AL_SD_vs_AR1_check.png",
       plot = p.al.chk, width = 8, height = 4.5, dpi = 300)
message("Saved: Figures/AL_SD_vs_AR1_check.png")

# --- Figure 2: MLD, one panel per ecosystem (GOA top, EBS bottom) ---
p.goa.chk <- sd_ar1_panel(mld.chk, "GOA.sd", "GOA.ar1",
                          "Gulf of Alaska MLD", "MLD SD (m)", "MLD AR(1)", "darkorange4")
p.ebs.chk <- sd_ar1_panel(mld.chk, "EBS.sd", "EBS.ar1",
                          "Eastern Bering Sea MLD", "MLD SD (m)", "MLD AR(1)", "steelblue4")

p.mld.checks <- (p.goa.chk / p.ebs.chk) +
  plot_annotation(
    title = "MLD checks: 15-yr rolling SD vs AR(1) by ecosystem",
    theme = theme(plot.title = element_text(size = 12, hjust = 0.5)))

ggsave("./Figures/MLD_SD_vs_AR1_check.png",
       plot = p.mld.checks, width = 8, height = 9, dpi = 300)
message("Saved: Figures/MLD_SD_vs_AR1_check.png")

# ============================================================
# SECTION 8C: CHECKS — AL SD vs MLD SD; SST AR(1) ~ MLD SD + MLD AR(1)
# ============================================================
# Check 1: 15-yr rolling AL SLP SD vs MLD SD, one panel per ecosystem,
#   same dual-axis format as Section 7 (AL SD vs MLD AR(1)). Uses al$AL.sd
#   (Section 5) and the rolling MLD SD from mld.chk (Section 8B).
# Check 2: multiple regression of SST AR(1) (winter.roll, Section 3) on
#   MLD SD and MLD AR(1) for GOA and EBS, GLS with AR(1) residuals
#   (nlme::gls + corAR1). Disentangles whether SST reddening tracks MLD
#   reddening (AR(1)) or MLD variability (SD) once both are included.
#   Predictors are z-scored so coefficients are comparable. Saves the
#   coefficient table (Output) and a forest plot.

# ---- Check 1: AL SD vs MLD SD dual-axis, one panel per ecosystem ----
al.mldsd <- mld.chk %>% select(year, GOA.sd, EBS.sd) %>%
  left_join(al %>% select(year, AL.sd), by = "year") %>%
  filter(!is.na(AL.sd))

mldsd.al.range  <- range(al.mldsd$AL.sd,  na.rm = TRUE)
goa.mldsd.range <- range(al.mldsd$GOA.sd, na.rm = TRUE)
ebs.mldsd.range <- range(al.mldsd$EBS.sd, na.rm = TRUE)

dg.sd <- al.mldsd %>% filter(!is.na(GOA.sd), !is.na(AL.sd))
de.sd <- al.mldsd %>% filter(!is.na(EBS.sd), !is.na(AL.sd))
chel.goa.sd <- chelton_test(dg.sd$GOA.sd, dg.sd$AL.sd)
chel.ebs.sd <- chelton_test(de.sd$EBS.sd, de.sd$AL.sd)

make_mldsd_dual_plot <- function(dat, sd.col, sd.color, region.label, sd.range, chel) {
  k <- diff(sd.range) / diff(mldsd.al.range)
  b <- sd.range[1] - k * mldsd.al.range[1]
  dat <- dat %>% filter(!is.na(.data[[sd.col]]), !is.na(AL.sd))
  lbl <- sprintf("r = %.2f, N* = %.1f\nt = %.2f, p = %.3f",
                 chel$r.obs, chel$n.eff, chel$t.obs, chel$pval)
  ggplot(dat, aes(x = year)) +
    geom_line(aes(y = .data[[sd.col]], color = "MLD SD"), linewidth = 0.7) +
    geom_point(aes(y = .data[[sd.col]], color = "MLD SD"), size = 1.5) +
    geom_line(aes(y = k * AL.sd + b, color = "AL SLP SD"),
              linewidth = 0.7, linetype = "dashed") +
    geom_point(aes(y = k * AL.sd + b, color = "AL SLP SD"), size = 1.5, shape = 1) +
    annotate("text", x = Inf, y = -Inf, label = lbl,
             hjust = 1.1, vjust = -0.4, size = 3.5, lineheight = 1.3) +
    scale_y_continuous(name = "MLD SD (m)",
                       sec.axis = sec_axis(~ (. - b) / k, name = "AL SLP SD (z-scored)")) +
    scale_color_manual(values = c("MLD SD" = sd.color, "AL SLP SD" = "gray30")) +
    labs(title = region.label, x = "Year", color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position    = "bottom",
          axis.title.y.left  = element_text(color = sd.color),
          axis.text.y.left   = element_text(color = sd.color),
          axis.title.y.right = element_text(color = "gray30"),
          axis.text.y.right  = element_text(color = "gray30"))
}

p.goa.mldsd <- make_mldsd_dual_plot(al.mldsd, "GOA.sd", "darkorange4", "GOA",
                                    goa.mldsd.range, chel.goa.sd)
p.ebs.mldsd <- make_mldsd_dual_plot(al.mldsd, "EBS.sd", "steelblue4",  "EBS",
                                    ebs.mldsd.range, chel.ebs.sd)

p.mldsd.dual <- p.goa.mldsd / p.ebs.mldsd +
  plot_annotation(
    title = "15-yr Rolling AL SLP SD vs. MLD SD — Winter (Nov–Mar)",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 12)))

ggsave("./Figures/AL_SD_MLD_SD_15yr_dual_axis.png",
       plot = p.mldsd.dual, width = 7, height = 7, dpi = 300)
message("Saved: Figures/AL_SD_MLD_SD_15yr_dual_axis.png")

# ---- Check 2: SST AR(1) ~ MLD SD + MLD AR(1) (GLS, AR(1) residuals) ----
chk2 <- winter.roll %>%
  select(year, SST.GOA = GOA.ar1, SST.EBS = EBS.ar1) %>%
  inner_join(mld.chk %>% select(year, GOA.sd, GOA.ar1, EBS.sd, EBS.ar1), by = "year")

fit_mld_reg <- function(sst.col, sd.col, ar1.col, sys) {
  d <- chk2 %>%
    transmute(year,
              y       = .data[[sst.col]],
              MLD.sd  = as.numeric(scale(.data[[sd.col]])),   # z-scored
              MLD.ar1 = as.numeric(scale(.data[[ar1.col]]))) %>%
    filter(complete.cases(.))
  fit <- gls(y ~ MLD.sd + MLD.ar1, data = d,
             correlation = corAR1(form = ~ 1), method = "ML")
  s <- summary(fit)$tTable
  data.frame(system = sys, term = rownames(s),
             beta = s[, 1], se = s[, 2], t = s[, 3], p = s[, 4],
             n = nrow(d), row.names = NULL)
}

reg.mld <- bind_rows(
  fit_mld_reg("SST.GOA", "GOA.sd", "GOA.ar1", "GOA"),
  fit_mld_reg("SST.EBS", "EBS.sd", "EBS.ar1", "EBS"))

write.csv(reg.mld, "./Output/sst_ar1_on_mld_sd_ar1_regression.csv", row.names = FALSE)
message("Saved: Output/sst_ar1_on_mld_sd_ar1_regression.csv")
message("SST AR(1) ~ MLD SD + MLD AR(1) (GLS AR(1); standardized predictors):")
print(reg.mld)

p.reg.mld <- reg.mld %>%
  filter(term != "(Intercept)") %>%
  mutate(term = recode(term, "MLD.sd" = "MLD SD", "MLD.ar1" = "MLD AR(1)")) %>%
  ggplot(aes(x = beta, y = term, color = system)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(xmin = beta - 1.96 * se, xmax = beta + 1.96 * se),
                  position = position_dodge(width = 0.4), linewidth = 0.7) +
  scale_color_manual(values = c("GOA" = "darkorange4", "EBS" = "steelblue4")) +
  labs(title = "SST AR(1) ~ MLD SD + MLD AR(1)  (GLS, AR(1) residuals)",
       subtitle = "Standardized predictors; bars = 95% CI",
       x = "Standardized coefficient (per 1 SD predictor)", y = NULL, color = "System") +
  theme_bw(base_size = 11)

ggsave("./Figures/SST_AR1_on_MLD_SD_AR1_coefs.png",
       plot = p.reg.mld, width = 7, height = 4, dpi = 300)
message("Saved: Figures/SST_AR1_on_MLD_SD_AR1_coefs.png")

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
    roll_right(x, width, FUN = safe_ar1)
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
    roll_right(x, width, FUN = safe_ar1)
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
    roll_right(x, width, FUN = safe_ar1)
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
    roll_right(x, width, FUN = safe_mean)
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
    roll_right(x, width, FUN = safe_mean)
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
    roll_right(x, width, FUN = safe_mean)
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
    roll_right(x, width, FUN = safe_ar1)
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
  coord_rect(xlim = c(110, 250), ylim = c(20, 66)) +
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
    roll_right(x, width, FUN = safe_ar1)
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

# --- Regrid to 0.5 deg for plotting ---
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

# --- Polar stereographic rendering (GRL manuscript style) ------------------
# Layout, projection and domain all live in ps_regression_map() so this panel
# and the other manuscript maps cannot drift apart.

p.sst.ar1.reg <- ps_regression_map(
  sst.ar1.reg.grid,
  fill.name = expression(beta~"(AR(1) per SD)"),
  col.lim   = col.lim.sst)

ggsave("./Figures/SST_AR1_AL_SLP_SD_regression.png",
       plot = p.sst.ar1.reg, width = grl.width.full, height = 4.8,
       dpi = 400, bg = "white")
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
sv.sst   <- partial_svd(sweep(sst.ok, 1, w.sst, "*"), nv = 3)

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

# ============================================================
# SECTION 21: MANUSCRIPT FIGURES 1-2 — study site map + SST time series
# ============================================================
# Fig. 1: a standalone study-site map, drawn through ps_regression_map() so
# it shares the projection and domain of Figs. 4a/5a. It shows Section 4's
# SLP EOF1 loadings with the EBS and GOA polygons and the AL Index box.
# Fig. 2: GOA and EBS winter SST and its rolling AR(1).

# --- Figure 1: study site map -----------------------------------------------
# Section 4's EOF1 is on the 0.25-deg 120-250E grid. It is thinned to the
# 0.5-deg cells of Figs. 4a/5a and cut to the map domain; the field is smooth,
# so thinning loses nothing visible. SLP is an atmospheric field, so land is
# left unfilled and the loadings show over it.
# Colours follow Findings_summary Fig. 4: RdBu stretched over the data range
# (about -0.001 to 0.011), so white sits near mid-range rather than at zero.
# Blue therefore marks weak loadings, not negative ones -- the field is
# one-signed apart from the far southwest -- which the caption should state.
fig1.eof <- eof.df.slp %>%
  filter(lon >= 160, lon <= 250, lat >= 20, lat <= 66,
         lon %% 0.5 == 0, lat %% 0.5 == 0) %>%
  rename(beta = EOF1)

# Share of the (sqrt-cos-weighted, detrended) monthly SLP variance in EOF1,
# for the caption
fig1.var <- sv.slp$d[1]^2 / sum(rowSums(detrend.slp^2) * w.slp^2)
message(sprintf("Fig. 1: SLP EOF1 explains %.1f%% of variance", 100 * fig1.var))

p.fig1 <- ps_regression_map(fig1.eof,
                            fill.name  = "SLP EOF1 loading",
                            col.lim    = max(abs(fig1.eof$beta)),
                            label.size = 6 / .pt,
                            land.fill  = NA,
                            fill.scale = scale_fill_distiller(
                              palette = "RdBu", direction = -1,
                              name = "SLP EOF1 loading"))

# GRL single-column width, so it prints unscaled. Height is the map's
# projected aspect at that width plus room for the colour bar, so the page
# carries no slack around the map.
fig1.map.w <- grl.width.single - 0.16   # less the 2-mm plot margins
fig1.map.h <- fig1.map.w * diff(p.fig1$coordinates$limits$y) /
                           diff(p.fig1$coordinates$limits$x)
for (ext in c("png", "pdf"))
  ggsave(paste0("./Figures/Fig_1.", ext), plot = p.fig1,
         width = grl.width.single, height = fig1.map.h + 0.55,
         dpi = 400, bg = "white")
message("Saved: Figures/Fig_1.png and .pdf")

# --- Figure 2: GOA and EBS winter SST and rolling AR(1) --------------------
# Region means are recomputed here rather than reused from Section 2, which
# masks with sf::st_within. This uses a dependency-free ray-casting mask and
# reads only the GOA/EBS bounding box out of the ERA5 file, so the whole
# panel set builds from a small cache.

fig2.sst.cache <- "./Output/era5_sst_goa_ebs_monthly.csv"

if (file.exists(fig2.sst.cache)) {
  message("Loading cached GOA/EBS monthly SST: ", fig2.sst.cache)
  fig2.monthly <- read.csv(fig2.sst.cache)
} else {
  sst.f <- if (file.exists("./Data/era5_sst_NP_1950_2026.nc"))
    "./Data/era5_sst_NP_1950_2026.nc" else "./Data/era5_sst_NP_1950_2025.nc"
  message("Computing GOA/EBS monthly SST means from: ", sst.f)
  nc.f  <- nc_open(sst.f)
  lon.f <- ncvar_get(nc.f, "longitude"); lat.f <- ncvar_get(nc.f, "latitude")
  tv    <- intersect(c("valid_time", "time"),
                     c(names(nc.f$var), names(nc.f$dim)))[1]
  t.f   <- ncvar_get(nc.f, tv)

  # Read only the block covering both polygons
  li <- which(lon.f >= 182 & lon.f <= 235)
  ai <- which(lat.f >= 49  & lat.f <= 66)
  sub <- ncvar_get(nc.f, "sst",
                   start = c(min(li), min(ai), 1),
                   count = c(length(li), length(ai), -1)) - 273.15
  nc_close(nc.f)

  # as.vector: ncvar_get returns 1-D arrays, and a dim attribute on the
  # weights makes the matrix arithmetic below non-conformable.
  lon.s <- as.vector(lon.f[li]); lat.s <- as.vector(lat.f[ai])
  g   <- expand.grid(lon = lon.s, lat = lat.s)   # lon varies fastest
  m.goa <- inpoly_vec(g$lon, g$lat, goa.poly$lon, goa.poly$lat)
  m.ebs <- inpoly_vec(g$lon, g$lat, ebs.poly$lon, ebs.poly$lat)
  message(sprintf("  mask cells: GOA %d, EBS %d", sum(m.goa), sum(m.ebs)))

  mat <- matrix(sub, nrow = length(lon.s) * length(lat.s))
  rm(sub); gc()
  wt  <- as.vector(cos(g$lat * pi / 180))

  # Note the parentheses around (!is.na(x)): `!` binds looser than `*` in R,
  # so `!is.na(x) * w` would negate the product and give a wrong denominator.
  area_mean <- function(mask) {
    x <- mat[mask, , drop = FALSE]; w <- wt[mask]
    colSums(x * w, na.rm = TRUE) / colSums((!is.na(x)) * w)
  }
  dts <- as.Date(as.POSIXct(t.f, origin = "1970-01-01", tz = "UTC"))
  fig2.monthly <- data.frame(year  = as.integer(format(dts, "%Y")),
                             month = as.integer(format(dts, "%m")),
                             GOA   = area_mean(m.goa),
                             EBS   = area_mean(m.ebs))
  rm(mat); gc()
  write.csv(fig2.monthly, fig2.sst.cache, row.names = FALSE)
  message("Saved: ", fig2.sst.cache)
}

# --- Winter (Nov-Mar, year = January) mean SST, in degrees C ------------
# Raw SST rather than anomalies: panels (a)/(c) show absolute winter means.
# AR(1) is unaffected by the choice -- a winter anomaly differs from a raw
# winter mean by a constant (the same five climatological months every
# year), and detrending removes constants -- so panels (b)/(d) are identical
# either way. Years without all five months are dropped, so the partial
# winter of 1950 and the partial 2026 do not masquerade as full seasons.
fig2.winter <- fig2.monthly %>%
  mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
  filter(month %in% c(11, 12, 1, 2, 3)) %>%
  group_by(year = win.year) %>%
  filter(n() == 5) %>%
  summarise(GOA = mean(GOA), EBS = mean(EBS), .groups = "drop")

# --- Rolling AR(1) on the detrended winter series (Section 3 method) -----
fig2.ar1 <- function(d, col) {
  x <- detrend_lin(d[[col]])
  data.frame(year = d$year, ar1 = roll_ar1_right(x))
}

fig2.ts <- fig2.winter %>%
  tidyr::pivot_longer(c(GOA, EBS), names_to = "region", values_to = "sst")

fig2.roll <- bind_rows(
  cbind(fig2.ar1(fig2.winter, "GOA"), region = "GOA"),
  cbind(fig2.ar1(fig2.winter, "EBS"), region = "EBS"))

# Colour encodes region, matching the GOA/EBS colours used in Figs. 3-5.
# Series colour is set directly rather than mapped, so the figure carries no
# legend; the row headers name the ecosystems. Tags are set per panel with
# labs(tag=) so the spanning row headers stay untagged.
fig2.zero <- function() geom_hline(yintercept = 0, colour = "grey70",
                                   linewidth = 0.25)
# Tags take the default top-left slot above each panel, as in Figs. 4-5;
# inside the plot area they collide with the top y tick label
fig2.panel.theme <- theme(plot.margin = margin(0.5, 1, 0.5, 1, "mm"),
                          plot.tag    = element_text(size = 9, face = "bold"))

ts_panel <- function(reg, tag) {
  ggplot(filter(fig2.ts, region == reg), aes(x = year, y = sst)) +
    geom_line(linewidth = 0.4, colour = area.colors[[reg]]) +
    geom_point(size = 0.5, colour = area.colors[[reg]]) +
    labs(x = NULL, y = "SST (°C)", tag = tag) +
    theme_grl() + fig2.panel.theme
}

ar1_panel <- function(reg, tag) {
  ggplot(filter(fig2.roll, region == reg, !is.na(ar1)), aes(x = year, y = ar1)) +
    fig2.zero() +
    geom_line(linewidth = 0.4, colour = area.colors[[reg]]) +
    geom_point(size = 0.5, colour = area.colors[[reg]]) +
    labs(x = "Year", y = "AR(1), 15-yr window", tag = tag) +
    theme_grl() + fig2.panel.theme
}

p.fig2a <- ts_panel("GOA", "(a)")
p.fig2b <- ar1_panel("GOA", "(b)")
p.fig2c <- ts_panel("EBS", "(c)")
p.fig2d <- ar1_panel("EBS", "(d)")

# Row headers span both columns, so each ecosystem name sits centred over its
# pair of panels rather than over the left-hand one only.
fig2.hdr <- function(txt)
  wrap_elements(full = grid::textGrob(txt, gp = grid::gpar(fontsize = 9)))

# GOA on the top pair and EBS on the bottom, each under its own spanning
# header. Built at GRL's full page width so it prints unscaled; the height
# keeps the panels' earlier ~0.75 height:width.
fig2.design <- "
EE
AB
FF
CD
"

p.fig2 <- p.fig2a + p.fig2b + p.fig2c + p.fig2d +
  fig2.hdr("Gulf of Alaska") + fig2.hdr("Eastern Bering Sea") +
  plot_layout(design = fig2.design, heights = c(0.04, 1, 0.04, 1))

for (ext in c("png", "pdf"))
  ggsave(paste0("./Figures/Fig_2.", ext), plot = p.fig2,
         width = grl.width.full, height = 6.1, dpi = 400, bg = "white")
message("Saved: Figures/Fig_2.png and .pdf")

# ============================================================
# SECTION 22: BOREAL LME WINTER SST AR(1) — DYNAMIC FACTOR ANALYSIS
# ============================================================
# Is the reddening seen in the GOA and EBS shared across boreal LMEs?
# Re-runs the LME DFA (Figures/lme.AR1.DFA.png) with the manuscript AR(1)
# in place of annual means and centred windows: complete Nov-Mar winters
# (year = January), linear detrend, then 15-yr right-aligned rolling AR(1),
# i.e. the Section 21 (Fig. 2) method. As in that figure, the DFA is fit to
# the LMEs alone, so comparing the GOA and EBS series with the shared trend
# is not circular. The LMEs are those in ERA5_processing_LMEs.R: LMEs north
# of 40N less six high-Arctic LMEs (which the original figure still had),
# and less the Black Sea, dropped on 2026-10-24.
#
# Ice check: ERA5 holds SST at 271.46 K (-1.69 C) under sea ice, so winter
# AR(1) in ice-covered LMEs partly tracks ice cover. The DFA is re-run
# without the LMEs that extend north of the EBS polygon (65N), and each
# series' winter ice fraction (share of ocean area at 271.46 K) is listed.
#
# ./Data/era5_sst_LME_1950_2026.nc (downloads.R Section 5) and the LME
# shapefile are read only when the monthly cache below is missing. GOA/EBS
# SST comes from the Section 21 cache, so their AR(1) matches Fig. 2.

lme.numbers <- c(1, 2, 8, 9, 18, 19, 20, 21, 22, 23, 24, 51, 52, 53, 54,
                 55, 56, 57, 58, 59, 60, 62, 63, 64)
lme.drop    <- c("Laptev Sea", "Kara Sea", "Canadian Eastern Arctic - West Greenland",
                 "Beaufort Sea", "Central Arctic", "East Siberian Sea",
                 "East Bering Sea", "Gulf of Alaska",   # EBS/GOA: study polygons
                 "Black Sea")                           # dropped 2026-10-24
ice.sst.K   <- 271.465           # ERA5 sea-ice SST is 271.46 K; margin for float noise
lme.north   <- max(ebs.poly$lat) # 65N, the northern edge of the EBS polygon

# Minimal polygon shapefile reader (no sf): one data frame of vertices per
# record, with a ring id so multi-part LMEs and holes can be masked.
read_shp_rings <- function(shp) {
  con <- file(shp, "rb"); on.exit(close(con))
  invisible(readBin(con, "raw", 100))                          # file header
  out <- list()
  repeat {
    hdr <- readBin(con, "integer", 2, size = 4, endian = "big")
    if (length(hdr) < 2) break
    rc   <- rawConnection(readBin(con, "raw", 2 * hdr[2]))    # length in 16-bit words
    type <- readBin(rc, "integer", 1, size = 4, endian = "little")
    if (type == 5) {                                           # polygon
      invisible(readBin(rc, "double", 4, size = 8, endian = "little"))  # bbox
      n     <- readBin(rc, "integer", 2, size = 4, endian = "little")   # parts, points
      parts <- readBin(rc, "integer", n[1], size = 4, endian = "little")
      xy    <- readBin(rc, "double", 2 * n[2], size = 8, endian = "little")
      out[[length(out) + 1]] <- data.frame(x = xy[c(TRUE, FALSE)],
                                           y = xy[c(FALSE, TRUE)],
                                           ring = findInterval(seq_len(n[2]) - 1, parts))
    } else {
      out[length(out) + 1] <- list(NULL)
    }
    close(rc)
  }
  out
}

# Even-odd rule across all rings, so a point in a hole sits inside two rings
in_rings <- function(x, y, rings) {
  inside <- rep(FALSE, length(x))
  for (r in split(rings, rings$ring)) {
    cand <- x >= min(r$x) & x <= max(r$x) & y >= min(r$y) & y <= max(r$y)
    if (any(cand))
      inside[cand] <- xor(inside[cand], inpoly_vec(x[cand], y[cand], r$x, r$y))
  }
  inside
}

# --- Monthly area-weighted SST and ice fraction per region (cached) ------
lme.sst.cache <- "./Output/era5_sst_lme_monthly.csv"

if (file.exists(lme.sst.cache)) {
  message("Loading cached LME monthly SST: ", lme.sst.cache)
  lme.monthly <- read.csv(lme.sst.cache)
} else {
  message("Computing LME monthly SST from ./Data/era5_sst_LME_1950_2026.nc")
  lme.dbf   <- foreign::read.dbf("./Data/LME shapefiles/LMEs66.dbf", as.is = TRUE)
  lme.rings <- read_shp_rings("./Data/LME shapefiles/LMEs66.shp")
  lme.keep  <- which(lme.dbf$LME_NUMBER %in% lme.numbers &
                     !lme.dbf$LME_NAME %in% lme.drop)

  nc.l  <- nc_open("./Data/era5_sst_LME_1950_2026.nc")
  lon.l <- as.vector(ncvar_get(nc.l, "longitude"))   # -180 to 179.5, as the shapefile
  lat.l <- as.vector(ncvar_get(nc.l, "latitude"))
  t.l   <- ncvar_get(nc.l, "valid_time")
  mat.l <- matrix(ncvar_get(nc.l, "sst"), nrow = length(lon.l) * length(lat.l))
  nc_close(nc.l)

  g.l  <- expand.grid(lon = lon.l, lat = lat.l)      # lon varies fastest, as sst
  wt.l <- cos(g.l$lat * pi / 180)
  dt.l <- as.Date(as.POSIXct(t.l, origin = "1970-01-01", tz = "UTC"))

  # Same NA-safe weighting as Section 21's area_mean(); ice is the share of
  # ocean area sitting at ERA5's sea-ice SST
  region_stats <- function(mask, region, lat.max) {
    x <- mat.l[mask, , drop = FALSE]; w <- wt.l[mask]
    ocean <- colSums((!is.na(x)) * w)
    data.frame(region = region, lat.max = lat.max, n.cells = sum(!is.na(x[, 1])),
               year  = as.integer(format(dt.l, "%Y")),
               month = as.integer(format(dt.l, "%m")),
               sst   = colSums(x * w, na.rm = TRUE) / ocean - 273.15,
               ice   = colSums((x <= ice.sst.K) * w, na.rm = TRUE) / ocean)
  }

  lme.monthly <- bind_rows(
    lapply(lme.keep, function(k)
      region_stats(in_rings(g.l$lon, g.l$lat, lme.rings[[k]]),
                   lme.dbf$LME_NAME[k], max(lme.rings[[k]]$y))),
    # Study polygons on the same 0.5-deg grid: their ice fraction, and a
    # check of this pipeline against Fig. 2's AR(1)
    region_stats(inpoly_vec(to360(g.l$lon), g.l$lat, ebs.poly$lon, ebs.poly$lat),
                 "EBS", max(ebs.poly$lat)),
    region_stats(inpoly_vec(to360(g.l$lon), g.l$lat, goa.poly$lon, goa.poly$lat),
                 "GOA", max(goa.poly$lat)))
  rm(mat.l); gc()

  write.csv(lme.monthly, lme.sst.cache, row.names = FALSE)
  message("Saved: ", lme.sst.cache)
}

# The drop list is applied again here so that it also holds for a cache
# written before an LME was dropped (the cache still carries the Black Sea)
lme.monthly <- filter(lme.monthly, !region %in% lme.drop)
lme.meta    <- distinct(lme.monthly, region, lat.max, n.cells)

# --- Winter means and rolling AR(1), as in Section 21 --------------------
winter_means <- function(d) {
  d %>%
    filter(month %in% c(11, 12, 1, 2, 3)) %>%
    mutate(year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
    group_by(region, year) %>%
    filter(n() == 5) %>%                     # complete winters only
    summarise(across(any_of(c("sst", "ice")), mean), .groups = "drop")
}

winter_ar1 <- function(d) {
  d %>%
    group_by(region) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(ar1 = roll_ar1_right(detrend_lin(sst))) %>%
    ungroup()
}

lme.winter <- winter_means(lme.monthly)

goa.ebs.winter <- read.csv("./Output/era5_sst_goa_ebs_monthly.csv") %>%
  pivot_longer(c(GOA, EBS), names_to = "region", values_to = "sst") %>%
  winter_means()

# DFA input is the LMEs only; GOA/EBS AR(1), from Fig. 2's series, is
# compared with the fitted trend afterwards
lme.ar1   <- winter_ar1(filter(lme.winter, !region %in% c("EBS", "GOA")))
study.ar1 <- winter_ar1(goa.ebs.winter)

# Check: the 0.5-deg pipeline should reproduce Fig. 2's GOA/EBS AR(1)
lme.pipe.check <- inner_join(
  winter_ar1(filter(lme.winter, region %in% c("EBS", "GOA"))),
  study.ar1,
  by = c("region", "year"), suffix = c(".lme", ".fig2")) %>%
  filter(!is.na(ar1.lme), !is.na(ar1.fig2)) %>%
  group_by(region) %>%
  summarise(r = cor(ar1.lme, ar1.fig2),
            max.abs.diff = max(abs(ar1.lme - ar1.fig2)), .groups = "drop")
message("Check, 0.5-deg LME pipeline vs Fig. 2 GOA/EBS AR(1):")
print(lme.pipe.check)

# --- Winter ice fraction per region --------------------------------------
lme.ice <- lme.winter %>%
  group_by(region) %>%
  summarise(ice.mean  = mean(ice),
            ice.51.80 = mean(ice[year <= 1980]),
            ice.97.26 = mean(ice[year >= 1997]), .groups = "drop") %>%
  left_join(lme.meta, by = "region") %>%
  mutate(north.of.ebs = lat.max > lme.north) %>%
  arrange(desc(ice.mean))
write.csv(lme.ice, "./Output/LME_winter_ice_fraction.csv", row.names = FALSE)
message("Winter ice fraction (share of ocean area at ERA5 sea-ice SST):")
print(as.data.frame(mutate(lme.ice, across(starts_with("ice"), ~ round(.x, 3)))))

# --- One-trend DFA, all LMEs and without those north of 65N --------------
ar1.wide <- lme.ar1 %>%
  filter(!is.na(ar1)) %>%
  select(region, year, ar1) %>%
  pivot_wider(names_from = year, values_from = ar1)
ar1.mat <- as.matrix(ar1.wide[, -1])
rownames(ar1.mat) <- ar1.wide$region

# Same model, R-matrix candidates and EM settings as calculate.ar1.sd.LME.R
dfa.control <- list(minit = 200, maxit = 30000, allow.degen = FALSE,
                    conv.test.slope.tol = 0.1, abstol = 0.0001)

fit_dfa_1trend <- function(y, run) {
  R.forms <- c("diagonal and equal", "diagonal and unequal", "equalvarcov")
  fits <- lapply(R.forms, function(R)
    MARSS::MARSS(y, model = list(A = "zero", R = R, m = 1), form = "dfa",
                 z.score = TRUE, control = dfa.control, silent = TRUE))
  sel <- data.frame(run = run, R = R.forms,
                    AICc = sapply(fits, function(f) f$AICc),
                    converged = sapply(fits, function(f) f$convergence == 0),
                    min.eig.R = sapply(fits, function(f)
                      min(eigen(coef(f, type = "matrix")$R, symmetric = TRUE,
                                only.values = TRUE)$values)))
  sel$dAICc <- sel$AICc - min(sel$AICc)
  # Boundary fits pin the trend to the data instead of estimating it: in
  # "diagonal and unequal" one series' error variance goes to ~0 and the
  # trend becomes that series; in "equalvarcov" the error correlation can
  # reach its floor of -1/(n-1), making the cross-series mean noise-free.
  # Either way the error covariance R is (near) singular and the trend SE
  # collapses, and such a trend cannot show whether reddening is shared.
  # AICc is therefore compared only among converged fits whose R has all
  # eigenvalues above 0.01 (data are z-scored).
  sel$usable <- sel$converged & sel$min.eig.R > 0.01
  if (!any(sel$usable)) stop(run, ": no usable DFA fit")
  best.i <- which(sel$usable)[which.min(sel$AICc[sel$usable])]
  sel$selected <- seq_along(R.forms) == best.i
  best <- fits[[best.i]]
  ci <- MARSS::MARSSparamCIs(best)
  # With one trend, loadings and trend share an arbitrary sign: orient the
  # trend with the majority of series (loadings summing positive)
  s  <- if (sum(ci$par$Z) < 0) -1 else 1
  lo <- s * as.vector(ci$par.lowCI$Z); up <- s * as.vector(ci$par.upCI$Z)
  list(selection = sel,
       loadings  = data.frame(run = run, region = rownames(y),
                              loading = s * as.vector(ci$par$Z),
                              lower = pmin(lo, up), upper = pmax(lo, up)),
       trend     = data.frame(run = run, year = as.integer(colnames(y)),
                              trend = s * as.vector(best$states),
                              se = as.vector(best$states.se)))
}

run.labels <- c(all   = "All LMEs",
                south = "Excluding LMEs extending north of 65°N")
lme.south  <- lme.meta$region[lme.meta$lat.max <= lme.north]

message("Fitting DFA: all LMEs (", nrow(ar1.mat), " series)...")
dfa.all <- fit_dfa_1trend(ar1.mat, run.labels[["all"]])
message("Fitting DFA: without LMEs north of 65N (",
        sum(rownames(ar1.mat) %in% lme.south), " series)...")
dfa.south <- fit_dfa_1trend(ar1.mat[rownames(ar1.mat) %in% lme.south, ],
                            run.labels[["south"]])

# A trend's sign is arbitrary, and with loadings of both signs the
# sum-of-loadings rule can orient the two runs oppositely, so the ice-check
# run is aligned with the all-LMEs run before they are compared
flip_dfa <- function(d) {
  lo <- d$loadings$lower
  d$loadings$lower   <- -d$loadings$upper
  d$loadings$upper   <- -lo
  d$loadings$loading <- -d$loadings$loading
  d$trend$trend      <- -d$trend$trend
  d
}
if (cor(dfa.all$trend$trend, dfa.south$trend$trend) < 0)
  dfa.south <- flip_dfa(dfa.south)

r.trend <- cor(dfa.all$trend$trend, dfa.south$trend$trend)
message(sprintf("Shared trend, all LMEs vs without north of 65N: r = %.2f", r.trend))
print(bind_rows(dfa.all$selection, dfa.south$selection))

# --- GOA/EBS AR(1) vs the shared LME trend -------------------------------
# Modified Chelton effective sample size from Section 5, since both series
# are heavily autocorrelated 15-yr rolling statistics. The test is
# two-sided because the trend's sign is arbitrary.
study_vs_trend <- function(dfa) {
  study.ar1 %>%
    inner_join(dfa$trend, by = "year") %>%
    group_by(run, region) %>%
    summarise(ct = list(chelton_test(ar1, trend)), .groups = "drop") %>%
    mutate(r     = sapply(ct, `[[`, "r.obs"),
           n.eff = sapply(ct, `[[`, "n.eff"),
           p     = sapply(ct, function(x) 2 * pt(-abs(x$t.obs), x$df))) %>%
    select(-ct)
}
lme.study.cor <- bind_rows(study_vs_trend(dfa.all), study_vs_trend(dfa.south))
message("GOA/EBS AR(1) vs shared LME trend (modified Chelton, two-sided):")
print(as.data.frame(lme.study.cor))
write.csv(lme.study.cor, "./Output/LME_winter_AR1_DFA_EBS_GOA_correlation.csv",
          row.names = FALSE)

write.csv(bind_rows(study.ar1, lme.ar1), "./Output/LME_winter_AR1_series.csv",
          row.names = FALSE)
write.csv(bind_rows(dfa.all$selection, dfa.south$selection),
          "./Output/LME_winter_AR1_DFA_model_selection.csv", row.names = FALSE)
write.csv(bind_rows(dfa.all$loadings, dfa.south$loadings),
          "./Output/LME_winter_AR1_DFA_loadings.csv", row.names = FALSE)
write.csv(bind_rows(dfa.all$trend, dfa.south$trend),
          "./Output/LME_winter_AR1_DFA_trend.csv", row.names = FALSE)
message("Saved: Output/LME_winter_AR1_*.csv, Output/LME_winter_ice_fraction.csv")

# --- Figures ---------------------------------------------------------------
# LMEs are neutral grey; GOA/EBS keep their Fig. 2 colours.
lme_label <- function(r) recode(r, EBS = "Eastern Bering Sea", GOA = "Gulf of Alaska")
lme.order <- lme_label(dfa.all$loadings$region[order(dfa.all$loadings$loading)])
lme.panel.theme <- theme(plot.tag = element_text(size = 9, face = "bold"))

# DFA trend units are arbitrary, so trends are shown standardized, with the
# SE scaled by the same factor
std_trend <- function(tr) {
  tr %>%
    group_by(run) %>%
    mutate(k = sd(trend), z = (trend - mean(trend)) / k,
           lo = z - 1.96 * se / k, hi = z + 1.96 * se / k) %>%
    ungroup()
}
fmt_p <- function(p) ifelse(p < 0.01, "p < 0.01", sprintf("p = %.2f", p))

dfa_loadings_panel <- function(ld, tag) {
  ld %>%
    mutate(label = factor(lme_label(region), levels = lme.order)) %>%
    ggplot(aes(x = label, y = loading)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
    geom_col(fill = "grey65", width = 0.7) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, linewidth = 0.3) +
    coord_flip() +
    labs(x = NULL, y = "Loading on shared trend", tag = tag) +
    theme_grl() + lme.panel.theme
}

# Shared LME trend with the GOA/EBS AR(1) series laid over it, all
# standardized; the legend carries each series' r with the trend
dfa_overlay_panel <- function(dfa, tag) {
  tr <- std_trend(dfa$trend)
  st <- study.ar1 %>%
    filter(year %in% tr$year) %>%
    group_by(region) %>%
    mutate(z = (ar1 - mean(ar1)) / sd(ar1)) %>%
    ungroup() %>%
    mutate(series = lme_label(region))
  cr   <- filter(lme.study.cor, run == dfa$trend$run[1])
  cols <- c("Shared LME trend"   = "black",
            "Eastern Bering Sea" = area.colors[["EBS"]],
            "Gulf of Alaska"     = area.colors[["GOA"]])
  labs.r <- setNames(sprintf("%s (r = %.2f, %s)", lme_label(cr$region), cr$r,
                             fmt_p(cr$p)), lme_label(cr$region))
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
    geom_ribbon(data = tr, aes(x = year, ymin = lo, ymax = hi), fill = "grey85") +
    geom_line(data = st, aes(x = year, y = z, colour = series), linewidth = 0.4) +
    geom_line(data = tr, aes(x = year, y = z, colour = "Shared LME trend"),
              linewidth = 0.6) +
    scale_colour_manual(values = cols, breaks = names(cols), name = NULL,
                        labels = c("Shared LME trend" = "Shared LME trend (±1.96 SE)",
                                   labs.r[names(cols)[-1]])) +
    labs(x = "Year (end of 15-yr window)", y = "Winter SST AR(1), standardized",
         tag = tag) +
    theme_grl() + lme.panel.theme
}

p.lme.dfa <- (dfa_loadings_panel(dfa.all$loadings, "(a)") |
                dfa_overlay_panel(dfa.all, "(b)")) +
  plot_layout(widths = c(1, 1.15), guides = "collect") &
  theme(legend.position = "bottom")

ggsave("./Figures/LME_winter_AR1_DFA.png", plot = p.lme.dfa,
       width = grl.width.full, height = 3.4, dpi = 400, bg = "white")
message("Saved: Figures/LME_winter_AR1_DFA.png")

# Ice check: both runs' loadings and trends on shared axes
run.colors <- setNames(c("#2a78d6", "#eb6834"), run.labels)

p.ice.a <- bind_rows(dfa.all$loadings, dfa.south$loadings) %>%
  mutate(label = factor(lme_label(region), levels = lme.order),
         run   = factor(run, levels = run.labels)) %>%
  ggplot(aes(x = label, y = loading, colour = run)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.6),
                  size = 0.12, linewidth = 0.3) +
  scale_colour_manual(values = run.colors, name = NULL, guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "Loading on shared trend", tag = "(a)") +
  theme_grl() + lme.panel.theme

p.ice.b <- std_trend(bind_rows(dfa.all$trend, dfa.south$trend)) %>%
  mutate(run = factor(run, levels = run.labels)) %>%
  ggplot(aes(x = year, y = z, colour = run, fill = run)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.4) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.6, size = 7 / .pt,
           label = sprintf("r = %.2f", r.trend)) +
  scale_colour_manual(values = run.colors, name = NULL) +
  scale_fill_manual(values = run.colors, name = NULL) +
  labs(x = "Year (end of 15-yr window)", y = "Shared trend, standardized",
       tag = "(b)") +
  theme_grl() + lme.panel.theme

p.lme.ice <- (p.ice.a | p.ice.b) +
  plot_layout(widths = c(1, 1.15), guides = "collect") &
  theme(legend.position = "bottom")

ggsave("./Figures/LME_winter_AR1_DFA_ice_check.png", plot = p.lme.ice,
       width = grl.width.full, height = 3.8, dpi = 400, bg = "white")
message("Saved: Figures/LME_winter_AR1_DFA_ice_check.png")

# Diagnostic: GOA/EBS, then the DFA inputs ordered by loading, each with its
# winter ice fraction
wrap_name <- function(x, width = 22)
  vapply(strwrap(x, width, simplify = FALSE), paste, "", collapse = "\n")

lme.series.dat <- bind_rows(study.ar1, lme.ar1) %>%
  filter(!is.na(ar1)) %>%
  left_join(select(lme.ice, region, ice.mean, north.of.ebs), by = "region") %>%
  mutate(strip = paste0(wrap_name(lme_label(region)), "\nice ",
                        round(100 * ice.mean), "%",
                        ifelse(north.of.ebs, ", >65°N", "")),
         col   = ifelse(region %in% names(area.colors),
                        area.colors[region], "grey20"))
lme.series.dat$strip <- factor(
  lme.series.dat$strip,
  levels = unique(lme.series.dat$strip[order(match(
    lme_label(lme.series.dat$region),
    c("Eastern Bering Sea", "Gulf of Alaska", rev(lme.order))))]))

p.lme.series <- ggplot(lme.series.dat, aes(x = year, y = ar1)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_line(aes(colour = col), linewidth = 0.35) +
  scale_colour_identity() +
  facet_wrap(~ strip, ncol = 6) +
  labs(x = "Year (end of 15-yr window)", y = "Winter SST AR(1), 15-yr window") +
  theme_grl(base_size = 7)

ggsave("./Figures/LME_winter_AR1_series.png", plot = p.lme.series,
       width = grl.width.full, height = 5.0, dpi = 400, bg = "white")
message("Saved: Figures/LME_winter_AR1_series.png")

# ============================================================
# SECTION 23: CMIP6 WINTER SST AR(1) — SSP2-4.5 vs PREINDUSTRIAL CONTROL
# ============================================================
# Is increased reddening an expression of climate change? Replots
# CMIP6SST.AR1.png (CMIP6.ar1.sd.R) with the manuscript AR(1): Nov-Mar means
# (year = January) over this study's GOA and EBS polygons, cos(lat)
# weighted, then 15-yr right-aligned rolling AR(1). Each of the 20 SSP2-4.5
# models is paired with its own piControl run (stored in its SSP5-8.5 file)
# and both are processed identically; models are averaged with the
# region-specific skill weights from CMIP6-attribution.
#
# As for the observed series (Section 21), each model series is linearly
# detrended over its full length (winters 1851-2099) before the rolling
# AR(1). Forced warming is curved over 250 years, but what a straight line
# leaves within 15-yr windows barely matters: the loess detrend used in
# CMIP6.ar1.sd.R changes the weighted mean AR(1) by at most 0.02.
# piControl years (1850-2099) are labels only.
#
# ./Data/CMIP6 (downloads.R Section 6) is read only when the monthly cache
# below is missing. The GOA polygon differs from the smaller coastal GOA box
# used to derive the model weights; the EBS polygons are identical.

cmip6.dir   <- "./Data/CMIP6"
cmip6.cache <- "./Output/cmip6_sst_ebs_goa_monthly.csv"

if (file.exists(cmip6.cache)) {
  message("Loading cached CMIP6 GOA/EBS monthly SST: ", cmip6.cache)
  cmip6.monthly <- read.csv(cmip6.cache)
} else {
  message("Computing CMIP6 GOA/EBS monthly SST from ", cmip6.dir)
  cmip6.models <- sub("_245\\.nc$", "", list.files(file.path(cmip6.dir, "ssp245"),
                                                    pattern = "_245\\.nc$"))
  cmip6.monthly <- bind_rows(lapply(cmip6.models, function(model) {
    runs <- c(hist_ssp245 = file.path(cmip6.dir, "ssp245", paste0(model, "_245.nc")),
              piControl   = file.path(cmip6.dir, "ssp585", paste0(model, ".nc")))
    bind_rows(lapply(names(runs), function(ex) {
      nc  <- nc_open(runs[[ex]])
      j   <- match(ex, ncvar_get(nc, "experiment"))
      if (is.na(j)) stop(ex, " not found in ", runs[[ex]])
      lon <- as.vector(ncvar_get(nc, "lon")); lat <- as.vector(ncvar_get(nc, "lat"))
      dts <- as.Date(ncvar_get(nc, "time"), origin = "1970-01-01")
      tos <- ncvar_get(nc, "tos", start = c(j, 1, 1, 1), count = c(1, -1, -1, -1))
      nc_close(nc)
      if (max(tos, na.rm = TRUE) > 200) tos <- tos - 273.15   # K -> C

      g   <- expand.grid(lon = lon, lat = lat)   # lon varies fastest, as tos
      mat <- matrix(tos, nrow = nrow(g))
      wt  <- cos(g$lat * pi / 180)
      # Same NA-safe area weighting as Section 21's area_mean()
      area_mean <- function(mask) {
        x <- mat[mask, , drop = FALSE]; w <- wt[mask]
        colSums(x * w, na.rm = TRUE) / colSums((!is.na(x)) * w)
      }
      data.frame(model = model, experiment = ex,
                 year  = as.integer(format(dts, "%Y")),
                 month = as.integer(format(dts, "%m")),
                 EBS = area_mean(inpoly_vec(g$lon, g$lat, ebs.poly$lon, ebs.poly$lat)),
                 GOA = area_mean(inpoly_vec(g$lon, g$lat, goa.poly$lon, goa.poly$lat)))
    }))
  }))
  write.csv(cmip6.monthly, cmip6.cache, row.names = FALSE)
  message("Saved: ", cmip6.cache)
}

# Check: every run should be 3000 consecutive months, Jan 1850 - Dec 2099
cmip6.span <- cmip6.monthly %>%
  group_by(model, experiment) %>%
  summarise(ok = n() == 3000 && all(diff(12 * year + month) == 1) &&
              min(year) == 1850 && max(year) == 2099, .groups = "drop")
if (!all(cmip6.span$ok))
  stop("CMIP6 runs with unexpected time axes: ",
       paste(with(filter(cmip6.span, !ok), paste(model, experiment)), collapse = ", "))
message(sprintf("CMIP6: %d models x %d experiments, all Jan 1850 - Dec 2099",
                n_distinct(cmip6.monthly$model), n_distinct(cmip6.monthly$experiment)))

# --- Winter means, linear detrend, rolling AR(1) ----------------------------
cmip6.ar1 <- cmip6.monthly %>%
  pivot_longer(c(EBS, GOA), names_to = "region", values_to = "sst") %>%
  filter(month %in% c(11, 12, 1, 2, 3)) %>%
  mutate(year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
  group_by(model, experiment, region, year) %>%
  filter(n() == 5) %>%                           # complete winters only
  summarise(sst = mean(sst), .groups = "drop_last") %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(ar1 = roll_ar1_right(detrend_lin(sst))) %>%
  ungroup()

# --- Skill-weighted mean across models --------------------------------------
cmip6.wts <- read.csv(file.path(cmip6.dir, "normalized_CMIP6_weights.csv")) %>%
  filter(region %in% c("Eastern_Bering_Sea", "Gulf_of_Alaska")) %>%
  mutate(region = recode(region, Eastern_Bering_Sea = "EBS", Gulf_of_Alaska = "GOA")) %>%
  select(region, model, weight = normalized_weight)

missing.wts <- setdiff(unique(cmip6.ar1$model), cmip6.wts$model)
if (length(missing.wts)) stop("No weights for: ", paste(missing.wts, collapse = ", "))

cmip6.wtd <- cmip6.ar1 %>%
  filter(!is.na(ar1)) %>%
  inner_join(cmip6.wts, by = c("region", "model")) %>%
  group_by(experiment, region, year) %>%
  summarise(ar1.wmean = weighted.mean(ar1, weight),
            ar1.wsd   = sqrt(sum(weight * (ar1 - weighted.mean(ar1, weight))^2) /
                               sum(weight)),
            n.models  = n(), .groups = "drop")

# Forced minus control, per model: mean AR(1) over each period of the forced
# run against the model's own mean control AR(1), then weighted across models
cmip6.ctrl <- cmip6.ar1 %>%
  filter(experiment == "piControl", !is.na(ar1)) %>%
  group_by(model, region) %>%
  summarise(control = mean(ar1), .groups = "drop")

cmip6.change <- bind_rows(lapply(list(c(1965, 2026), c(2012, 2026), c(2070, 2099)),
                                 function(p) {
  cmip6.ar1 %>%
    filter(experiment == "hist_ssp245", year >= p[1], year <= p[2]) %>%
    group_by(model, region) %>%
    summarise(forced = mean(ar1), .groups = "drop") %>%
    inner_join(cmip6.ctrl, by = c("model", "region")) %>%
    inner_join(cmip6.wts,  by = c("model", "region")) %>%
    group_by(region) %>%
    summarise(period          = paste(p, collapse = "-"),
              forced.wmean    = weighted.mean(forced, weight),
              control.wmean   = weighted.mean(control, weight),
              diff.wmean      = weighted.mean(forced - control, weight),
              n.models.higher = sum(forced > control),
              n.models        = n(), .groups = "drop")
}))
message("CMIP6 winter AR(1), SSP2-4.5 vs each model's own control (weighted):")
print(as.data.frame(mutate(cmip6.change, across(ends_with("wmean"), ~ round(.x, 3)))))

write.csv(cmip6.wtd, "./Output/CMIP6_winter_AR1_weighted_mean.csv", row.names = FALSE)
write.csv(cmip6.change, "./Output/CMIP6_winter_AR1_ssp245_vs_piControl.csv",
          row.names = FALSE)

# --- Figure -------------------------------------------------------------------
# Weighted model mean, with a +/- 1 weighted SD band across models
cmip6.cols <- c(piControl = "#2a78d6", hist_ssp245 = "#eb6834")
cmip6.labs <- c(piControl = "Preindustrial control", hist_ssp245 = "Historical + SSP2-4.5")

p.cmip6 <- cmip6.wtd %>%
  mutate(region     = factor(region, levels = c("GOA", "EBS"),
                             labels = c("Gulf of Alaska", "Eastern Bering Sea")),
         experiment = factor(experiment, levels = names(cmip6.cols))) %>%
  ggplot(aes(x = year, y = ar1.wmean, colour = experiment, fill = experiment)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_ribbon(aes(ymin = ar1.wmean - ar1.wsd, ymax = ar1.wmean + ar1.wsd),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~ region, nrow = 1) +
  scale_colour_manual(values = cmip6.cols, labels = cmip6.labs, name = NULL) +
  scale_fill_manual(values = cmip6.cols, labels = cmip6.labs, name = NULL) +
  labs(x = "Year (end of 15-yr window)",
       y = "Winter SST AR(1), weighted model mean") +
  theme_grl() +
  theme(legend.position = "bottom")

ggsave("./Figures/CMIP6_winter_AR1_ssp245_piControl.png", plot = p.cmip6,
       width = grl.width.full, height = 3.0, dpi = 400, bg = "white")
message("Saved: Figures/CMIP6_winter_AR1_ssp245_piControl.png")

# ============================================================
# SECTION 24: MANUSCRIPT FIGURE 3 — LME DFA (a-b) + CMIP6 AR(1) (c-d)
# ============================================================
# Top row, Section 22: (a) boreal-LME loadings on the shared winter AR(1)
# trend, (b) that trend with the GOA/EBS AR(1) laid over it. Bottom row,
# Section 23: CMIP6 winter AR(1), SSP2-4.5 vs preindustrial control, for
# (c) GOA and (d) EBS on a common y axis, with the observed (ERA5, Fig. 2)
# AR(1) for each region added. Each row keeps its own legend, so the rows
# are assembled separately and then stacked.

fig3.legend.theme <- theme(legend.position    = "bottom",
                           legend.margin      = margin(0, 0, 0, 0),
                           legend.box.spacing = unit(1, "mm"))

fig3.row1 <- (dfa_loadings_panel(dfa.all$loadings, "(a)") |
                dfa_overlay_panel(dfa.all, "(b)")) +
  plot_layout(widths = c(1, 1.15), guides = "collect") &
  fig3.legend.theme

fig3.obs <- filter(study.ar1, !is.na(ar1))

fig3.cmip6.ylim <- range(c(cmip6.wtd$ar1.wmean - cmip6.wtd$ar1.wsd,
                           cmip6.wtd$ar1.wmean + cmip6.wtd$ar1.wsd, fig3.obs$ar1))

# Observed series join the colour legend; the fill scale carries the same
# levels (observed unfilled) so the two legends merge into one
fig3.levels <- c(names(cmip6.cols), "observed")
fig3.cols   <- c(cmip6.cols, observed = "black")
fig3.labs   <- c(cmip6.labs, observed = "Observed (ERA5)")

# One region per panel, so each carries its own tag; the region name sits in
# a strip as in the standalone CMIP6 figure. `observed = FALSE` leaves out the
# ERA5 series (and its legend key) on the same axes, for the first step of
# the slide build in Section 27.
cmip6_panel <- function(reg, reg.label, tag, y.title = TRUE, observed = TRUE) {
  levs <- if (observed) fig3.levels else names(cmip6.cols)
  cmip6.wtd %>%
    filter(region == reg) %>%
    mutate(experiment = factor(experiment, levels = names(cmip6.cols)),
           panel      = reg.label) %>%
    ggplot(aes(x = year, y = ar1.wmean, colour = experiment, fill = experiment)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
    geom_ribbon(aes(ymin = ar1.wmean - ar1.wsd, ymax = ar1.wmean + ar1.wsd),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.4) +
    { if (observed)
        geom_line(data = filter(fig3.obs, region == reg),
                  aes(x = year, y = ar1, colour = "observed"),
                  inherit.aes = FALSE, linewidth = 0.5) } +
    facet_wrap(~ panel) +
    scale_colour_manual(values = fig3.cols, limits = levs,
                        labels = fig3.labs[levs], name = NULL) +
    scale_fill_manual(values = c(cmip6.cols, observed = NA), limits = levs,
                      labels = fig3.labs[levs], name = NULL) +
    coord_cartesian(ylim = fig3.cmip6.ylim) +
    labs(x = "Year (end of 15-yr window)",
         y = if (y.title) "Winter SST AR(1), 15-yr window" else NULL,
         tag = tag) +
    theme_grl() + lme.panel.theme
}

fig3.row2 <- (cmip6_panel("GOA", "Gulf of Alaska", "(c)") |
                cmip6_panel("EBS", "Eastern Bering Sea", "(d)", y.title = FALSE)) +
  plot_layout(guides = "collect") &
  fig3.legend.theme

p.fig3 <- wrap_elements(full = fig3.row1) / wrap_elements(full = fig3.row2) +
  plot_layout(heights = c(1.15, 1))

ggsave("./Figures/Fig_3.png", plot = p.fig3,
       width = grl.width.full, height = 6.2, dpi = 400, bg = "white")
message("Saved: Figures/Fig_3.png")

# ============================================================
# SECTION 25: MANUSCRIPT FIGURE 4 — SST AR(1) and AL SLP SD
# ============================================================
# (a) Section 17's cellwise regression of SST AR(1) on AL SLP SD. (b) GOA
# and (c) EBS: 15-yr rolling AL SLP SD and SST AR(1) on dual axes, as in
# Section 5 but from complete winters only, like Figs. 2-3. Both panels use
# one AR(1) range and one AL SD range, so the axes match across the row and
# each axis title is drawn once: AR(1) at left, AL SD at right, year beneath.
# The figure is built at GRL's full page width (190 mm), so it prints with
# no scaling and 8-pt text stays 8 pt. The map row is taller than the
# time-series row; the map keeps its projection's shape and is centred.
# Saved as PNG and as a PDF at the same physical size, for checking font
# sizes at 100% zoom.

# Section 17's map, redrawn with 6-pt area labels: at this panel size the
# default 7-pt labels run over the EBS polygon and the AL box
p.fig4a <- ps_regression_map(sst.ar1.reg.grid,
                             fill.name  = expression(beta~"(AR(1) per SD)"),
                             col.lim    = col.lim.sst,
                             label.size = 6 / .pt) +
  labs(tag = "(a)") +
  theme(plot.tag = element_text(size = 9, face = "bold"))

# --- Complete winters only (all five Nov-Mar months) ------------------------
# SST AR(1) is Fig. 2's series (Section 21). The AL series drops the partial
# winter of 1950 before the rolling SD, so both series start in 1965.
fig4.al <- al.monthly %>%
  mutate(year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
  filter(month %in% c(11, 12, 1, 2, 3)) %>%
  group_by(year) %>%
  filter(n() == 5) %>%
  summarise(SLP = mean(SLP), .groups = "drop") %>%
  arrange(year) %>%
  mutate(AL.sd = as.numeric(scale(roll_sd(SLP))))

fig4.dat <- fig2.roll %>%
  pivot_wider(names_from = region, values_from = ar1, names_glue = "{region}.ar1") %>%
  inner_join(select(fig4.al, year, AL.sd), by = "year") %>%
  filter(!is.na(GOA.ar1), !is.na(EBS.ar1), !is.na(AL.sd))

fig4.chel.goa <- chelton_test(fig4.dat$GOA.ar1, fig4.dat$AL.sd)
fig4.chel.ebs <- chelton_test(fig4.dat$EBS.ar1, fig4.dat$AL.sd)
message(sprintf("Fig. 4, complete winters %d-%d: GOA r = %.2f, p = %.3f; EBS r = %.2f, p = %.3f",
                min(fig4.dat$year), max(fig4.dat$year),
                fig4.chel.goa$r.obs, fig4.chel.goa$pval,
                fig4.chel.ebs$r.obs, fig4.chel.ebs$pval))

fig4.al.sd.range <- range(fig4.dat$AL.sd)
fig4.ar1.range   <- range(c(fig4.dat$GOA.ar1, fig4.dat$EBS.ar1))
fig4.k <- diff(fig4.ar1.range) / diff(fig4.al.sd.range)   # AL SD onto the AR(1) axis
fig4.b <- fig4.ar1.range[1] - fig4.k * fig4.al.sd.range[1]

# --- Shared by Figs. 4-5 -----------------------------------------------------
# One dual-axis panel: a region's 15-yr AR(1) on the left axis and AL SLP SD
# on the right, mapped onto the AR(1) axis by k and b. `d` holds year,
# GOA.ar1, EBS.ar1 and AL.sd; `series` names the GOA, EBS and AL series.
# Both panels of a figure carry all three series in their scales, so their
# legends are identical and collect into one.
dual_ar1_panel <- function(d, reg, reg.label, chel, tag, series, k, b, ar1.range,
                           left.title = NULL, right.title = NULL) {
  cols   <- setNames(c(area.colors[["GOA"]], area.colors[["EBS"]], "grey30"), series)
  ltys   <- setNames(c("solid", "solid", "dashed"), series)
  shapes <- setNames(c(16, 16, 1), series)
  d.long <- bind_rows(
    data.frame(year = d$year, y = d[[paste0(reg, ".ar1")]], series = series[[reg]]),
    data.frame(year = d$year, y = k * d$AL.sd + b, series = series[["AL"]])) %>%
    mutate(panel = reg.label)
  lbl <- sprintf("r = %.2f, p = %.3f\nN* = %.1f", chel$r.obs, chel$pval, chel$n.eff)

  # show.legend = TRUE: ggplot otherwise keys only the series present in
  # each panel, and the two legends would then differ and not collect
  ggplot(d.long, aes(x = year, y = y, colour = series, linetype = series,
                     shape = series)) +
    geom_line(linewidth = 0.4, show.legend = TRUE) +
    geom_point(size = 0.8, show.legend = TRUE) +
    annotate("text", x = min(d.long$year), y = -Inf, label = lbl, hjust = 0,
             vjust = -0.35, size = 7 / .pt, lineheight = 1) +
    facet_wrap(~ panel) +
    scale_colour_manual(values = cols, limits = unname(series), name = NULL) +
    scale_linetype_manual(values = ltys, limits = unname(series), name = NULL) +
    scale_shape_manual(values = shapes, limits = unname(series), name = NULL) +
    scale_y_continuous(name = left.title, limits = ar1.range,
                       sec.axis = sec_axis(~ (. - b) / k, name = right.title)) +
    labs(x = "Year (end of 15-yr window)", tag = tag) +
    theme_grl() +
    theme(plot.tag = element_text(size = 9, face = "bold"))
}

# Assemble and save a three-panel figure at GRL's full page width (190 mm),
# so it prints unscaled: (a) map on top, centred; (b)-(c) below with one
# legend and shared axis titles. Panel sizes are in inches: (b)-(c) are as
# wide as two dual-axis panels allow within 190 mm and 0.71 as tall as
# wide; the map is taller and keeps its projection's aspect. Saved as PNG
# and as a PDF at the same physical size, for checking fonts at 100% zoom.
save_map_dual_figure <- function(p.map, p.b, p.c, file.stem,
                                 map.h = 2.6, dual.w = 2.75, dual.h = 1.95) {
  map.w <- map.h * diff(p.map$coordinates$limits$x) /
                   diff(p.map$coordinates$limits$y)
  row1 <- p.map +
    plot_layout(widths = unit(map.w, "in"), heights = unit(map.h, "in"))
  row2 <- (p.b | p.c) +
    plot_layout(widths  = unit(rep(dual.w, 2), "in"),
                heights = unit(dual.h, "in"),
                guides = "collect", axis_titles = "collect_x") &
    theme(legend.position    = "bottom",
          legend.margin      = margin(0, 0, 0, 0),
          legend.box.spacing = unit(1, "mm"))
  # Each row: its panel plus the space its labels, strip, colour bar and
  # legend need around it
  row.h <- c(map.h + 0.80, dual.h + 1.00)
  fig <- wrap_elements(full = row1) / wrap_elements(full = row2) +
    plot_layout(heights = row.h)
  # Standard pdf() device: cairo_pdf needs XQuartz, which this machine lacks
  for (ext in c("png", "pdf"))
    ggsave(paste0(file.stem, ".", ext), plot = fig, width = grl.width.full,
           height = sum(row.h), dpi = 400, bg = "white")
  message(sprintf("Saved: %s.png and .pdf (%.0f x %.0f mm)", file.stem,
                  grl.width.full * 25.4, sum(row.h) * 25.4))
  invisible(fig)
}

fig4.series <- c(GOA = "Gulf of Alaska SST AR(1)",
                 EBS = "Eastern Bering Sea SST AR(1)",
                 AL  = "AL SLP SD")

p.fig4 <- save_map_dual_figure(
  p.fig4a,
  dual_ar1_panel(fig4.dat, "GOA", "Gulf of Alaska", fig4.chel.goa, "(b)",
                 fig4.series, fig4.k, fig4.b, fig4.ar1.range,
                 left.title = "SST AR(1)"),
  dual_ar1_panel(fig4.dat, "EBS", "Eastern Bering Sea", fig4.chel.ebs, "(c)",
                 fig4.series, fig4.k, fig4.b, fig4.ar1.range,
                 right.title = "AL SLP SD (z-scored)"),
  "./Figures/Fig_4")

# ============================================================
# SECTION 26: MANUSCRIPT FIGURE 5 — MLD and AL SLP
# ============================================================
# Same layout and size as Fig. 4. (a) Section 6's cellwise regression of
# winter MLD anomaly on winter AL SLP anomaly. (b) GOA and (c) EBS: 15-yr
# rolling MLD AR(1), from Section 7's detrended winter MLD anomalies, with
# the complete-winter AL SLP SD used in Fig. 4.
# Complete winters only: ORAS5 starts in Jan 1958, so winter 1958 (Jan-Mar)
# is dropped and MLD AR(1) starts in 1973. Winter 2026 is kept: the ORAS5
# update in downloads.R Section 4 requested Jan-Mar 2026, but the cached
# winter series cannot confirm that March is present.

p.fig5a <- ps_regression_map(mld.reg.grid,
                             fill.name  = expression(beta~"(m / hPa)"),
                             col.lim    = col.lim.mld,
                             label.size = 6 / .pt) +
  labs(tag = "(a)") +
  theme(plot.tag = element_text(size = 9, face = "bold"))

fig5.dat <- mld.win %>%
  filter(year > 1958) %>%
  arrange(year) %>%
  mutate(GOA.ar1 = roll_ar1_right(GOA),
         EBS.ar1 = roll_ar1_right(EBS)) %>%
  select(year, GOA.ar1, EBS.ar1) %>%
  inner_join(select(fig4.al, year, AL.sd), by = "year") %>%
  filter(!is.na(GOA.ar1), !is.na(EBS.ar1), !is.na(AL.sd))

fig5.chel.goa <- chelton_test(fig5.dat$GOA.ar1, fig5.dat$AL.sd)
fig5.chel.ebs <- chelton_test(fig5.dat$EBS.ar1, fig5.dat$AL.sd)
message(sprintf("Fig. 5, complete winters %d-%d: GOA r = %.2f, p = %.3f; EBS r = %.2f, p = %.3f",
                min(fig5.dat$year), max(fig5.dat$year),
                fig5.chel.goa$r.obs, fig5.chel.goa$pval,
                fig5.chel.ebs$r.obs, fig5.chel.ebs$pval))

fig5.ar1.range <- range(c(fig5.dat$GOA.ar1, fig5.dat$EBS.ar1))
fig5.k <- diff(fig5.ar1.range) / diff(range(fig5.dat$AL.sd))
fig5.b <- fig5.ar1.range[1] - fig5.k * min(fig5.dat$AL.sd)
fig5.series <- c(GOA = "Gulf of Alaska MLD AR(1)",
                 EBS = "Eastern Bering Sea MLD AR(1)",
                 AL  = "AL SLP SD")

p.fig5 <- save_map_dual_figure(
  p.fig5a,
  dual_ar1_panel(fig5.dat, "GOA", "Gulf of Alaska", fig5.chel.goa, "(b)",
                 fig5.series, fig5.k, fig5.b, fig5.ar1.range,
                 left.title = "MLD AR(1)"),
  dual_ar1_panel(fig5.dat, "EBS", "Eastern Bering Sea", fig5.chel.ebs, "(c)",
                 fig5.series, fig5.k, fig5.b, fig5.ar1.range,
                 right.title = "AL SLP SD (z-scored)"),
  "./Figures/Fig_5")

# ============================================================
# SECTION 27: PRESENTATION FIGURES — slide versions of Figs. 1-5
# ============================================================
# Versions of the manuscript figures for a PowerPoint talk, written to
# ./presentation_figs. Each is drawn at the size it will be shown on a 16:9
# slide (13.33 x 7.5 in): 5.75 in tall, which fills the slide below a title,
# and at most 12.3 in wide. Inserted at 100%, the text is the point size it
# shows at -- 18-pt axis titles, 16-pt tick labels -- legible from the back
# of a room. Compared with the manuscript versions:
#   - Multi-panel figures are split so each slide carries one result, and
#     panel tags are dropped, as there is no caption for them to refer to.
#   - Strokes are heavier and the CMIP6 bands darker, to survive a projector;
#     facet strips lose their grey boxes and read as panel titles.
#   - Fig. 3 opens with a map of the DFA's LMEs, and its CMIP6 panels come
#     as a two-step build: models first, then the observed series added on
#     identical axes, so the step doesn't jump.
# The panels come from the manuscript builders (Sections 21-26), restyled,
# so data, colours and statistics stay identical to the manuscript.
# PNG at 300 dpi, which PowerPoint inserts at its true size.
# Requires Sections 1-7, 17 and 21-26.

talk.dir <- "./presentation_figs"
dir.create(talk.dir, showWarnings = FALSE)

talk.w <- 12.3    # in: slide width less margins
talk.h <- 5.75    # in: slide height less a title
talk.regions <- c(GOA = "Gulf of Alaska", EBS = "Eastern Bering Sea")

theme_talk <- function(base_size = 18) {
  theme_grl(base_size) +
    theme(axis.text          = element_text(size = base_size - 2, colour = "black"),
          axis.ticks         = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks.length  = unit(2, "mm"),
          legend.text        = element_text(size = base_size - 2),
          legend.key.width   = unit(12, "mm"),
          legend.position    = "bottom",
          legend.margin      = margin(0, 0, 0, 0),
          legend.box.spacing = unit(2, "mm"),
          panel.grid.major   = element_line(colour = "grey90", linewidth = 0.4),
          panel.border       = element_rect(colour = "black", linewidth = 0.8, fill = NA),
          plot.title         = element_text(size = base_size + 2, hjust = 0.5),
          strip.background   = element_blank(),
          strip.text         = element_text(size = base_size + 2,
                                            margin = margin(0, 0, 2, 0, "mm")),
          plot.margin        = margin(2, 4, 2, 4, "mm"))
}

# Restyle a manuscript panel for a slide: thicker strokes and points,
# annotation text scaled from the manuscript's 7 pt to 16 pt, darker bands,
# no panel tag and the slide theme. Layers are changed in place, so pass a
# panel built for the slide, never one a manuscript figure still uses.
talkify <- function(p, line.k = 2.5, text.k = 16 / 7, alpha.k = 1.7) {
  for (l in p$layers) {
    ap <- l$aes_params
    if (inherits(l$geom, "GeomText")) {
      if (!is.null(ap$size)) ap$size <- ap$size * text.k
    } else {
      for (a in intersect(c("linewidth", "size"), names(ap))) ap[[a]] <- ap[[a]] * line.k
      if (!is.null(ap$alpha)) ap$alpha <- min(1, ap$alpha * alpha.k)
    }
    l$aes_params <- ap
  }
  p + labs(tag = NULL) + theme_talk()
}

talk_save <- function(p, name, width = talk.w, height = talk.h) {
  f <- file.path(talk.dir, paste0(name, ".png"))
  ggsave(f, plot = p, width = width, height = height, dpi = 300, bg = "white")
  message(sprintf("Saved: %s (%.2f x %.2f in)", f, width, height))
}

# --- Maps: Figs. 1, 4a, 5a ---------------------------------------------------
# ps_regression_map() at slide scale, with the colour bar upright at the
# right so the map takes the full height. The figure width is solved from
# the map's fixed aspect and the measured width of everything else (labels,
# colour bar), so the saved figure has no slack beside the map.
talk_map <- function(grid, fill.name, col.lim, ...) {
  ps_regression_map(grid, fill.name = fill.name, col.lim = col.lim,
                    label.size = 14 / .pt, grat.size = 15 / .pt,
                    line.scale = 2.4, ...) +
    guides(fill = guide_colourbar(theme = theme(
      legend.title.position = "top",
      legend.key.width      = unit(6, "mm"),
      legend.key.height     = unit(90, "mm"),
      legend.frame          = element_rect(colour = "black", linewidth = 0.5),
      legend.ticks          = element_line(colour = "black", linewidth = 0.5),
      legend.ticks.length   = unit(1.5, "mm")))) +
    theme(legend.position  = "right",
          legend.direction = "vertical",
          legend.title     = element_text(size = 16),
          legend.text      = element_text(size = 16),
          legend.box.spacing = unit(4, "mm"),
          plot.margin      = margin(2, 2, 2, 2, "mm"))
}

talk_map_width <- function(p, height = talk.h) {
  gt   <- ggplotGrob(p)
  fixd <- function(u, conv) sum(conv(u[grid::unitType(u) != "null"], "in",
                                     valueOnly = TRUE))
  asp  <- diff(p$coordinates$limits$x) / diff(p$coordinates$limits$y)
  fixd(gt$widths, grid::convertWidth) +
    (height - fixd(gt$heights, grid::convertHeight)) * asp
}

save_talk_map <- function(p, name) talk_save(p, name, width = talk_map_width(p))

# Fig. 1: colours as in the manuscript (RdBu stretched over the data range),
# so blue marks weak positive loadings, not negative ones
save_talk_map(talk_map(fig1.eof, "SLP EOF1\nloading", max(abs(fig1.eof$beta)),
                       land.fill  = NA,
                       fill.scale = scale_fill_distiller(
                         palette = "RdBu", direction = -1,
                         name = "SLP EOF1\nloading")),
              "Fig1_study_site_map")

save_talk_map(talk_map(sst.ar1.reg.grid, "β (AR(1)\nper SD)", col.lim.sst),
              "Fig4_1_SST_AR1_regression_map")

save_talk_map(talk_map(mld.reg.grid, "β\n(m / hPa)", col.lim.mld),
              "Fig5_1_MLD_regression_map")

# --- Fig. 2: winter SST, then its rolling AR(1) -------------------------------
# One variable per slide, GOA and EBS side by side, or stacked with GOA on
# top (`stack = TRUE`). Each y axis is scaled to its own panel as in the
# manuscript; stacked panels share the year axis, drawn once at the bottom.
talk_region_pair <- function(panel_fn, ..., stack = FALSE) {
  p.goa <- talkify(panel_fn("GOA", NULL) + labs(title = talk.regions[["GOA"]], ...))
  p.ebs <- talkify(panel_fn("EBS", NULL) + labs(title = talk.regions[["EBS"]], ...))
  if (stack)
    (p.goa / p.ebs) + plot_layout(axis_titles = "collect", axes = "collect_x")
  else
    (p.goa | p.ebs) + plot_layout(axis_titles = "collect")
}

# Stacked, the SST slide is taller and narrower than the others; 6.5 in
# still fits a 7.5-in slide under a one-line title
talk_save(talk_region_pair(ts_panel, x = "Year", y = "Winter SST (°C)",
                           stack = TRUE),
          "Fig2_1_winter_SST", width = 9, height = 6.5)
talk_save(talk_region_pair(ar1_panel, x = "Year (end of 15-yr window)",
                           y = "Winter SST AR(1), 15-yr window"),
          "Fig2_2_SST_AR1")

# --- Fig. 3 context: map of the boreal LMEs ---------------------------------
# The LMEs of the Section 22 DFA, outlined on the polar stereographic
# projection of Figs. 1, 4 and 5 (Pacific at the bottom), with the EBS and GOA
# study polygons shaded in their figure colours. Land is drawn over the
# shading, so only the polygons' ocean is coloured. Labels sit in two columns
# beside the map, each joined to its region by a leader line: the European
# LMEs are too small to carry their own.
lme.map.lat <- 38       # southern edge of the map; the LMEs reach 40.8N
lme.map.R   <- -ps_project(0, lme.map.lat, lon0 = 0)$y  # its radius
lme.dbf   <- foreign::read.dbf("./Data/LME shapefiles/LMEs66.dbf", as.is = TRUE)
lme.rings <- read_shp_rings("./Data/LME shapefiles/LMEs66.shp")
lme.keep  <- which(lme.dbf$LME_NUMBER %in% lme.numbers &
                   !lme.dbf$LME_NAME %in% lme.drop)

# Outlines as projected paths, broken wherever an edge runs along the 180
# meridian: the shapefile cuts the Northern Bering - Chukchi Seas there, and
# drawing the cut would put a false boundary through it
lme_outline <- function(k) {
  bind_rows(lapply(split(lme.rings[[k]], lme.rings[[k]]$ring), function(g) {
    seam <- abs(g$x) > 179.999
    runs <- split(g, cumsum(c(FALSE, seam[-1] & seam[-nrow(g)])))
    bind_rows(lapply(names(runs), function(j) {
      if (nrow(runs[[j]]) < 2) return(NULL)
      d <- densify_lonlat(runs[[j]]$x, runs[[j]]$y)
      cbind(ps_project(d$lon, d$lat), grp = paste(k, g$ring[1], j))
    }))
  }))
}

# Label anchor: the in-LME point of a coarse grid nearest the median of all
# such points, so it falls inside even a concave or multi-part LME
lme_anchor <- function(k) {
  r <- lme.rings[[k]]
  g <- expand.grid(x = seq(min(r$x), max(r$x), length.out = 60),
                   y = seq(min(r$y), max(r$y), length.out = 60))
  g <- g[in_rings(g$x, g$y, r), ]
  g[which.min((g$x - median(g$x))^2 + (g$y - median(g$y))^2), ]
}

lme.lines <- bind_rows(lapply(lme.keep, lme_outline))
lme.labs  <- bind_rows(
  lapply(lme.keep, function(k)
    cbind(lme_anchor(k), lab = lme.dbf$LME_NAME[k], face = "plain")),
  # Study polygons: fixed points in their open ocean. GOA's is at the
  # polygon's western end, so its leader passes below the EBS polygon.
  data.frame(x = c(186, 200), y = c(58.5, 52),
             lab = unname(talk.regions[c("EBS", "GOA")]), face = "bold"))
lme.labs[c("x", "y")] <- ps_project(lme.labs$x, lme.labs$y)
lme.labs$lab[lme.labs$lab == "Northern Bering - Chukchi Seas"] <-
  "Northern Bering -\nChukchi Seas"
# Left column: the anchors left of the pole, plus the Baltic (just right of
# the pole) to even the columns up, and GOA so the two study systems sit
# together, GOA under EBS
lme.labs$left <- lme.labs$x < 0 |
  lme.labs$lab %in% c("Baltic Sea", talk.regions[["GOA"]])

lme.study <- bind_rows(lapply(c("EBS", "GOA"), function(a) {
  poly <- if (a == "EBS") ebs.poly else goa.poly
  cbind(ps_path(c(poly$lon, poly$lon[1]), c(poly$lat, poly$lat[1]), n = 500),
        area = a)
}))
lme.coast  <- ps_coast(mapWorld.clean, c(0, 360), c(lme.map.lat, 90))
lme.edge   <- ps_path(c(0, 360), rep(lme.map.lat, 2), n = 1000)
lme.grat   <- bind_rows(
  bind_rows(lapply(c(40, 60, 80), function(P)
    cbind(ps_path(c(0, 360), c(P, P), n = 600), grp = paste0("p", P)))),
  bind_rows(lapply(seq(0, 330, 30), function(L)
    cbind(ps_path(c(L, L), c(lme.map.lat, 80)), grp = paste0("m", L)))))

# Panel limits: the disc plus a label column each side, both as wide as the
# longest label, converted from inches at the scale the 5.75-in height sets.
# Equal columns put the pole at the centre of the figure.
lme.lab.size <- 14
lab_width <- function(lab, face)
  max(mapply(function(l, f) grid::convertWidth(grid::grobWidth(grid::textGrob(
    l, gp = grid::gpar(fontsize = lme.lab.size, fontface = f))), "in",
    valueOnly = TRUE), lab, face))
lme.seg.gap <- 0.06                                # disc to leader-line end
lme.txt.pad <- 0.012                               # leader-line end to text
lme.gap     <- lme.seg.gap + lme.txt.pad
lme.ylim    <- c(-1, 1) * (lme.map.R + 0.02)
lme.in.per.unit <- (talk.h - 0.16) / diff(lme.ylim) # less the 2-mm margins
lme.xlim  <- c(-1, 1) * (lme.map.R + lme.gap + 0.02 +
                         lab_width(lme.labs$lab, lme.labs$face) / lme.in.per.unit)

# Each column keeps its labels in the order of their anchors, pushing
# overlapping neighbours apart evenly until none overlap, inside the panel.
# A leader line runs straight from each anchor to its label.
lme.line.h <- lme.lab.size / 72 * 1.2 / lme.in.per.unit   # one text line
lme_label_column <- function(d, left) {
  d  <- d[order(-d$y), ]
  n  <- nrow(d)
  h  <- lme.line.h * (1 + lengths(regmatches(d$lab, gregexpr("\n", d$lab))))
  ly <- d$y
  for (it in 1:1000) {
    need <- (h[-n] + h[-1]) / 2 - (ly[-n] - ly[-1])   # overlap of each pair
    if (all(need <= 1e-9)) break
    i <- which(need > 1e-9)
    ly[i]     <- ly[i]     + need[i] / 2
    ly[i + 1] <- ly[i + 1] - need[i] / 2
    ly <- pmin(pmax(ly, lme.ylim[1] + h / 2), lme.ylim[2] - h / 2)
  }
  side <- if (left) -1 else 1
  d$ly <- ly
  d$seg.x <- side * (lme.map.R + lme.seg.gap)
  d$txt.x <- side * (lme.map.R + lme.gap)
  list(geom_segment(data = d, aes(x = x, y = y, xend = seg.x, yend = ly),
                    colour = "grey35", linewidth = 0.4),
       geom_text(data = d, aes(x = txt.x, y = ly, label = lab, fontface = face),
                 hjust = if (left) 1 else 0, size = lme.lab.size / .pt,
                 lineheight = 0.9))
}

p.talk.lme.map <- ggplot() +
  geom_polygon(data = lme.study, aes(x = x, y = y, group = area, fill = area),
               alpha = 0.6) +
  # Land outlined in its own fill colour: world2 cuts land at 0/360, and a
  # contrasting outline would draw that cut through Britain and France
  geom_polygon(data = lme.coast, aes(x = x, y = y, group = group),
               fill = "grey85", colour = "grey85", linewidth = 0.2) +
  geom_path(data = lme.grat, aes(x = x, y = y, group = grp),
            colour = "grey55", linewidth = 0.3, linetype = "dotted") +
  geom_path(data = lme.lines, aes(x = x, y = y, group = grp),
            colour = "grey20", linewidth = 0.45) +
  geom_path(data = lme.study, aes(x = x, y = y, group = area, colour = area),
            linewidth = 1) +
  geom_path(data = lme.edge, aes(x = x, y = y), colour = "black", linewidth = 0.8) +
  geom_point(data = lme.labs, aes(x = x, y = y), size = 1.3, colour = "grey15") +
  lme_label_column(filter(lme.labs, left), TRUE) +
  lme_label_column(filter(lme.labs, !left), FALSE) +
  scale_fill_manual(values = area.colors, guide = "none") +
  scale_colour_manual(values = area.colors, guide = "none") +
  coord_equal(xlim = lme.xlim, ylim = lme.ylim, expand = FALSE) +
  theme_grl_map() +
  theme(plot.margin = margin(2, 2, 2, 2, "mm"))

talk_save(p.talk.lme.map, "Fig3_1_LME_map",
          width = diff(lme.xlim) * lme.in.per.unit + 0.16)

# --- Fig. 3: LME DFA, then CMIP6 vs observed --------------------------------
# Collected legends take their position from the patchwork, not the panels,
# hence the `& theme()` on each multi-panel slide.
talk.collect <- theme(legend.position = "bottom")

# The trend legend sits in its own cell under the trend panel, one entry per
# row (each carries its r and p, too long to share a row at 16 pt), while
# the loadings panel spans the full height beside both. The trend panel is
# then too short for its y title on one line.
p.talk.lme <- talkify(dfa_loadings_panel(dfa.all$loadings, NULL)) +
  (talkify(dfa_overlay_panel(dfa.all, NULL)) +
     labs(y = "Winter SST AR(1),\nstandardized") +
     guides(colour = guide_legend(ncol = 1))) +
  guide_area() +
  plot_layout(design = "AB\nAC", widths = c(1, 1.5),
              heights = unit(c(1, 0.95), c("null", "in")), guides = "collect")
talk_save(p.talk.lme, "Fig3_2_LME_DFA")

talk_cmip6 <- function(observed)
  (talkify(cmip6_panel("GOA", talk.regions[["GOA"]], NULL, observed = observed)) |
     talkify(cmip6_panel("EBS", talk.regions[["EBS"]], NULL, y.title = FALSE,
                         observed = observed))) +
  plot_layout(guides = "collect", axis_titles = "collect_x") & talk.collect

talk_save(talk_cmip6(observed = FALSE), "Fig3_3_CMIP6_models")
talk_save(talk_cmip6(observed = TRUE),  "Fig3_4_CMIP6_models_observed")

# --- Figs. 4-5: AR(1) vs AL SLP SD, dual axes --------------------------------
talk_dual <- function(d, chel.goa, chel.ebs, series, k, b, ar1.range, left.title)
  (talkify(dual_ar1_panel(d, "GOA", talk.regions[["GOA"]], chel.goa, NULL,
                          series, k, b, ar1.range, left.title = left.title)) |
     talkify(dual_ar1_panel(d, "EBS", talk.regions[["EBS"]], chel.ebs, NULL,
                            series, k, b, ar1.range,
                            right.title = "AL SLP SD (z-scored)"))) +
  plot_layout(guides = "collect", axis_titles = "collect_x") & talk.collect

talk_save(talk_dual(fig4.dat, fig4.chel.goa, fig4.chel.ebs, fig4.series,
                    fig4.k, fig4.b, fig4.ar1.range, "SST AR(1)"),
          "Fig4_2_SST_AR1_AL_SD")
talk_save(talk_dual(fig5.dat, fig5.chel.goa, fig5.chel.ebs, fig5.series,
                    fig5.k, fig5.b, fig5.ar1.range, "MLD AR(1)"),
          "Fig5_2_MLD_AR1_AL_SD")

# --- Context: Alaska fishery disasters since 2014, declared and pending -------
# Federal fishery resource disaster determinations for Alaska, from NOAA's
# Fishery Resource Disaster Online System (disasters.fisheries.noaa.gov,
# read 2026-10-24; all 38 Alaska requests are saved in
# Data/noaa_fishery_disasters_alaska.csv). Mapped: the 22 approved
# (declared) disasters, all for fishing years 2016-2023 -- none was declared
# for 2014 or 2015 -- and the 14 requests still pending (2023-2025), which
# are labelled the same way with their seasons starred. Not mapped: the 2
# not approved. A request can cover several fisheries; each is placed at its
# fishery area, and areas share one label listing their fisheries and
# seasons. The GOA Pacific cod and Bering Sea crab markers stand for
# fisheries spread across those seas.
ak.dis.noaa <- read.csv("./Data/noaa_fishery_disasters_alaska.csv")

# number = NOAA disaster number; area keys index ak.dis.areas below
ak.dis <- data.frame(rbind(
  c(72,  "kod",   "pink salmon",     "2016"),
  c(72,  "pws",   "pink salmon",     "2016"),
  c(72,  "chig",  "pink salmon",     "2016"),
  c(72,  "lci",   "pink salmon",     "2016"),
  c(72,  "sap",   "pink salmon",     "2016"),
  c(72,  "sea",   "pink salmon",     "2016"),
  c(72,  "yak",   "pink salmon",     "2016"),
  c(81,  "cod",   "Pacific cod",     "2018"),
  c(85,  "chig",  "sockeye salmon",  "2018"),
  c(96,  "nort",  "red king crab",   "2019"),
  c(106, "uci",   "salmon",          "2018"),   # East Side Set Net
  c(106, "uci",   "salmon",          "2020"),
  c(107, "pws",   "salmon",          "2018"),   # Copper River
  c(107, "pws",   "salmon",          "2020"),
  c(108, "tann",  "Tanner crab",     "2019/20"),
  c(109, "cod",   "Pacific cod",     "2020"),
  c(110, "nort",  "salmon",          "2020"),
  c(110, "yuk",   "salmon",          "2020"),
  c(110, "chig",  "salmon",          "2020"),
  c(110, "kusk",  "salmon",          "2020"),
  c(110, "sea",   "salmon",          "2020"),
  c(110, "yuk",   "salmon",          "2021"),
  c(115, "chig",  "salmon",          "2021"),
  c(116, "bbrkc", "red king crab",   "2021/22"),
  c(116, "snow",  "snow crab",       "2021/22"),
  c(118, "nort",  "salmon",          "2021"),
  c(119, "nort",  "red king crab",   "2020"),
  c(119, "nort",  "red king crab",   "2021"),
  c(121, "kusk",  "salmon",          "2021"),
  c(123, "pws",   "salmon",          "2020"),   # Copper River/PWS coho and pink
  c(126, "bbrkc", "red king crab",   "2022/23"),
  c(126, "snow",  "snow crab",       "2022/23"),
  c(127, "yuk",   "salmon",          "2022"),
  c(136, "kusk",  "salmon",          "2022"),
  c(138, "uci",   "salmon",          "2021"),   # East Side Setnet
  c(138, "uci",   "salmon",          "2022"),
  c(138, "nel",   "salmon",          "2022"),
  c(140, "snow",  "snow crab",       "2023/24"),
  c(141, "chig",  "salmon",          "2022"),
  c(142, "uci",   "salmon",          "2023")),  # East Side Setnet
  stringsAsFactors = FALSE)
names(ak.dis) <- c("number", "area", "fishery", "season")
ak.dis$status <- "declared"

# Pending requests. The Alaska Peninsula sockeye requests (158, 169) do not
# split the peninsula, so like 72 they are listed under both its areas; the
# North Peninsula area also holds the declared 2022 Nelson Lagoon disaster.
ak.dis.pending <- data.frame(rbind(
  c(145, "nel",   "salmon",          "2023"),   # North Peninsula drift gillnet
  c(151, "kotz",  "chum salmon",     "2024"),
  c(154, "pws",   "salmon",          "2024"),   # pink and chum
  c(155, "lci",   "pink salmon",     "2024"),
  c(156, "kod",   "pink salmon",     "2024"),
  c(157, "chig",  "salmon",          "2024"),
  c(158, "nel",   "salmon",          "2024"),   # Alaska Peninsula sockeye
  c(158, "sap",   "salmon",          "2024"),   # ... and South Alaska pink
  c(164, "uci",   "salmon",          "2024"),   # East Side Setnet
  c(166, "yuk",   "salmon",          "2024"),
  c(167, "snow",  "snow crab",       "2024/25"),
  c(169, "nel",   "salmon",          "2025"),   # Alaska Peninsula sockeye
  c(169, "sap",   "salmon",          "2025"),
  c(171, "pws",   "salmon",          "2025"),   # set gillnet sockeye
  c(172, "yuk",   "salmon",          "2023"),   # Algaaciq Tribe request
  c(172, "yuk",   "salmon",          "2024"),
  c(172, "yuk",   "salmon",          "2025"),
  c(173, "sea",   "salmon",          "2025")),  # Metlakatla, Annette Island
  stringsAsFactors = FALSE)
names(ak.dis.pending) <- c("number", "area", "fishery", "season")
ak.dis.pending$status <- "pending"
ak.dis <- rbind(ak.dis, ak.dis.pending)
ak.dis$number <- as.integer(ak.dis$number)
# The mapped requests must be exactly NOAA's approved and pending ones
stopifnot(
  setequal(ak.dis$number[ak.dis$status == "declared"],
           ak.dis.noaa$number[ak.dis.noaa$determination == "Approved"]),
  setequal(ak.dis$number[ak.dis$status == "pending"],
           ak.dis.noaa$number[ak.dis.noaa$determination == "Pending"]))

# Area points (lon/lat, 0-360) and their labels. Labels are placed by hand
# in projected map units: the Bering Sea and Yukon-Kuskokwim areas in a
# right-aligned column over the western Bering Sea, the Alaska Peninsula and
# Gulf areas along the open North Pacific, Cook Inlet and Prince William
# Sound over interior Alaska. hjust 1 puts a label left of its position, 0
# right of it; the leader line ends at that position.
ak.dis.areas <- data.frame(rbind(
  c("kotz",  "Kotzebue Sound",                      197.4, 66.9, -0.040, -0.400, 0),
  c("nort",  "Norton Sound",                        196.5, 64.0, -0.205, -0.415, 1),
  c("yuk",   "Yukon River",                         197.0, 62.0, -0.205, -0.465, 1),
  c("kusk",  "Kuskokwim River",                     198.2, 60.8, -0.205, -0.500, 1),
  c("snow",  "Bering Sea",                          188.5, 59.0, -0.205, -0.535, 1),
  c("tann",  "Bering Sea",                          193.5, 56.8, -0.205, -0.578, 1),
  c("bbrkc", "Bristol Bay",                         197.5, 57.4, -0.205, -0.612, 1),
  c("nel",   "North Alaska Peninsula",              198.8, 56.0, -0.080, -0.650, 1),
  c("sap",   "South Alaska Peninsula",              199.6, 55.2, -0.080, -0.685, 1),
  c("chig",  "Chignik",                             201.6, 56.3, -0.030, -0.685, 0),
  c("kod",   "Kodiak",                              207.3, 57.5,  0.030, -0.655, 0),
  c("cod",   "Gulf of Alaska",                      211.0, 57.0,  0.165, -0.655, 0),
  c("lci",   "Lower Cook Inlet",                    208.3, 59.4,  0.075, -0.565, 0),
  c("uci",   "Upper Cook Inlet",                    208.7, 60.6,  0.000, -0.440, 0),
  c("pws",   "Prince William Sound / Copper River", 214.0, 60.5,  0.090, -0.475, 0),
  c("yak",   "Yakutat",                             220.3, 59.5,  0.190, -0.510, 0),
  c("sea",   "Southeast Alaska",                    225.5, 57.0,  0.160, -0.598, 0)),
  stringsAsFactors = FALSE)
names(ak.dis.areas) <- c("area", "place", "lon", "lat", "lx", "ly", "hjust")
ak.dis.areas[3:7] <- lapply(ak.dis.areas[3:7], as.numeric)

# Seasons as items with runs of consecutive years as ranges: 2020, 2021,
# 2022 -> 2020-22; 2021/22, 2022/23 -> 2021/22-2022/23
season_items <- function(s) {
  if (!length(s)) return(character(0))
  s <- sort(unique(s)); split.yr <- grepl("/", s[1])
  y <- as.integer(substr(s, 1, 4))
  runs <- split(seq_along(y), cumsum(c(1, diff(y) != 1)))
  unname(vapply(runs, function(i) {
    if (length(i) == 1) s[i]
    else if (split.yr) paste0(s[i[1]], "\u2013", s[i[length(i)]])
    else sprintf("%d\u2013%02d", y[i[1]], y[i[length(i)]] %% 100)
  }, ""))
}
# Declared seasons, then pending ones starred (in every area the pending
# seasons postdate the declared ones). wrap_line() breaks a label line
# longer than `width` characters after a comma.
season_list <- function(s, pending) {
  d <- unique(s[!pending]); p <- setdiff(unique(s[pending]), d)
  c(season_items(d), if (length(p)) paste0(season_items(p), "*"))
}
wrap_line <- function(x, width = 34) {
  if (nchar(x) <= width) return(x)
  cut <- gregexpr(", ", x)[[1]]
  cut <- cut[cut > 0 & cut + 1 <= width]
  if (!length(cut)) return(x)
  k <- max(cut)
  paste0(substr(x, 1, k), "\n", substr(x, k + 2, nchar(x)))
}

# One label per area: place name, then a line per fishery. An area's salmon
# disasters share one line, naming the species only when they all agree.
ak.dis.labs <- ak.dis %>%
  mutate(group = ifelse(grepl("salmon", fishery), "salmon", fishery)) %>%
  group_by(area, group) %>%
  summarise(fishery = if (n_distinct(fishery) == 1) fishery[1] else "salmon",
            seasons = paste(season_list(season, status == "pending"),
                            collapse = ", "), .groups = "drop") %>%
  mutate(line = vapply(paste(fishery, seasons), wrap_line, "")) %>%
  group_by(area) %>%
  summarise(detail = paste(line, collapse = "\n"), .groups = "drop") %>%
  inner_join(ak.dis.areas, by = "area")
ak.dis.labs[c("x", "y")] <- ps_project(ak.dis.labs$lon, ak.dis.labs$lat)
ak.dis.labs$seg.x <- ak.dis.labs$lx + ifelse(ak.dis.labs$hjust == 1, 0.004, -0.004)

# A rectangle of the polar stereographic projection, so the map fills the
# panel; land, graticule and study polygons are clipped by the panel
ak.xlim  <- c(-0.36, 0.30); ak.ylim <- c(-0.705, -0.385)
ak.coast <- ps_coast(mapWorld.clean, c(150, 250), c(45, 75))
ak.grat  <- bind_rows(
  bind_rows(lapply(seq(50, 70, 5), function(P)
    cbind(ps_path(c(150, 250), c(P, P), n = 400), grp = paste0("p", P)))),
  bind_rows(lapply(seq(160, 240, 10), function(L)
    cbind(ps_path(c(L, L), c(45, 75)), grp = paste0("m", L)))))

p.talk.ak.dis <- ggplot() +
  geom_polygon(data = lme.study, aes(x = x, y = y, group = area, fill = area),
               alpha = 0.3) +
  geom_polygon(data = ak.coast, aes(x = x, y = y, group = group),
               fill = "grey88", colour = "grey55", linewidth = 0.3) +
  geom_path(data = ak.grat, aes(x = x, y = y, group = grp),
            colour = "grey60", linewidth = 0.3, linetype = "dotted") +
  geom_path(data = lme.study, aes(x = x, y = y, group = area, colour = area),
            linewidth = 0.8) +
  geom_segment(data = ak.dis.labs, aes(x = x, y = y, xend = seg.x, yend = ly),
               colour = "grey30", linewidth = 0.45) +
  geom_point(data = ak.dis.labs, aes(x = x, y = y), shape = 21, size = 3,
             fill = "grey10", colour = "white", stroke = 0.8) +
  # Haloed, as the study-polygon outlines pass under some labels
  geom_text_halo(transmute(ak.dis.labs, x = lx, y = ly, lab = place, hjust),
                 size = 14 / .pt, r = 0.0012, n = 16, vjust = -0.15) +
  geom_text_halo(transmute(ak.dis.labs, x = lx, y = ly, lab = detail, hjust),
                 size = 13 / .pt, r = 0.0012, n = 16, fontface = "plain",
                 vjust = 1.15, lineheight = 0.95) +
  annotate("text", x = ak.xlim[1] + 0.006, y = ak.ylim[1] + 0.006,
           label = "* request pending", hjust = 0, vjust = 0, size = 13 / .pt) +
  scale_fill_manual(values = area.colors, guide = "none") +
  scale_colour_manual(values = area.colors, guide = "none") +
  coord_equal(xlim = ak.xlim, ylim = ak.ylim, expand = FALSE) +
  theme_grl_map() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        plot.margin  = margin(2, 2, 2, 2, "mm"))

talk_save(p.talk.ak.dis, "AK_fishery_disasters_map",
          width = talk_map_width(p.talk.ak.dis))
