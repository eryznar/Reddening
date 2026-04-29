# PURPOSE: Cellwise regression of SST AR(1) on AL SLP SD for the
# CESM2 FCM and MDM ensembles, mirroring the ERA5 analysis in
# analysis.R Section 13. Produces a 3-panel map: ERA5 | FCM | MDM.
#
# Per member workflow (matching observations Section 13):
#   1. SLP file -> AL-box (45-55N, 192.5-207.5E) area-weighted mean
#      monthly series -> NDJFM annual means (year = Jan year, full
#      5-month coverage only) -> 15-yr right-aligned rolling SD ->
#      z-scored AL.sd.
#   2. SST file (region lon 160-250E, lat 20-66N) -> per-cell monthly
#      climatology over the full record -> monthly anomalies -> per-
#      cell linear detrend -> NDJFM annual means (n == 5) -> 15-yr
#      rolling AR(1) per cell.
#   3. Per cell: GLS(y = SST AR(1), x = AL.sd, corAR1) -> beta.
#
# Per-member beta grids are cached in Output/CESM_member_betas/.
# Per-cell ensemble mean and SD across members are written to
# Output/cesm_{fcm,mdm}_sst_ar1_al_slp_sd_regression.csv.
#
# Time domain: full available range per ensemble (1850-2014).
# Significance overlay: deferred (this script only plots beta.mean).

source("./Scripts/load.libs.functions.R")
library(ncdf4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
library(nlme)
library(patchwork)
library(maps)

# ============================================================
# CONFIG
# ============================================================

fcm.sst.dir <- "./Data/CESM2 ensemble/SST/FCM"
fcm.slp.dir <- "./Data/CESM2 ensemble/SLP/FCM"
mdm.sst.dir <- "./Data/CESM2 ensemble/SST/MDM"
mdm.slp.dir <- "./Data/CESM2 ensemble/SLP/MDM"

mem.beta.dir <- "./Output/CESM_member_betas"
fcm.agg.cache <- "./Output/cesm_fcm_sst_ar1_al_slp_sd_regression.csv"
mdm.agg.cache <- "./Output/cesm_mdm_sst_ar1_al_slp_sd_regression.csv"
obs.cache     <- "./Output/sst_ar1_al_slp_sd_regression.csv"

fig.out <- "./Figures/SST_AR1_AL_SLP_SD_regression_obs_FCM_MDM.png"

dir.create(mem.beta.dir, showWarnings = FALSE, recursive = TRUE)

# Spatial domain (matches analysis.R Section 13)
lon.range <- c(160, 250)
lat.range <- c(20,  66)

# Aleutian Low box (Litzow 2020 PNAS)
al.lon <- c(192.5, 207.5)
al.lat <- c(45,    55)

# Rolling window
win.width <- 15

# CESM2-LE2 fixed time domain (T2 verified all files share length/calendar)
cesm.year.start <- 1850
cesm.year.end   <- 2014

# ============================================================
# HELPERS (member ID parsing, donor-coord fallback, etc.)
# ============================================================

fcm.id.re <- "LE2-[0-9]{4}\\.[0-9]{3}"

parse_id <- function(basename, ensemble) {
  if (ensemble == "FCM") {
    m <- regmatches(basename, regexpr(fcm.id.re, basename))
  } else {
    toks <- strsplit(basename, "\\.")[[1]]
    toks <- toks[grepl("^[0-9]{2}$", toks)]
    m <- if (length(toks)) toks[1] else character(0)
  }
  if (length(m) == 1 && nzchar(m)) m else NA_character_
}

list_members <- function(sst.dir, slp.dir, ensemble) {
  sst.files <- list.files(sst.dir, "\\.nc$", full.names = TRUE)
  slp.files <- list.files(slp.dir, "\\.nc$", full.names = TRUE)
  sst.ids <- vapply(basename(sst.files), parse_id, character(1), ensemble = ensemble)
  slp.ids <- vapply(basename(slp.files), parse_id, character(1), ensemble = ensemble)
  ids <- sort(intersect(sst.ids, slp.ids))
  data.frame(
    ensemble = ensemble,
    member   = ids,
    sst.file = sst.files[match(ids, sst.ids)],
    slp.file = slp.files[match(ids, slp.ids)],
    stringsAsFactors = FALSE
  )
}

# Donor-axis fallback: 2 FCM SST files (LE2-1301.019 / LE2-1301.020) ship
# with all-zero lon/lat coord arrays. Detect bad axes and substitute from
# a known-good sibling.
coord_looks_bad <- function(v, name) {
  if (is.null(v) || length(v) < 10) return(TRUE)
  if (all(v == 0, na.rm = TRUE)) return(TRUE)
  if (name == "lon" && max(v, na.rm = TRUE) < 180) return(TRUE)
  if (name == "lat" && max(v, na.rm = TRUE) < 60)  return(TRUE)
  FALSE
}

donor.coords <- NULL

read_coord <- function(nc, name) {
  v <- nc$dim[[name]]$vals
  if (coord_looks_bad(v, name)) {
    v <- tryCatch(as.numeric(ncvar_get(nc, name)), error = function(e) NULL)
  }
  if (coord_looks_bad(v, name)) {
    if (is.null(donor.coords))
      stop("read_coord: no donor axis available yet for '", name, "'")
    return(donor.coords[[name]])
  }
  as.numeric(v)
}

seed_donor_coords <- function(dirs) {
  for (d in dirs) {
    fs <- list.files(d, "\\.nc$", full.names = TRUE)
    for (f in fs) {
      nc <- tryCatch(nc_open(f), error = function(e) NULL); if (is.null(nc)) next
      lon <- nc$dim$lon$vals; lat <- nc$dim$lat$vals
      nc_close(nc)
      if (!coord_looks_bad(lon, "lon") && !coord_looks_bad(lat, "lat")) {
        return(list(lon = as.numeric(lon), lat = as.numeric(lat)))
      }
    }
  }
  stop("No netCDF file with valid lon/lat coordinates found.")
}

build_start_count <- function(nc, varname, sel) {
  dims  <- sapply(nc$var[[varname]]$dim, function(d) d$name)
  start <- integer(length(dims)); count <- integer(length(dims))
  for (k in seq_along(dims)) {
    nm <- dims[k]
    if (!is.null(sel[[nm]])) {
      start[k] <- sel[[nm]]$start
      count[k] <- sel[[nm]]$count
    } else {
      start[k] <- 1
      count[k] <- nc$dim[[nm]]$len
    }
  }
  list(start = start, count = count, dims = dims)
}

safe_ar1 <- function(v) {
  if (any(is.na(v))) return(NA_real_)
  suppressWarnings(acf(v, lag.max = 1, plot = FALSE)$acf[2])
}
roll_ar1_cell <- function(x, width = win.width) {
  rollapply(x, width = width, fill = NA, align = "right", FUN = safe_ar1)
}
roll_sd_safe <- function(x, width = win.width) {
  rollapply(x, width = width, fill = NA, align = "right",
            FUN = function(v) if (any(is.na(v))) NA_real_ else sd(v))
}

# ============================================================
# Per-member SLP -> AL.sd (15-yr rolling SD, z-scored)
# ============================================================

extract_al_slp_sd <- function(slp.file) {
  nc  <- nc_open(slp.file); on.exit(nc_close(nc), add = TRUE)
  lon <- read_coord(nc, "lon"); lat <- read_coord(nc, "lat")
  i.lon <- which(lon >= al.lon[1] & lon <= al.lon[2])
  i.lat <- which(lat >= al.lat[1] & lat <= al.lat[2])
  if (length(i.lon) == 0 || length(i.lat) == 0)
    stop("AL box selection empty for ", basename(slp.file))
  sel <- list(
    lon = list(start = min(i.lon), count = length(i.lon)),
    lat = list(start = min(i.lat), count = length(i.lat))
  )
  sc  <- build_start_count(nc, "PSL", sel)
  psl <- ncvar_get(nc, "PSL", start = sc$start, count = sc$count,
                   collapse_degen = FALSE)
  psl <- drop(psl)
  stopifnot(length(dim(psl)) == 3)
  w <- cos(lat[i.lat] * pi / 180)
  W <- matrix(rep(w, each = length(i.lon)), length(i.lon), length(i.lat))
  al.month <- apply(psl, 3, function(sl) {
    m <- !is.na(sl); if (!any(m)) NA_real_ else sum(sl[m] * W[m]) / sum(W[m])
  })
  n.t <- length(al.month)
  yrs <- rep(cesm.year.start:cesm.year.end, each = 12)[seq_len(n.t)]
  mos <- rep(1:12, times = (cesm.year.end - cesm.year.start + 1))[seq_len(n.t)]

  winter <- data.frame(year = yrs, month = mos, slp = al.month) %>%
    filter(month %in% c(11, 12, 1, 2, 3)) %>%
    mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
    group_by(win.year) %>%
    summarise(AL = mean(slp, na.rm = TRUE), n = n(), .groups = "drop") %>%
    filter(n == 5) %>% rename(year = win.year) %>% arrange(year)

  winter %>%
    mutate(AL.sd = roll_sd_safe(AL),
           AL.sd = as.numeric(scale(AL.sd))) %>%
    select(year, AL.sd) %>%
    filter(!is.na(AL.sd))
}

# ============================================================
# Per-member SST -> per-cell beta map
# ============================================================

extract_sst_grid <- function(sst.file) {
  nc  <- nc_open(sst.file); on.exit(nc_close(nc), add = TRUE)
  lon <- read_coord(nc, "lon"); lat <- read_coord(nc, "lat")

  i.lon <- which(lon >= lon.range[1] & lon <= lon.range[2])
  i.lat <- which(lat >= lat.range[1] & lat <= lat.range[2])
  if (length(i.lon) == 0 || length(i.lat) == 0)
    stop("SST domain selection empty for ", basename(sst.file))

  sel <- list(
    lon  = list(start = min(i.lon), count = length(i.lon)),
    lat  = list(start = min(i.lat), count = length(i.lat)),
    z_t  = list(start = 1,          count = 1)   # surface only if present
  )
  sc  <- build_start_count(nc, "SST", sel)
  sst <- ncvar_get(nc, "SST", start = sc$start, count = sc$count,
                   collapse_degen = FALSE)
  sst[sst > 1e20] <- NA_real_
  sst <- drop(sst)   # (lon, lat, time)
  stopifnot(length(dim(sst)) == 3)

  list(lon = lon[i.lon], lat = lat[i.lat], sst = sst)
}

process_member <- function(sst.file, slp.file, ensemble, member) {
  cache.file <- file.path(mem.beta.dir,
                          sprintf("%s_%s.csv", ensemble, member))
  if (file.exists(cache.file)) return(read.csv(cache.file))

  cat(sprintf("  [%s/%s] AL.sd ... ", ensemble, member))
  al.df <- extract_al_slp_sd(slp.file)
  cat(sprintf("years=%d-%d (%d) | SST grid ... ",
              min(al.df$year), max(al.df$year), nrow(al.df)))

  sst.g <- extract_sst_grid(sst.file)
  lon.vec <- sst.g$lon; lat.vec <- sst.g$lat
  nx <- length(lon.vec); ny <- length(lat.vec); nt <- dim(sst.g$sst)[3]
  yrs <- rep(cesm.year.start:cesm.year.end, each = 12)[seq_len(nt)]
  mos <- rep(1:12, times = (cesm.year.end - cesm.year.start + 1))[seq_len(nt)]

  sst.mat <- matrix(sst.g$sst, nrow = nx * ny, ncol = nt)
  rm(sst.g); gc()

  # Monthly climatology per cell over full record
  clim <- matrix(NA_real_, nrow = nx * ny, ncol = 12)
  for (m in 1:12) {
    cols <- which(mos == m)
    clim[, m] <- rowMeans(sst.mat[, cols, drop = FALSE], na.rm = TRUE)
  }
  anom <- sst.mat - clim[, mos]
  rm(sst.mat); gc()

  # NDJFM annual means with full 5-month coverage
  is.win  <- mos %in% c(11, 12, 1, 2, 3)
  win.yr  <- ifelse(mos %in% c(11, 12), yrs + 1L, yrs)
  uy      <- sort(unique(win.yr[is.win]))
  ann     <- matrix(NA_real_, nrow = nx * ny, ncol = length(uy))
  n.mo    <- integer(length(uy))
  for (k in seq_along(uy)) {
    cols <- which(is.win & win.yr == uy[k])
    n.mo[k] <- length(cols)
    if (length(cols) == 5)
      ann[, k] <- rowMeans(anom[, cols, drop = FALSE], na.rm = TRUE)
  }
  uy.full <- uy[n.mo == 5]
  ann     <- ann[, n.mo == 5, drop = FALSE]
  rm(anom); gc()

  # Per-cell linear detrend of the annual NDJFM series
  detr <- ann
  for (k in seq_len(nrow(detr))) {
    y <- detr[k, ]
    if (sum(!is.na(y)) < 10) { detr[k, ] <- NA_real_; next }
    detr[k, ] <- residuals(lm(y ~ seq_along(y), na.action = na.exclude))
  }
  rm(ann); gc()

  # Per-cell 15-yr rolling AR(1)
  ar1.mat <- matrix(NA_real_, nrow = nrow(detr), ncol = ncol(detr))
  for (k in seq_len(nrow(detr))) {
    if (any(is.na(detr[k, ]))) next
    ar1.mat[k, ] <- roll_ar1_cell(detr[k, ])
  }
  rm(detr); gc()

  # Align to AL.sd years
  yr.use <- intersect(uy.full, al.df$year)
  if (length(yr.use) < win.width) {
    cat("INSUFFICIENT YEARS\n")
    return(NULL)
  }
  x.vec  <- al.df$AL.sd[match(yr.use, al.df$year)]
  ar1.al <- ar1.mat[, match(yr.use, uy.full), drop = FALSE]
  rm(ar1.mat); gc()

  good.cells <- which(rowSums(!is.na(ar1.al)) >= 10)
  cat(sprintf("GLS @ %d cells ... ", length(good.cells)))

  beta.vec <- rep(NA_real_, nx * ny)
  df.reg   <- data.frame(y = NA_real_, x = x.vec)
  for (k in good.cells) {
    y <- ar1.al[k, ]
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
  }

  lonlat <- expand.grid(lon = lon.vec, lat = lat.vec,
                        KEEP.OUT.ATTRS = FALSE)
  out <- data.frame(lon = lonlat$lon, lat = lonlat$lat, beta = beta.vec) %>%
    filter(!is.na(beta))

  write.csv(out, cache.file, row.names = FALSE)
  cat("done.\n")
  out
}

# ============================================================
# RUN
# ============================================================

cat("Seeding coordinate-axis donor...\n")
donor.coords <- seed_donor_coords(c(fcm.sst.dir, mdm.sst.dir,
                                    fcm.slp.dir, mdm.slp.dir))
cat("Donor lon len:", length(donor.coords$lon),
    " lat len:", length(donor.coords$lat), "\n")

cat.df <- bind_rows(
  list_members(fcm.sst.dir, fcm.slp.dir, "FCM"),
  list_members(mdm.sst.dir, mdm.slp.dir, "MDM")
)
cat("FCM members:", sum(cat.df$ensemble == "FCM"),
    " | MDM members:", sum(cat.df$ensemble == "MDM"), "\n")

aggregate_ensemble <- function(ensemble, agg.cache) {
  if (file.exists(agg.cache)) {
    cat(sprintf("Loading cached %s aggregate: %s\n", ensemble, agg.cache))
    return(read.csv(agg.cache))
  }
  rows <- cat.df %>% filter(ensemble == !!ensemble)
  cat(sprintf("\n=== Processing %s (%d members) ===\n", ensemble, nrow(rows)))
  member.list <- list()
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    res <- tryCatch(process_member(r$sst.file, r$slp.file,
                                   ensemble, r$member),
                    error = function(e) {
                      cat("ERROR:", conditionMessage(e), "\n"); NULL
                    })
    if (!is.null(res)) {
      res$member <- r$member
      member.list[[length(member.list) + 1]] <- res
    }
  }
  all.df <- bind_rows(member.list)
  agg <- all.df %>%
    group_by(lon, lat) %>%
    summarise(beta.mean  = mean(beta, na.rm = TRUE),
              beta.sd    = sd(beta,   na.rm = TRUE),
              n.members  = sum(!is.na(beta)),
              .groups    = "drop")
  write.csv(agg, agg.cache, row.names = FALSE)
  cat(sprintf("Saved: %s (%d cells)\n", agg.cache, nrow(agg)))
  agg
}

fcm.agg <- aggregate_ensemble("FCM", fcm.agg.cache)
mdm.agg <- aggregate_ensemble("MDM", mdm.agg.cache)

# ============================================================
# PLOT
# ============================================================

cat("\n=== Building 3-panel figure ===\n")
if (!file.exists(obs.cache))
  stop("Observation regression cache not found at ", obs.cache,
       ". Run analysis.R Section 13 first.")
obs.df <- read.csv(obs.cache) %>% filter(lon >= 160)

# Two color scales: observations on their own (larger magnitude) and
# FCM+MDM sharing a common (smaller magnitude) CESM scale.
col.lim.obs  <- max(abs(obs.df$beta), na.rm = TRUE)
col.lim.cesm <- max(c(abs(fcm.agg$beta.mean),
                      abs(mdm.agg$beta.mean)), na.rm = TRUE)

mapWorld.clean <- map_data("world2")
al.box.df <- data.frame(
  x = c(192.5, 207.5, 207.5, 192.5, 192.5),
  y = c(45,    45,    55,    55,    45)
)

lon_label <- function(x) {
  ifelse(x <= 180, paste0(x, "\u00b0E"),
         paste0(360 - x, "\u00b0W"))
}

make_panel <- function(df, beta.col, title.text, col.lim,
                       legend.name = "\u03b2") {
  ggplot() +
    geom_tile(data = df, aes(x = lon, y = lat,
                              fill = .data[[beta.col]])) +
    geom_polygon(data = mapWorld.clean,
                 aes(x = long, y = lat, group = group),
                 fill = "gray85", color = "gray30", linewidth = 0.25) +
    geom_path(data = al.box.df, aes(x = x, y = y),
              color = "black", linewidth = 0.7, linetype = "dashed",
              inherit.aes = FALSE) +
    scale_fill_gradient2(low = "steelblue4", mid = "white",
                         high = "darkorange4", midpoint = 0,
                         limits = c(-col.lim, col.lim),
                         name = legend.name) +
    coord_cartesian(xlim = c(160, 250), ylim = c(20, 66)) +
    scale_x_continuous(breaks = seq(160, 250, 20), labels = lon_label) +
    scale_y_continuous(breaks = seq(20, 60, 10),
                       labels = function(y) paste0(y, "\u00b0N")) +
    labs(title = title.text, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
          legend.position  = "right")
}

# Observations: own scale (legend name distinguishes it from CESM panels
# so patchwork keeps both legends when guides are collected).
p.obs <- make_panel(obs.df, "beta", "ERA5 (obs)",
                    col.lim = col.lim.obs,
                    legend.name = "\u03b2 (obs)")

# FCM + MDM: shared CESM scale (identical legend name -> patchwork
# collapses to a single shared legend).
p.fcm <- make_panel(fcm.agg, "beta.mean",
                    sprintf("FCM ensemble mean (n=%d)",
                            max(fcm.agg$n.members, na.rm = TRUE)),
                    col.lim = col.lim.cesm,
                    legend.name = "\u03b2 (CESM2)")
p.mdm <- make_panel(mdm.agg, "beta.mean",
                    sprintf("MDM ensemble mean (n=%d)",
                            max(mdm.agg$n.members, na.rm = TRUE)),
                    col.lim = col.lim.cesm,
                    legend.name = "\u03b2 (CESM2)")

p.all <- (p.obs | p.fcm | p.mdm) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "SST AR(1) regressed on AL SLP SD (15-yr NDJFM rolling)",
    subtitle = "Per-member GLS AR(1); ensemble mean \u03b2; obs and CESM2 on separate scales",
    theme    = theme(plot.title    = element_text(hjust = 0.5, size = 12),
                     plot.subtitle = element_text(hjust = 0.5, size = 9))
  )

ggsave(fig.out, plot = p.all, width = 14, height = 5, dpi = 300)
cat("Saved:", fig.out, "\n")

# ============================================================
# FCM vs MDM significance test
# ------------------------------------------------------------
# Three-panel plot:
#   1) FCM ensemble mean beta
#   2) MDM ensemble mean beta
#   3) FCM - MDM difference, with cells where p <= 0.05 (after
#      Benjamini-Hochberg FDR) outlined.
#
# Test: cellwise two-sample Welch's t-test using the cached
# ensemble mean, SD, and n at each cell. No autocorrelation
# correction (treats members as the units of replication).
# ============================================================

cat("\n=== FCM vs MDM significance test ===\n")

fig.diff <- "./Figures/SST_AR1_AL_SLP_SD_FCM_MDM_difference.png"

cmp <- inner_join(
  fcm.agg %>% select(lon, lat,
                     fcm.mean = beta.mean,
                     fcm.sd   = beta.sd,
                     fcm.n    = n.members),
  mdm.agg %>% select(lon, lat,
                     mdm.mean = beta.mean,
                     mdm.sd   = beta.sd,
                     mdm.n    = n.members),
  by = c("lon", "lat")
) %>%
  mutate(diff   = fcm.mean - mdm.mean,
         se     = sqrt(fcm.sd^2 / fcm.n + mdm.sd^2 / mdm.n),
         t.stat = diff / se,
         df     = (fcm.sd^2 / fcm.n + mdm.sd^2 / mdm.n)^2 /
                  ((fcm.sd^2 / fcm.n)^2 / (fcm.n - 1) +
                   (mdm.sd^2 / mdm.n)^2 / (mdm.n - 1)),
         p.raw  = 2 * pt(-abs(t.stat), df = df))

cmp$p.fdr <- p.adjust(cmp$p.raw, method = "BH")
cmp$sig   <- !is.na(cmp$p.fdr) & cmp$p.fdr <= 0.05

cat(sprintf("  cells compared: %d | sig (FDR q<=0.05): %d (%.1f%%)\n",
            nrow(cmp), sum(cmp$sig, na.rm = TRUE),
            100 * mean(cmp$sig, na.rm = TRUE)))

col.lim.diff <- max(abs(cmp$diff), na.rm = TRUE)

p.fcm2 <- make_panel(fcm.agg, "beta.mean",
                     sprintf("FCM ensemble mean (n=%d)",
                             max(fcm.agg$n.members, na.rm = TRUE)),
                     col.lim = col.lim.cesm,
                     legend.name = "\u03b2 (CESM2)")
p.mdm2 <- make_panel(mdm.agg, "beta.mean",
                     sprintf("MDM ensemble mean (n=%d)",
                             max(mdm.agg$n.members, na.rm = TRUE)),
                     col.lim = col.lim.cesm,
                     legend.name = "\u03b2 (CESM2)")

p.diff <- make_panel(cmp, "diff", "FCM \u2212 MDM",
                     col.lim = col.lim.diff,
                     legend.name = "\u0394\u03b2") +
  geom_point(data = filter(cmp, sig),
             aes(x = lon, y = lat),
             shape = 16, size = 0.25, color = "black",
             inherit.aes = FALSE)

p.cmp <- (p.fcm2 | p.mdm2 | p.diff) +
  plot_annotation(
    title    = "FCM vs MDM: SST AR(1) ~ AL SLP SD ensemble means",
    subtitle = "Difference panel: stippled cells significant at FDR q \u2264 0.05 (Welch t-test, no AR correction)",
    theme    = theme(plot.title    = element_text(hjust = 0.5, size = 12),
                     plot.subtitle = element_text(hjust = 0.5, size = 9))
  )

ggsave(fig.diff, plot = p.cmp, width = 14, height = 5, dpi = 300)
cat("Saved:", fig.diff, "\n")

# ============================================================
# System comparison: per-member area-weighted mean beta in the
# GOA and EBS polygons; Welch's t-test FCM vs MDM for each.
# ============================================================

cat("\n=== System comparison (GOA / EBS area means) ===\n")

# Polygons (lon in 0-360E to match CESM grid)
goa.x.poly <- c(201, 201, 205, 208, 225, 231, 201)
goa.y.poly <- c(55,  56.5, 59,  61,  61,  55,  55)
ebs.x.poly <- c(183, 183, 203, 203, 191)
ebs.y.poly <- c(53,  65,  65,  57.5, 53)

area_weighted_mean <- function(df, poly.x, poly.y) {
  inside <- sp::point.in.polygon(df$lon, df$lat, poly.x, poly.y) > 0
  d <- df[inside, ]
  good <- !is.na(d$beta) & is.finite(d$beta)
  if (!any(good)) return(NA_real_)
  w <- cos(d$lat[good] * pi / 180)
  sum(d$beta[good] * w) / sum(w)
}

mem.files <- list.files(mem.beta.dir, "\\.csv$", full.names = TRUE)
sys.df <- lapply(mem.files, function(f) {
  bn    <- tools::file_path_sans_ext(basename(f))
  parts <- strsplit(bn, "_", fixed = TRUE)[[1]]
  ens   <- parts[1]
  mem   <- paste(parts[-1], collapse = "_")
  d     <- read.csv(f)
  data.frame(ensemble = ens, member = mem,
             GOA = area_weighted_mean(d, goa.x.poly, goa.y.poly),
             EBS = area_weighted_mean(d, ebs.x.poly, ebs.y.poly))
}) %>% bind_rows()

cat(sprintf("  members: FCM=%d, MDM=%d\n",
            sum(sys.df$ensemble == "FCM"),
            sum(sys.df$ensemble == "MDM")))

t.goa <- t.test(GOA ~ ensemble, data = sys.df)
t.ebs <- t.test(EBS ~ ensemble, data = sys.df)
cat("\nGOA Welch t-test (FCM vs MDM):\n"); print(t.goa)
cat("\nEBS Welch t-test (FCM vs MDM):\n"); print(t.ebs)

sys.long <- sys.df %>%
  pivot_longer(c(GOA, EBS), names_to = "region", values_to = "beta") %>%
  mutate(region = factor(region, levels = c("EBS", "GOA")))

t.lab <- data.frame(
  region = factor(c("EBS", "GOA"), levels = c("EBS", "GOA")),
  label  = c(sprintf("t = %.2f, df = %.1f, p = %.3f",
                     t.ebs$statistic, t.ebs$parameter, t.ebs$p.value),
             sprintf("t = %.2f, df = %.1f, p = %.3f",
                     t.goa$statistic, t.goa$parameter, t.goa$p.value))
)

p.sys <- ggplot(sys.long, aes(ensemble, beta)) +
  geom_boxplot(outlier.shape = NA, fill = "gray90", width = 0.5) +
  geom_jitter(width = 0.12, alpha = 0.55, size = 1.1) +
  facet_wrap(~ region, scales = "free_y") +
  geom_text(data = t.lab,
            aes(x = 1.5, y = Inf, label = label),
            vjust = 1.6, size = 3.4, inherit.aes = FALSE) +
  labs(x = NULL, y = "Area-weighted mean \u03b2",
       title    = "Per-member SST AR(1) ~ AL SLP SD: GOA & EBS area means",
       subtitle = "Cos(lat)-weighted mean of cellwise \u03b2 inside each polygon") +
  theme_bw(base_size = 11) +
  theme(plot.title    = element_text(hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

fig.sys <- "./Figures/SST_AR1_AL_SLP_SD_FCM_MDM_system_comparison.png"
ggsave(fig.sys, plot = p.sys, width = 8, height = 4.5, dpi = 300)
write.csv(sys.df,
          "./Output/cesm_fcm_mdm_system_area_means.csv",
          row.names = FALSE)
cat("Saved:", fig.sys, "\n")
