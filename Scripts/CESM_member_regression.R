# PURPOSE: CESM2 per-member SST-SLP linkage verification (front-end).
#
#   Before doing any per-member SLP -> SST regression with these files, we
#   must be confident that the SST and SLP netCDF files sharing a member ID
#   actually come from the same integration. Filenames encode IDs
#     FCM: LE2-####.###   (macro ocean-IC year . micro atmos-seed)
#     MDM: ##             (01-20)
#   but filename-based pairing could in principle be wrong.
#
#   This script runs three verification tests. The member-level regression
#   machinery should only be added *after* all three pass.
#
#   T1 — Metadata cross-check: parse internal netCDF global attributes
#        (history, case) and confirm the embedded ID matches the filename ID.
#   T2 — Time-axis alignment: confirm paired SST/SLP share identical time
#        length and calendar.
#   T3 — ENSO fingerprint permutation test: within each ensemble, the real
#        SST-SLP pairing should give a strongly negative Nino3.4 vs SOI
#        correlation (FCM, where wind-ocean coupling is active) or no
#        systematic signal (MDM, where the ocean sees climatological wind
#        stress). Random shuffled pairings should decorrelate in FCM.
#
# Outputs:
#   Output/CESM_member_T1_metadata_check.csv
#   Output/CESM_member_T2_time_axis_check.csv
#   Output/CESM_member_T3_tropical_indices.csv
#   Output/CESM_member_T3_ENSO_permutation.csv
#   Figures/CESM_member_T3_ENSO_permutation.png

source("./Scripts/load.libs.functions.R")
library(ncdf4)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# CONFIG
# ============================================================

fcm.sst.dir <- "./Data/CESM2 ensemble/SST/FCM"
fcm.slp.dir <- "./Data/CESM2 ensemble/SLP/FCM"
mdm.sst.dir <- "./Data/CESM2 ensemble/SST/MDM"
mdm.slp.dir <- "./Data/CESM2 ensemble/SLP/MDM"

t1.out <- "./Output/CESM_member_T1_metadata_check.csv"
t2.out <- "./Output/CESM_member_T2_time_axis_check.csv"
t3.cache <- "./Output/CESM_member_T3_tropical_indices.csv"
t3.out   <- "./Output/CESM_member_T3_ENSO_permutation.csv"
t3.fig   <- "./Figures/CESM_member_T3_ENSO_permutation.png"

n.perm <- 999
set.seed(42)

# Nino-3.4 box (5N-5S, 170W-120W = 190-240E)
nino.lon <- c(190, 240)
nino.lat <- c(-5,   5)

# SOI stations — nearest grid point to Tahiti and Darwin
tahiti.lonlat <- c(210.6, -17.65)
darwin.lonlat <- c(130.85, -12.46)

# ============================================================
# HELPERS
# ============================================================

fcm.id.re <- "LE2-[0-9]{4}\\.[0-9]{3}"
mdm.id.re <- "(?<![0-9])[0-9]{2}(?![0-9])"

parse_id <- function(basename, ensemble) {
  if (ensemble == "FCM") {
    m <- regmatches(basename, regexpr(fcm.id.re, basename))
  } else {
    # MDM files: take the first 2-digit group that is not inside a date or extension
    # e.g. cesm2.MD.postP.01.SST.185001-201412.nc  -> 01
    #      cesm2.MD.01.PSL.185001-201412.nc        -> 01
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

get_global_att <- function(nc, name) {
  a <- ncatt_get(nc, 0, name)
  if (isTRUE(a$hasatt)) a$value else NA_character_
}

# ============================================================
# Catalog
# ============================================================

cat("== Cataloguing files ==\n")
cat.df <- bind_rows(
  list_members(fcm.sst.dir, fcm.slp.dir, "FCM"),
  list_members(mdm.sst.dir, mdm.slp.dir, "MDM")
)
cat("FCM pairs:", sum(cat.df$ensemble == "FCM"), "\n")
cat("MDM pairs:", sum(cat.df$ensemble == "MDM"), "\n")

# ============================================================
# T1 — metadata cross-check
# ============================================================
# For each file, open, pull history + case attributes, and verify the
# filename member ID appears verbatim in the internal attribute string.
# FCM SST files have no global attributes (post-processed), so for those
# we record "no metadata; filename-only" — NOT a failure, but flagged.

cat("\n== T1: metadata cross-check ==\n")
check_metadata <- function(file, ensemble, mem.id) {
  nc <- nc_open(file)
  on.exit(nc_close(nc), add = TRUE)
  hist <- get_global_att(nc, "history")
  case <- get_global_att(nc, "case")
  attrs <- paste(na.omit(c(hist, case)), collapse = " | ")
  if (!nzchar(attrs)) {
    return(list(status = "no_metadata", attrs = ""))
  }
  pattern <- if (ensemble == "FCM") mem.id else paste0("[^0-9]", mem.id, "[^0-9]")
  ok <- grepl(pattern, attrs)
  list(status = if (ok) "match" else "MISMATCH", attrs = substr(attrs, 1, 400))
}

t1.rows <- list()
for (i in seq_len(nrow(cat.df))) {
  r <- cat.df[i, ]
  sst.chk <- check_metadata(r$sst.file, r$ensemble, r$member)
  slp.chk <- check_metadata(r$slp.file, r$ensemble, r$member)
  t1.rows[[i]] <- data.frame(
    ensemble = r$ensemble, member = r$member,
    sst.status = sst.chk$status, slp.status = slp.chk$status,
    sst.attrs = sst.chk$attrs, slp.attrs = slp.chk$attrs,
    stringsAsFactors = FALSE
  )
}
t1.df <- bind_rows(t1.rows)

t1.summary <- t1.df %>%
  group_by(ensemble, sst.status, slp.status) %>%
  summarise(n = n(), .groups = "drop")
cat("T1 summary:\n"); print(t1.summary)
if (any(t1.df$sst.status == "MISMATCH") || any(t1.df$slp.status == "MISMATCH")) {
  warning("T1 FAIL: at least one file has internal metadata that disagrees with its filename ID.")
} else {
  cat("T1 PASS: no filename-vs-metadata mismatches.\n")
}
write.csv(t1.df, t1.out, row.names = FALSE)

# ============================================================
# T2 — time-axis alignment
# ============================================================

cat("\n== T2: time-axis alignment ==\n")
get_time <- function(file) {
  nc <- nc_open(file); on.exit(nc_close(nc), add = TRUE)
  t.vals <- tryCatch(ncvar_get(nc, "time"), error = function(e) NA_real_)
  cal    <- tryCatch(ncatt_get(nc, "time", "calendar")$value, error = function(e) NA_character_)
  units  <- tryCatch(ncatt_get(nc, "time", "units")$value,    error = function(e) NA_character_)
  list(n = length(t.vals), calendar = cal, units = units,
       first = t.vals[1], last = t.vals[length(t.vals)])
}

t2.rows <- list()
for (i in seq_len(nrow(cat.df))) {
  r <- cat.df[i, ]
  s <- get_time(r$sst.file); p <- get_time(r$slp.file)
  ok.n   <- isTRUE(s$n == p$n)
  ok.cal <- isTRUE(s$calendar == p$calendar)
  t2.rows[[i]] <- data.frame(
    ensemble = r$ensemble, member = r$member,
    sst.n = s$n, slp.n = p$n, same.length = ok.n,
    sst.cal = s$calendar, slp.cal = p$calendar, same.cal = ok.cal,
    sst.units = s$units, slp.units = p$units,
    stringsAsFactors = FALSE
  )
}
t2.df <- bind_rows(t2.rows)
cat("T2 summary (length + calendar):\n")
print(t2.df %>% group_by(ensemble, same.length, same.cal) %>% summarise(n = n(), .groups = "drop"))
if (!all(t2.df$same.length) || !all(t2.df$same.cal)) {
  warning("T2 FAIL: some SST-SLP pairs disagree on time-axis length or calendar.")
} else {
  cat("T2 PASS: all pairs share length and calendar.\n")
}
write.csv(t2.df, t2.out, row.names = FALSE)

# ============================================================
# T3 — ENSO fingerprint permutation test
# ============================================================
# Build per-member annual Nino-3.4 (from SST) and SOI-proxy (Tahiti-Darwin
# standardised PSL anomalies, from SLP). ENSO guarantees that same-integration
# Nino-3.4 and SOI are strongly anti-correlated at annual resolution in FCM.

cat("\n== T3: tropical-index extraction + permutation test ==\n")

nearest_idx <- function(coord, target) which.min(abs(coord - target))

# Read a 1-D coordinate axis. A couple of CESM2-LE FCM SST files have their
# lon/lat arrays stored as all-zeros; detect that and fall back to a donor
# axis taken from a known-good sibling file (all FCM/MDM SST/SLP files share
# the fv0.9x1.25 grid).
coord_looks_bad <- function(v, name) {
  if (is.null(v) || length(v) < 10) return(TRUE)
  if (all(v == 0, na.rm = TRUE)) return(TRUE)
  if (name == "lon" && max(v, na.rm = TRUE) < 180) return(TRUE)
  if (name == "lat" && max(v, na.rm = TRUE) < 60)  return(TRUE)
  FALSE
}

donor.coords <- NULL   # populated lazily, first good file wins

read_coord <- function(nc, name) {
  v <- nc$dim[[name]]$vals
  if (coord_looks_bad(v, name)) {
    v <- tryCatch(as.numeric(ncvar_get(nc, name)), error = function(e) NULL)
  }
  if (coord_looks_bad(v, name)) {
    if (is.null(donor.coords)) stop("read_coord: no donor axis available yet for '", name, "'")
    return(donor.coords[[name]])
  }
  as.numeric(v)
}

# Walk the ensemble directories once at startup; cache the first file whose
# lon and lat both pass coord_looks_bad() as the donor.
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

cat("Seeding coordinate-axis donor...\n")
donor.coords <- seed_donor_coords(c(fcm.sst.dir, mdm.sst.dir, fcm.slp.dir, mdm.slp.dir))
cat("Donor lon len:", length(donor.coords$lon),
    " lat len:", length(donor.coords$lat), "\n")

# Build start/count vectors from the file's actual dim order, so the code is
# robust if a file stores dims in a different order than expected.
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

extract_nino34 <- function(sst.file) {
  nc <- nc_open(sst.file); on.exit(nc_close(nc), add = TRUE)
  lon <- as.numeric(read_coord(nc, "lon"))
  lat <- as.numeric(read_coord(nc, "lat"))
  i.lon <- which(lon >= nino.lon[1] & lon <= nino.lon[2])
  i.lat <- which(lat >= nino.lat[1] & lat <= nino.lat[2])
  if (length(i.lon) == 0 || length(i.lat) == 0) {
    stop("Nino-3.4 selection empty. lon range=[",
         round(min(lon),2), ",", round(max(lon),2), "], lat range=[",
         round(min(lat),2), ",", round(max(lat),2), "]")
  }
  sel <- list(
    lon  = list(start = min(i.lon), count = length(i.lon)),
    lat  = list(start = min(i.lat), count = length(i.lat)),
    z_t  = list(start = 1,          count = 1)
  )
  sc  <- build_start_count(nc, "SST", sel)
  sst <- ncvar_get(nc, "SST", start = sc$start, count = sc$count,
                   collapse_degen = FALSE)
  sst[sst > 1e20] <- NA_real_
  sst <- drop(sst)   # strip z_t singleton; leaves (lon, lat, time)
  stopifnot(length(dim(sst)) == 3)
  w <- cos(lat[i.lat] * pi / 180)
  W <- matrix(rep(w, each = length(i.lon)), length(i.lon), length(i.lat))
  apply(sst, 3, function(sl) {
    m <- !is.na(sl)
    if (!any(m)) NA_real_ else sum(sl[m] * W[m]) / sum(W[m])
  })
}

extract_soi <- function(slp.file) {
  nc <- nc_open(slp.file); on.exit(nc_close(nc), add = TRUE)
  lon <- as.numeric(read_coord(nc, "lon"))
  lat <- as.numeric(read_coord(nc, "lat"))
  it <- nearest_idx(lon, tahiti.lonlat[1]); jt <- nearest_idx(lat, tahiti.lonlat[2])
  id <- nearest_idx(lon, darwin.lonlat[1]); jd <- nearest_idx(lat, darwin.lonlat[2])
  sel.t <- list(lon = list(start = it, count = 1),
                lat = list(start = jt, count = 1))
  sel.d <- list(lon = list(start = id, count = 1),
                lat = list(start = jd, count = 1))
  sc.t <- build_start_count(nc, "PSL", sel.t)
  sc.d <- build_start_count(nc, "PSL", sel.d)
  psl.t <- ncvar_get(nc, "PSL", start = sc.t$start, count = sc.t$count,
                     collapse_degen = FALSE)
  psl.d <- ncvar_get(nc, "PSL", start = sc.d$start, count = sc.d$count,
                     collapse_degen = FALSE)
  list(tahiti = as.numeric(psl.t), darwin = as.numeric(psl.d))
}

monthly_to_annual <- function(x, years) {
  data.frame(year = years, x = x) %>%
    group_by(year) %>% summarise(x = mean(x, na.rm = TRUE), .groups = "drop")
}

standardise_anom <- function(x, months) {
  clim <- tapply(x, months, mean, na.rm = TRUE)
  sdev <- tapply(x, months, sd,   na.rm = TRUE)
  (x - clim[as.character(months)]) / sdev[as.character(months)]
}

if (file.exists(t3.cache)) {
  cat("Loading cached T3 tropical indices: ", t3.cache, "\n")
  t3.df <- read.csv(t3.cache)
} else {
  cat("Extracting Nino-3.4 and SOI-proxy per member (this takes a few minutes)...\n")
  n.t <- 1980
  years.m  <- rep(1850:2014, each = 12)[seq_len(n.t)]
  months.m <- rep(1:12, times = 165)[seq_len(n.t)]

  rows <- list()
  for (i in seq_len(nrow(cat.df))) {
    r <- cat.df[i, ]
    cat("  ", r$ensemble, r$member, "\n")
    nino <- extract_nino34(r$sst.file)
    soi  <- extract_soi(r$slp.file)
    if (length(nino) != n.t) next
    if (length(soi$tahiti) != n.t) next
    t.an <- standardise_anom(soi$tahiti, months.m)
    d.an <- standardise_anom(soi$darwin, months.m)
    soi.idx <- t.an - d.an
    nino.an <- monthly_to_annual(nino, years.m)
    soi.an  <- monthly_to_annual(soi.idx, years.m)
    rows[[i]] <- data.frame(
      ensemble = r$ensemble, member = r$member,
      year = nino.an$year, nino34 = nino.an$x, soi = soi.an$x
    )
  }
  t3.df <- bind_rows(rows)
  write.csv(t3.df, t3.cache, row.names = FALSE)
  cat("Saved: ", t3.cache, "\n")
}

# Correlate Nino3.4 with SOI, 1950-2014 (avoid spin-up decades)
t3.use <- t3.df %>% filter(year >= 1950, year <= 2014)

real.corr <- t3.use %>%
  group_by(ensemble, member) %>%
  summarise(r = cor(nino34, soi, use = "complete.obs"), .groups = "drop")

cat("\nReal-pairing Nino3.4 vs SOI correlations:\n")
print(real.corr %>% group_by(ensemble) %>%
        summarise(n        = n(),
                  n.na     = sum(is.na(r)),
                  mean.r   = mean(r, na.rm = TRUE),
                  median.r = median(r, na.rm = TRUE),
                  min.r    = min(r, na.rm = TRUE),
                  max.r    = max(r, na.rm = TRUE),
                  .groups  = "drop"))

# Permutation test: within each ensemble, shuffle the SOI-member labels
# (i.e., re-pair SST of member A with SOI of a random other member),
# recompute the mean r across all pairings. Repeat 999 times.
# Statistic: mean correlation across re-paired members.
cat("\nRunning permutation test (", n.perm, "shuffles per ensemble)...\n")

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}

perm_stat <- function(df) {
  # df has columns member, year, nino34, soi
  members <- unique(df$member)
  nino.wide <- df %>% select(member, year, nino34) %>%
    pivot_wider(names_from = member, values_from = nino34)
  soi.wide  <- df %>% select(member, year, soi) %>%
    pivot_wider(names_from = member, values_from = soi)
  real.r <- sapply(members, function(m) safe_cor(nino.wide[[m]], soi.wide[[m]]))
  real.mean <- mean(real.r, na.rm = TRUE)
  perm.means <- replicate(n.perm, {
    p <- sample(members)
    r <- sapply(seq_along(members), function(k)
      safe_cor(nino.wide[[members[k]]], soi.wide[[p[k]]]))
    mean(r, na.rm = TRUE)
  })
  list(real = real.mean, perm = perm.means, real.r = real.r)
}

perm.out <- lapply(split(t3.use, t3.use$ensemble), perm_stat)

t3.pval <- data.frame(
  ensemble   = names(perm.out),
  real.mean  = sapply(perm.out, function(x) x$real),
  perm.mean  = sapply(perm.out, function(x) mean(x$perm)),
  perm.sd    = sapply(perm.out, function(x) sd(x$perm)),
  # one-sided: P(perm <= real)  (real expected strongly negative)
  p.one.side = sapply(perm.out, function(x) mean(x$perm <= x$real))
)
cat("\nT3 summary:\n"); print(t3.pval)
write.csv(t3.pval, t3.out, row.names = FALSE)

# Plot: permutation null vs real, per ensemble
t3.plot.df <- bind_rows(lapply(names(perm.out), function(e) {
  data.frame(ensemble = e, perm.mean = perm.out[[e]]$perm)
}))
real.line <- data.frame(ensemble = names(perm.out),
                        real = sapply(perm.out, function(x) x$real))

p.t3 <- ggplot(t3.plot.df, aes(x = perm.mean)) +
  geom_histogram(bins = 40, fill = "gray70", color = "gray30") +
  geom_vline(data = real.line, aes(xintercept = real),
             color = "red", linewidth = 0.9) +
  facet_wrap(~ ensemble, scales = "free") +
  labs(title = "T3 — ENSO fingerprint permutation test",
       subtitle = "Null = mean Nino3.4-SOI corr. under shuffled SST-SLP member pairing; red = real pairing",
       x = "Mean Nino3.4 vs SOI correlation across members",
       y = "Shuffle count") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9))

ggsave(t3.fig, p.t3, width = 9, height = 4, dpi = 300)
cat("Saved:", t3.fig, "\n")

cat("\nVerification complete. Review T1-T3 outputs before adding member-level regression machinery below.\n")

# ============================================================
# PER-MEMBER AL SLP SD  vs  SST AR(1) CORRELATIONS
# ============================================================
# Mirrors the ERA5 analysis in analysis.R:
#   - SLP: area-weighted winter (Nov-Mar) mean over Aleutian Low box,
#          then 15-yr right-aligned rolling SD.
#   - SST: area-weighted winter mean -> 1950-1979 monthly climatology ->
#          detrended anomaly -> Nov-Mar annual mean ->
#          15-yr right-aligned rolling AR(1).
# Reuses the SST AR(1) cache written by CESM_ensemble_tests.R.
# Then correlates AL.sd against GOA.ar1 and EBS.ar1 per member, and
# compares the resulting r distributions across FCM vs MDM.

library(zoo)

# Aleutian-Low box (Litzow 2020 PNAS convention)
al.lon <- c(192.5, 207.5)
al.lat <- c(45,    55)

# Rolling window + overlap with ERA5
win.width <- 15
yr.min    <- 1964
yr.max    <- 2014

al.slp.cache <- "./Output/CESM_member_AL_SLP_SD.csv"
sst.ar1.cache <- "./Output/CESM_ensemble_sst_ar1.csv"   # written by CESM_ensemble_tests.R
corr.out     <- "./Output/CESM_member_correlations.csv"
corr.fig     <- "./Figures/CESM_member_ALSD_SSTAR1_r_distributions.png"

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
  rollapply(v, width = width, fill = NA, align = "right",
            FUN = function(x) if (all(is.na(x))) NA_real_ else sd(x, na.rm = TRUE))
}

# ---- AL SLP extraction ----
process_al_slp_member <- function(slp.file) {
  nc  <- nc_open(slp.file); on.exit(nc_close(nc), add = TRUE)
  lon <- read_coord(nc, "lon"); lat <- read_coord(nc, "lat")
  nt  <- nc$dim$time$len
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
  psl <- drop(psl)   # (lon, lat, time)
  stopifnot(length(dim(psl)) == 3)
  w <- cos(lat[i.lat] * pi / 180)
  W <- matrix(rep(w, each = length(i.lon)), length(i.lon), length(i.lat))
  al.month <- apply(psl, 3, function(sl) {
    m <- !is.na(sl); if (!any(m)) NA_real_ else sum(sl[m] * W[m]) / sum(W[m])
  })
  n.t <- length(al.month)
  years.m  <- rep(1850:2014, each = 12)[seq_len(n.t)]
  months.m <- rep(1:12,      times = 165)[seq_len(n.t)]

  winter <- data.frame(year = years.m, month = months.m, slp = al.month) %>%
    filter(month %in% c(11, 12, 1, 2, 3)) %>%
    mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
    group_by(win.year) %>%
    summarise(AL = mean(slp, na.rm = TRUE), n = n(), .groups = "drop") %>%
    filter(n == 5) %>% rename(year = win.year) %>% arrange(year) %>% select(-n)

  winter %>% mutate(AL.sd = roll_sd(AL)) %>%
    filter(year >= yr.min, year <= yr.max) %>%
    select(year, AL.sd)
}

# ---- Run SLP extraction with caching ----
run_slp_ensemble <- function(dir, model.type) {
  files <- list.files(dir, "\\.nc$", full.names = TRUE)
  cat(model.type, ":", length(files), "SLP files\n")
  lapply(files, function(f) {
    cat("  ", basename(f), "\n")
    out <- tryCatch(process_al_slp_member(f),
                    error = function(e) { cat("    ERROR:", e$message, "\n"); NULL })
    if (!is.null(out)) out %>% mutate(member = parse_id(basename(f), model.type),
                                      model  = model.type)
  }) %>% bind_rows()
}

if (file.exists(al.slp.cache)) {
  message("Loading cached AL SLP SD: ", al.slp.cache)
  al.df <- read.csv(al.slp.cache)
} else {
  message("Processing CESM2 SLP ensembles for AL box...")
  al.df <- bind_rows(
    run_slp_ensemble(fcm.slp.dir, "FCM"),
    run_slp_ensemble(mdm.slp.dir, "MDM")
  )
  write.csv(al.df, al.slp.cache, row.names = FALSE)
  message("Saved: ", al.slp.cache)
}

# ---- Load SST AR(1) from ensemble_tests cache ----
if (!file.exists(sst.ar1.cache))
  stop("SST AR(1) cache not found: ", sst.ar1.cache,
       "\nRun Scripts/CESM_ensemble_tests.R first.")
sst.df <- read.csv(sst.ar1.cache)

# CESM_ensemble_tests.R stored MDM IDs as full basename stems
# (e.g. "cesm2.MD.postP.01.SST.185001-201412"), while this script's
# parse_id emits "01". Normalise both sides to a common ID.
normalise_ids <- function(df) {
  df$member <- mapply(function(m, e) parse_id(m, e), df$member, df$model)
  df
}
sst.df <- normalise_ids(sst.df)
al.df  <- normalise_ids(al.df)

# ---- Merge and correlate per member ----
merged <- inner_join(
  sst.df %>% select(model, member, year, GOA.ar1, EBS.ar1),
  al.df  %>% select(model, member, year, AL.sd),
  by = c("model", "member", "year")
)
cat("Merged rows: ", nrow(merged),
    "  FCM members: ", length(unique(merged$member[merged$model == "FCM"])),
    "  MDM members: ", length(unique(merged$member[merged$model == "MDM"])), "\n")

corr.df <- merged %>%
  group_by(model, member) %>%
  summarise(GOA = safe_cor(AL.sd, GOA.ar1),
            EBS = safe_cor(AL.sd, EBS.ar1),
            .groups = "drop") %>%
  pivot_longer(c(GOA, EBS), names_to = "region", values_to = "r") %>%
  mutate(region = factor(region, levels = c("GOA", "EBS")),
         model  = factor(model,  levels = c("FCM", "MDM")))

cat("\nPer-member r summary:\n")
print(corr.df %>% group_by(model, region) %>%
        summarise(n        = n(),
                  n.na     = sum(is.na(r)),
                  mean.r   = mean(r, na.rm = TRUE),
                  median.r = median(r, na.rm = TRUE),
                  min.r    = min(r, na.rm = TRUE),
                  max.r    = max(r, na.rm = TRUE),
                  .groups  = "drop"))
write.csv(corr.df, corr.out, row.names = FALSE)
cat("Saved:", corr.out, "\n")

# ---- 4-panel r distribution plot ----
ens.mean <- corr.df %>% group_by(model, region) %>%
  summarise(r.bar = mean(r, na.rm = TRUE), .groups = "drop")

p.corr <- ggplot(corr.df, aes(x = r, fill = model)) +
  geom_histogram(bins = 20, color = "gray30", alpha = 0.85,
                 position = "identity") +
  geom_vline(data = ens.mean, aes(xintercept = r.bar),
             color = "black", linewidth = 0.8, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.3) +
  facet_wrap(region ~ model, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("FCM" = "darkorange", "MDM" = "steelblue")) +
  labs(title    = "Per-member correlation: AL SLP SD vs SST AR(1)",
       subtitle = paste0("15-yr right-aligned rolling windows, winter (Nov-Mar), ",
                         yr.min, "-", yr.max,
                         ". Dashed = ensemble mean r."),
       x = "r  (AL SLP SD vs SST AR(1))", y = "Number of ensemble members") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "none",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        plot.subtitle    = element_text(size = 9))

ggsave(corr.fig, p.corr, width = 9, height = 6, dpi = 300)
cat("Saved:", corr.fig, "\n")

# ============================================================
# BASELINE AL SLP -> SST REGRESSION (annual, no rolling windows)
# ============================================================
# Sanity check that FCM preserves the expected AL-box SLP -> GOA/EBS SST
# relationship and that MDM does not. Per member, regress winter SST anomaly
# on winter AL SLP anomaly (both detrended) with GLS + corAR1() residuals.
# All years available (1850-2014). Four-panel histogram of slopes.

library(nlme)
library(sf)
library(zoo)

# GOA/EBS polygons (match analysis.R / CESM_ensemble_tests.R)
goa.x.360 <- c(191, 191, 203, 205, 208, 223, 234)
goa.y     <- c(50,  53,  57.5, 59,  61,  61,  50)
ebs.x.360 <- c(183, 183, 203, 203, 191)
ebs.y     <- c(53,  65,  65,  57.5, 53)
goa.x <- ifelse(goa.x.360 > 180, goa.x.360 - 360, goa.x.360)
ebs.x <- ifelse(ebs.x.360 > 180, ebs.x.360 - 360, ebs.x.360)

clim.yr.min <- 1950   # climatology reference for monthly anomalies
clim.yr.max <- 1979

sst.anom.cache  <- "./Output/CESM_member_winter_SST_anom.csv"
al.anom.cache   <- "./Output/CESM_member_winter_AL_anom.csv"
slope.out       <- "./Output/CESM_member_AL_SST_slopes.csv"
slope.fig       <- "./Figures/CESM_member_AL_SST_slope_distributions.png"

# ---- Pipeline: monthly area-mean -> climatology -> anomaly -> detrend -> winter mean ----
monthly_winter_anom <- function(monthly.vec, years.m, months.m,
                                clim.min = clim.yr.min, clim.max = clim.yr.max) {
  df <- data.frame(year = years.m, month = months.m, x = monthly.vec)
  clim <- df %>% filter(year >= clim.min, year <= clim.max) %>%
    group_by(month) %>% summarise(c = mean(x, na.rm = TRUE), .groups = "drop")
  df <- df %>% left_join(clim, by = "month") %>%
    mutate(anom = x - c)
  # detrend monthly anomaly
  ok <- is.finite(df$anom)
  if (sum(ok) < 3) return(data.frame(year = integer(), anom = numeric()))
  df$anom.d <- NA_real_
  df$anom.d[ok] <- residuals(lm(df$anom[ok] ~ seq_len(sum(ok))))
  winter <- df %>% filter(month %in% c(11, 12, 1, 2, 3)) %>%
    mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
    group_by(win.year) %>%
    summarise(anom = mean(anom.d, na.rm = TRUE), n = n(), .groups = "drop") %>%
    filter(n == 5) %>% rename(year = win.year) %>% select(year, anom) %>%
    arrange(year)
  winter
}

# ---- SST winter anomaly per member (GOA + EBS) ----
process_sst_winter_anom <- function(sst.file) {
  nc  <- nc_open(sst.file); on.exit(nc_close(nc), add = TRUE)
  lon <- read_coord(nc, "lon"); lat <- read_coord(nc, "lat")
  sc  <- build_start_count(nc, "SST", list(z_t = list(start = 1, count = 1)))
  sst <- ncvar_get(nc, "SST", start = sc$start, count = sc$count,
                   collapse_degen = FALSE)
  sst[sst > 1e20] <- NA_real_
  sst <- drop(sst)                   # (lon, lat, time)
  stopifnot(length(dim(sst)) == 3)
  n.t <- dim(sst)[3]
  years.m  <- rep(1850:2014, each = 12)[seq_len(n.t)]
  months.m <- rep(1:12,      times = 165)[seq_len(n.t)]

  goa.coords <- rbind(cbind(goa.x, goa.y), c(goa.x[1], goa.y[1]))
  ebs.coords <- rbind(cbind(ebs.x, ebs.y), c(ebs.x[1], ebs.y[1]))
  goa.sf <- st_sf(geometry = st_sfc(st_polygon(list(goa.coords)), crs = 4326))
  ebs.sf <- st_sf(geometry = st_sfc(st_polygon(list(ebs.coords)), crs = 4326))

  lons.sf <- ifelse(lon > 180, lon - 360, lon)
  grid    <- expand.grid(lon = lons.sf, lat = lat)
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

  goa.w <- monthly_winter_anom(wmean(in.goa), years.m, months.m) %>%
    rename(GOA.anom = anom)
  ebs.w <- monthly_winter_anom(wmean(in.ebs), years.m, months.m) %>%
    rename(EBS.anom = anom)
  full_join(goa.w, ebs.w, by = "year") %>% arrange(year)
}

# ---- AL SLP winter anomaly per member ----
process_al_winter_anom <- function(slp.file) {
  nc  <- nc_open(slp.file); on.exit(nc_close(nc), add = TRUE)
  lon <- read_coord(nc, "lon"); lat <- read_coord(nc, "lat")
  i.lon <- which(lon >= al.lon[1] & lon <= al.lon[2])
  i.lat <- which(lat >= al.lat[1] & lat <= al.lat[2])
  stopifnot(length(i.lon) > 0, length(i.lat) > 0)
  sel <- list(lon = list(start = min(i.lon), count = length(i.lon)),
              lat = list(start = min(i.lat), count = length(i.lat)))
  sc  <- build_start_count(nc, "PSL", sel)
  psl <- ncvar_get(nc, "PSL", start = sc$start, count = sc$count,
                   collapse_degen = FALSE)
  psl <- drop(psl)
  stopifnot(length(dim(psl)) == 3)
  n.t <- dim(psl)[3]
  years.m  <- rep(1850:2014, each = 12)[seq_len(n.t)]
  months.m <- rep(1:12,      times = 165)[seq_len(n.t)]
  w <- cos(lat[i.lat] * pi / 180)
  W <- matrix(rep(w, each = length(i.lon)), length(i.lon), length(i.lat))
  al.month <- apply(psl, 3, function(sl) {
    m <- !is.na(sl); if (!any(m)) NA_real_ else sum(sl[m] * W[m]) / sum(W[m])
  })
  monthly_winter_anom(al.month, years.m, months.m) %>% rename(AL.anom = anom)
}

# ---- Ensemble runners ----
run_anom_ensemble <- function(dir, model.type, proc.fn) {
  files <- list.files(dir, "\\.nc$", full.names = TRUE)
  cat(model.type, ":", length(files), "files\n")
  lapply(files, function(f) {
    cat("  ", basename(f), "\n")
    out <- tryCatch(proc.fn(f),
                    error = function(e) { cat("    ERROR:", e$message, "\n"); NULL })
    if (!is.null(out)) out %>% mutate(member = parse_id(basename(f), model.type),
                                      model  = model.type)
  }) %>% bind_rows()
}

if (file.exists(sst.anom.cache)) {
  message("Loading cached SST winter anom: ", sst.anom.cache)
  sst.anom.df <- read.csv(sst.anom.cache)
} else {
  message("Processing per-member SST winter anomalies...")
  sst.anom.df <- bind_rows(
    run_anom_ensemble(fcm.sst.dir, "FCM", process_sst_winter_anom),
    run_anom_ensemble(mdm.sst.dir, "MDM", process_sst_winter_anom)
  )
  write.csv(sst.anom.df, sst.anom.cache, row.names = FALSE)
}

if (file.exists(al.anom.cache)) {
  message("Loading cached AL winter anom: ", al.anom.cache)
  al.anom.df <- read.csv(al.anom.cache)
} else {
  message("Processing per-member AL SLP winter anomalies...")
  al.anom.df <- bind_rows(
    run_anom_ensemble(fcm.slp.dir, "FCM", process_al_winter_anom),
    run_anom_ensemble(mdm.slp.dir, "MDM", process_al_winter_anom)
  )
  write.csv(al.anom.df, al.anom.cache, row.names = FALSE)
}

# ---- Merge and fit GLS per member per region ----
reg.df <- inner_join(
  sst.anom.df %>% select(model, member, year, GOA.anom, EBS.anom),
  al.anom.df  %>% select(model, member, year, AL.anom),
  by = c("model", "member", "year")
)
cat("Regression-data rows: ", nrow(reg.df),
    "  FCM members: ", length(unique(reg.df$member[reg.df$model == "FCM"])),
    "  MDM members: ", length(unique(reg.df$member[reg.df$model == "MDM"])), "\n")

gls_slope <- function(y, x) {
  ok <- is.finite(y) & is.finite(x)
  if (sum(ok) < 20) return(NA_real_)
  d <- data.frame(y = y[ok], x = x[ok], t = seq_len(sum(ok)))
  fit <- tryCatch(
    gls(y ~ x, data = d, correlation = corAR1(form = ~ t),
        method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  unname(coef(fit)[["x"]])
}

groups <- reg.df %>% group_by(model, member) %>% group_split()
slope.df <- lapply(groups, function(d) {
  cat("  ", d$model[1], "-", d$member[1], "\n")
  data.frame(
    model  = d$model[1],
    member = d$member[1],
    GOA    = gls_slope(d$GOA.anom, d$AL.anom),
    EBS    = gls_slope(d$EBS.anom, d$AL.anom)
  )
}) %>% bind_rows() %>%
  pivot_longer(c(GOA, EBS), names_to = "region", values_to = "slope") %>%
  mutate(region = factor(region, levels = c("GOA", "EBS")),
         model  = factor(model,  levels = c("FCM", "MDM")))

cat("\nPer-member slope summary (units: SST anom per Pa of AL SLP anom):\n")
print(slope.df %>% group_by(model, region) %>%
        summarise(n           = n(),
                  n.na        = sum(is.na(slope)),
                  mean.slope  = mean(slope, na.rm = TRUE),
                  median.slope= median(slope, na.rm = TRUE),
                  p.neg       = mean(slope < 0, na.rm = TRUE),
                  .groups     = "drop"))
write.csv(slope.df, slope.out, row.names = FALSE)
cat("Saved:", slope.out, "\n")

# ---- 4-panel slope distribution plot ----
ens.mean.slope <- slope.df %>% group_by(model, region) %>%
  summarise(s.bar = mean(slope, na.rm = TRUE), .groups = "drop")

p.slope <- ggplot(slope.df, aes(x = slope, fill = model)) +
  geom_histogram(bins = 20, color = "gray30", alpha = 0.85,
                 position = "identity") +
  geom_vline(data = ens.mean.slope, aes(xintercept = s.bar),
             color = "black", linewidth = 0.8, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray50", linewidth = 0.3) +
  facet_grid(region ~ model, scales = "free") +
  scale_fill_manual(values = c("FCM" = "darkorange", "MDM" = "steelblue")) +
  labs(title    = "Per-member GLS slope: winter SST anomaly on winter AL SLP anomaly",
       subtitle = "Annual winter means (Nov-Mar), 1850-2014, detrended. GLS with corAR1(). Dashed = ensemble mean.",
       x = "Slope  (deg C per Pa)", y = "Number of ensemble members") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "none",
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        plot.subtitle    = element_text(size = 9))

ggsave(slope.fig, p.slope, width = 9, height = 6, dpi = 300)
cat("Saved:", slope.fig, "\n")
