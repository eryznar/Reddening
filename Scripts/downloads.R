# PURPOSE: Download all input data for the Reddening project
# Author: Mike Litzow
#
# Sections:
#   1. ERA5 SST        — 1950-2025, monthly, North Pacific (ecmwfr / CDS)
#   2. ERA5 SLP        — 1950-2025, monthly, North Pacific (ecmwfr / CDS)
#   3. ORAS5 MLD       — 1958-2025, monthly, North Pacific (cdsapi / CDS)
#   4. 2026 update     — Jan-Mar 2026 for ERA5 SST and SLP; appended to existing files
#                        ORAS5 MLD 2026 months appended via CDO
#   5. ERA5 SST LMEs   — Jan 1950-Mar 2026, monthly, boreal LMEs (40-84N, 0.5 deg)
#   6. CMIP6 SST       — historical+SSP2-4.5 and piControl, 20 models, from
#                        github.com/mikelitzow/CMIP6-attribution
#
# Prerequisites:
#   install.packages("ecmwfr")
#   pip install cdsapi   (in .venv)
#   brew install cdo
#   Create ~/.cdsapirc with your CDS API key.
#   Run wf_set_key() once to store your CDS credentials in the macOS Keychain.
#
# Skip logic:
#   Each section checks whether its output file already exists and skips the
#   download if so. Delete the relevant output file to force a re-download.

library(ecmwfr)

cdo     <- "/usr/local/bin/cdo"
py_path <- "/Users/MikeLitzow/Documents/github/Reddening/.venv/bin/python"

if (!file.exists(path.expand("~/.cdsapirc"))) {
  stop("~/.cdsapirc not found. Create it with your CDS API key.")
}

options(ecmwfr.sleep = 120)   # 120s between CDS status checks

# ============================================================
# SECTION 1: ERA5 SST 1950-2025
# ============================================================
sst.file <- "./Data/era5_sst_NP_1950_2025.nc"

if (file.exists(sst.file)) {
  message("ERA5 SST 1950-2025 already exists — skipping.")
} else {
  message("Downloading ERA5 SST 1950-2025...")
  wf_request(
    request = list(
      dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
      product_type       = "monthly_averaged_reanalysis",
      variable           = "sea_surface_temperature",
      year               = as.character(1950:2025),
      month              = sprintf("%02d", 1:12),
      time               = "00:00",
      area               = c(66, 110, 20, -110),   # N, W, S, E (250E = -110W)
      data_format        = "netcdf",
      download_format    = "unarchived",
      target             = "era5_sst_NP_1950_2025.nc"
    ),
    transfer = TRUE,
    path     = "./Data"
  )
  message("Done: ", sst.file)
}

# ============================================================
# SECTION 2: ERA5 SLP 1950-2025
# ============================================================
slp.file <- "./Data/era5_slp_NP_1950_2025.nc"

if (file.exists(slp.file)) {
  message("ERA5 SLP 1950-2025 already exists — skipping.")
} else {
  message("Downloading ERA5 SLP 1950-2025...")
  wf_request(
    request = list(
      dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
      product_type       = "monthly_averaged_reanalysis",
      variable           = "mean_sea_level_pressure",
      year               = as.character(1950:2025),
      month              = sprintf("%02d", 1:12),
      time               = "00:00",
      area               = c(70, 120, 20, -110),   # N, W, S, E (250E = -110W)
      data_format        = "netcdf",
      download_format    = "unarchived",
      target             = "era5_slp_NP_1950_2025.nc"
    ),
    transfer = TRUE,
    path     = "./Data"
  )
  message("Done: ", slp.file)
}

# ============================================================
# SECTION 3: ORAS5 MLD 1958-2025
# ============================================================
mld.file <- "./Data/oras5_mld_NP_1958_2025.nc"

if (file.exists(mld.file)) {
  message("ORAS5 MLD 1958-2025 already exists — skipping.")
} else {
  message("Downloading ORAS5 MLD 1958-2025 via Python/cdsapi...")

  py_script_mld <- sprintf('
import cdsapi, os, zipfile, subprocess, glob
from pathlib import Path

tmp_dir = Path("./Data/tmp_oras5")
tmp_dir.mkdir(parents=True, exist_ok=True)

cdo = "%s"
c = cdsapi.Client()

LAT_CROP = "sellonlatbox,-180,180,20,66"

nc_files = []

for year in range(1958, 2026):
    nc_file = tmp_dir / f"oras5_mld_NP_{year}.nc"
    if nc_file.exists():
        print(f"{year}: already done, skipping.")
        nc_files.append(str(nc_file))
        continue

    product_type = "consolidated" if year <= 2014 else "operational"
    zip_file = tmp_dir / f"oras5_mld_{year}.zip"

    if not zip_file.exists():
        print(f"{year}: downloading...")
        c.retrieve(
            "reanalysis-oras5",
            {
                "product_type": [product_type],
                "vertical_resolution": "single_level",
                "variable": "mixed_layer_depth_0_03",
                "year":  [str(year)],
                "month": [f"{m:02d}" for m in range(1, 13)],
            },
            str(zip_file)
        )

    print(f"{year}: unzipping...")
    with zipfile.ZipFile(zip_file, "r") as z:
        z.extractall(tmp_dir)
    monthly_files = sorted(glob.glob(str(tmp_dir / f"*_{year}??_*.nc")))
    if not monthly_files:
        raise RuntimeError(f"No NC files found after unzipping {zip_file}")

    merged_glo = str(tmp_dir / f"tmp_global_{year}.nc")
    subprocess.run([cdo, "mergetime"] + monthly_files + [merged_glo], check=True)
    subprocess.run([cdo, LAT_CROP, merged_glo, str(nc_file)], check=True)

    for f in monthly_files + [merged_glo]:
        os.remove(f)
    nc_files.append(str(nc_file))
    print(f"{year}: done.")

print("Merging all years...")
out = "./Data/oras5_mld_NP_1958_2025.nc"
subprocess.run([cdo, "cat"] + nc_files + [out], check=True)
print(f"Done: {out}")
', cdo)

  py_file <- tempfile(fileext = ".py")
  writeLines(py_script_mld, py_file)
  system(paste(py_path, py_file))
  message("Done: ", mld.file)
}

# ============================================================
# SECTION 4: 2026 UPDATE — Jan-Mar 2026
# ============================================================
# Downloads Jan-Mar 2026 for ERA5 SST and SLP as separate files,
# then appends them to the existing 1950-2025 files using CDO mergetime.
# Also downloads Jan-Feb 2026 ORAS5 MLD (Mar 2026 may not yet be available)
# and appends to the existing MLD file.
#
# Output files (overwrite existing):
#   ./Data/era5_sst_NP_1950_2026.nc
#   ./Data/era5_slp_NP_1950_2026.nc
#   ./Data/oras5_mld_NP_1958_2026.nc
#
# Skip logic: if the 2026 output file already exists, skip that variable.

# ---- ERA5 SST 2026 ----
sst.2026.dl   <- "./Data/era5_sst_NP_2026_update.nc"
sst.2026.out  <- "./Data/era5_sst_NP_1950_2026.nc"

if (file.exists(sst.2026.out)) {
  message("ERA5 SST 2026 update already exists — skipping.")
} else {
  message("Downloading ERA5 SST Jan-Mar 2026...")
  wf_request(
    request = list(
      dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
      product_type       = "monthly_averaged_reanalysis",
      variable           = "sea_surface_temperature",
      year               = "2026",
      month              = c("01", "02", "03"),
      time               = "00:00",
      area               = c(66, 110, 20, -110),
      data_format        = "netcdf",
      download_format    = "unarchived",
      target             = "era5_sst_NP_2026_update.nc"
    ),
    transfer = TRUE,
    path     = "./Data"
  )
  message("Appending ERA5 SST 2026 to existing file via CDO...")
  system(paste(cdo, "mergetime", sst.file, sst.2026.dl, sst.2026.out))
  message("Done: ", sst.2026.out)
}

# ---- ERA5 SLP 2026 ----
slp.2026.dl   <- "./Data/era5_slp_NP_2026_update.nc"
slp.2026.out  <- "./Data/era5_slp_NP_1950_2026.nc"

if (file.exists(slp.2026.out)) {
  message("ERA5 SLP 2026 update already exists — skipping.")
} else {
  message("Downloading ERA5 SLP Jan-Mar 2026...")
  wf_request(
    request = list(
      dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
      product_type       = "monthly_averaged_reanalysis",
      variable           = "mean_sea_level_pressure",
      year               = "2026",
      month              = c("01", "02", "03"),
      time               = "00:00",
      area               = c(70, 120, 20, -110),
      data_format        = "netcdf",
      download_format    = "unarchived",
      target             = "era5_slp_NP_2026_update.nc"
    ),
    transfer = TRUE,
    path     = "./Data"
  )
  message("Appending ERA5 SLP 2026 to existing file via CDO...")
  system(paste(cdo, "mergetime", slp.file, slp.2026.dl, slp.2026.out))
  message("Done: ", slp.2026.out)
}

# ---- ORAS5 MLD 2026 ----
# CDS delivers as ZIP of monthly NetCDF files; download Jan-Mar 2026 only.
# Note: 2026 is "operational" product type.
mld.2026.out <- "./Data/oras5_mld_NP_1958_2026.nc"

if (file.exists(mld.2026.out)) {
  message("ORAS5 MLD 2026 update already exists — skipping.")
} else {
  message("Downloading ORAS5 MLD Jan-Mar 2026 via Python/cdsapi...")

  py_script_mld2026 <- sprintf('
import cdsapi, os, zipfile, subprocess, glob
from pathlib import Path

tmp_dir = Path("./Data/tmp_oras5")
tmp_dir.mkdir(parents=True, exist_ok=True)

cdo = "%s"
c = cdsapi.Client()

LAT_CROP = "sellonlatbox,-180,180,20,66"

# Download Jan-Mar 2026 as a single ZIP
year = 2026
months_2026 = ["01", "02", "03"]
zip_file = tmp_dir / "oras5_mld_2026_update.zip"
nc_file  = tmp_dir / "oras5_mld_NP_2026_update.nc"

if not nc_file.exists():
    if not zip_file.exists():
        print("Downloading ORAS5 MLD Jan-Mar 2026...")
        c.retrieve(
            "reanalysis-oras5",
            {
                "product_type": ["operational"],
                "vertical_resolution": "single_level",
                "variable": "mixed_layer_depth_0_03",
                "year":  [str(year)],
                "month": months_2026,
            },
            str(zip_file)
        )

    print("Unzipping...")
    with zipfile.ZipFile(zip_file, "r") as z:
        z.extractall(tmp_dir)
    monthly_files = sorted(glob.glob(str(tmp_dir / f"*_{year}??_*.nc")))
    if not monthly_files:
        raise RuntimeError(f"No NC files found after unzipping {zip_file}")

    merged_glo = str(tmp_dir / "tmp_global_2026_update.nc")
    subprocess.run([cdo, "mergetime"] + monthly_files + [merged_glo], check=True)
    subprocess.run([cdo, LAT_CROP, merged_glo, str(nc_file)], check=True)

    for f in monthly_files + [merged_glo]:
        os.remove(f)
    print(f"Done: {nc_file}")
else:
    print("ORAS5 2026 update file already processed, skipping download.")

# Append to existing 1958-2025 file
existing = "./Data/oras5_mld_NP_1958_2025.nc"
out      = "./Data/oras5_mld_NP_1958_2026.nc"
print("Merging 2026 update with existing MLD file...")
subprocess.run([cdo, "mergetime", existing, str(nc_file), out], check=True)
print(f"Done: {out}")
', cdo)

  py_file2 <- tempfile(fileext = ".py")
  writeLines(py_script_mld2026, py_file2)
  system(paste(py_path, py_file2))
  message("Done: ", mld.2026.out)
}

# ============================================================
# SECTION 5: ERA5 SST 1950-2026 — BOREAL LMEs (0.5 deg)
# ============================================================
# SST for the boreal LME comparison (ERA5_processing_LMEs.R,
# calculate.ar1.sd.LME.R). The LMEs compared with the EBS and GOA span
# 40.8-83N at all longitudes, well outside the North Pacific domain above,
# so this is a separate request:
#   - area 40-84N, 180W-179.5E (the east edge stops at 179.5 so the
#     dateline column is not returned twice)
#   - regridded by CDS to 0.5 deg with the `grid` key (native is 0.25 deg)
#   - Jan 1950 - Mar 2026, the same months as era5_sst_NP_1950_2026.nc
# Downloaded in two pieces, then joined here in R (ncdf4) rather than CDO.
# Needs only library(ecmwfr) and a key stored with wf_set_key() -- no CDO,
# Python or ~/.cdsapirc -- so this section can be run on its own.
#
# Output: ./Data/era5_sst_LME_1950_2026.nc (sst in K, NaN over land;
#         dims longitude, latitude, valid_time as in the NP SST file)
#
# Skip logic: skips entirely if the joined file exists; otherwise skips
# whichever piece has already been downloaded.

lme.sst.dl   <- "./Data/era5_sst_LME_1950_2025.nc"
lme.sst.2026 <- "./Data/era5_sst_LME_2026_update.nc"
lme.sst.out  <- "./Data/era5_sst_LME_1950_2026.nc"

if (file.exists(lme.sst.out)) {
  message("ERA5 SST for boreal LMEs already exists — skipping.")
} else {
  if (!file.exists(lme.sst.dl)) {
    message("Downloading ERA5 SST 1950-2025 for boreal LMEs (0.5 deg)...")
    wf_request(
      request = list(
        dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
        product_type       = "monthly_averaged_reanalysis",
        variable           = "sea_surface_temperature",
        year               = as.character(1950:2025),
        month              = sprintf("%02d", 1:12),
        time               = "00:00",
        area               = c(84, -180, 40, 179.5),   # N, W, S, E
        grid               = c("0.5", "0.5"),
        data_format        = "netcdf",
        download_format    = "unarchived",
        target             = "era5_sst_LME_1950_2025.nc"
      ),
      transfer = TRUE,
      path     = "./Data",
      time_out = 4 * 3600   # default is 1 h; allow for a long CDS queue
    )
  }

  if (!file.exists(lme.sst.2026)) {
    message("Downloading ERA5 SST Jan-Mar 2026 for boreal LMEs (0.5 deg)...")
    wf_request(
      request = list(
        dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
        product_type       = "monthly_averaged_reanalysis",
        variable           = "sea_surface_temperature",
        year               = "2026",
        month              = c("01", "02", "03"),
        time               = "00:00",
        area               = c(84, -180, 40, 179.5),
        grid               = c("0.5", "0.5"),
        data_format        = "netcdf",
        download_format    = "unarchived",
        target             = "era5_sst_LME_2026_update.nc"
      ),
      transfer = TRUE,
      path     = "./Data",
      time_out = 4 * 3600
    )
  }

  if (!all(file.exists(c(lme.sst.dl, lme.sst.2026)))) {
    stop("LME SST download incomplete (the CDS job may still be queued). ",
         "Rerun this section; finished pieces are skipped.")
  }

  # ---- Join the two pieces ----
  message("Joining 1950-2025 and Jan-Mar 2026 LME SST...")
  nc.a <- ncdf4::nc_open(lme.sst.dl)
  nc.b <- ncdf4::nc_open(lme.sst.2026)

  # as.vector() drops the 1-d array dim that ncvar_get() attaches, which
  # would otherwise make the all.equal() spacing check below always fail
  lon <- as.vector(ncdf4::ncvar_get(nc.a, "longitude"))
  lat <- as.vector(ncdf4::ncvar_get(nc.a, "latitude"))
  if (!identical(lon, as.vector(ncdf4::ncvar_get(nc.b, "longitude"))) ||
      !identical(lat, as.vector(ncdf4::ncvar_get(nc.b, "latitude")))) {
    stop("The 1950-2025 and 2026 LME SST files are on different grids.")
  }
  if (!isTRUE(all.equal(abs(diff(lat)), rep(0.5, length(lat) - 1)))) {
    warning("LME SST grid spacing is not 0.5 deg; CDS may have ignored the grid key.")
  }

  # Months must run on without gaps or overlap across the two pieces
  t.a   <- ncdf4::ncvar_get(nc.a, "valid_time")
  t.all <- c(t.a, ncdf4::ncvar_get(nc.b, "valid_time"))
  dates <- as.Date(as.POSIXct(t.all, origin = "1970-01-01", tz = "UTC"))
  if (any(diff(12 * as.integer(format(dates, "%Y")) +
               as.integer(format(dates, "%m"))) != 1)) {
    stop("LME SST months have gaps or overlaps across the two pieces.")
  }

  keep <- !duplicated(lon %% 360)   # drop a repeated dateline column, if any

  # Write under a temporary name so an interrupted join is not mistaken
  # for a finished file by the skip logic
  tmp.out <- paste0(lme.sst.out, ".part")
  dim.lon <- ncdf4::ncdim_def("longitude", "degrees_east", lon[keep])
  dim.lat <- ncdf4::ncdim_def("latitude", "degrees_north", lat)
  dim.t   <- ncdf4::ncdim_def("valid_time", "seconds since 1970-01-01", t.all,
                              calendar = "proleptic_gregorian")
  var.sst <- ncdf4::ncvar_def("sst", "K", list(dim.lon, dim.lat, dim.t),
                              missval = NaN, longname = "Sea surface temperature",
                              prec = "float", compression = 4,
                              chunksizes = c(sum(keep), length(lat), 1))
  nc.out <- ncdf4::nc_create(tmp.out, var.sst, force_v4 = TRUE)
  ncdf4::ncvar_put(nc.out, var.sst,
                   ncdf4::ncvar_get(nc.a, "sst", collapse_degen = FALSE)[keep, , , drop = FALSE],
                   start = c(1, 1, 1),
                   count = c(sum(keep), length(lat), length(t.a)))
  ncdf4::ncvar_put(nc.out, var.sst,
                   ncdf4::ncvar_get(nc.b, "sst", collapse_degen = FALSE)[keep, , , drop = FALSE],
                   start = c(1, 1, length(t.a) + 1),
                   count = c(sum(keep), length(lat), length(t.all) - length(t.a)))
  ncdf4::ncatt_put(nc.out, 0, "source",
                   paste("ERA5 monthly means regridded to 0.5 deg by CDS; joined from",
                         basename(lme.sst.dl), "and", basename(lme.sst.2026)))
  ncdf4::nc_close(nc.out)
  ncdf4::nc_close(nc.a)
  ncdf4::nc_close(nc.b)
  file.rename(tmp.out, lme.sst.out)

  message("Done: ", lme.sst.out, " (", length(t.all), " months, ",
          format(dates[1], "%b %Y"), " - ", format(dates[length(dates)], "%b %Y"), ")")
}

# ============================================================
# SECTION 6: CMIP6 SST — github.com/mikelitzow/CMIP6-attribution
# ============================================================
# Monthly SST (tos), regridded to 1 deg over the North Pacific
# (109.5-250.5E, 10.5-66.5N), 250 years per run:
#   - ssp245/<model>_245.nc : historical + SSP2-4.5, 1850-2099
#   - ssp585/<model>.nc     : piControl and historical + SSP5-8.5. Only the
#                             piControl experiment is used, as the matching
#                             preindustrial run for each SSP2-4.5 model.
# Only the 20 models with SSP2-4.5 runs are downloaded (~1.6 GB), plus the
# region-specific model weights. Plain HTTPS; no CDS key needed, so this
# section can be run on its own.
#
# Output: ./Data/CMIP6/ssp245/*.nc, ./Data/CMIP6/ssp585/*.nc,
#         ./Data/CMIP6/normalized_CMIP6_weights.csv
# Skip logic: each file is skipped if already present.

cmip6.repo   <- "https://raw.githubusercontent.com/mikelitzow/CMIP6-attribution/master/CMIP6"
cmip6.models <- c("ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CAMS-CSM1-0",
                  "CanESM5-CanOE", "CanESM5", "CESM2-WACCM", "CESM2", "CIESM",
                  "CMCC-CM2-SR5", "CMCC-ESM2", "GISS-E2-1-G", "HadGEM3-GC31-LL",
                  "INM-CM4-8", "INM-CM5-0", "MIROC-ES2L", "MIROC6",
                  "NorESM2-LM", "NorESM2-MM", "UKESM1-0-LL")

options(timeout = max(1800, getOption("timeout")))   # files are 25-56 MB
for (sub in c("ssp245", "ssp585"))
  dir.create(file.path("./Data/CMIP6", sub), recursive = TRUE, showWarnings = FALSE)

cmip6.files <- c(
  "summaries/normalized_CMIP6_weights.csv",
  paste0("CMIP6_outputs/1850-2099_runs/ssp245/", cmip6.models, "_245.nc"),
  paste0("CMIP6_outputs/1850-2099_runs/ssp585/", cmip6.models, ".nc"))

for (f in cmip6.files) {
  dest <- file.path("./Data/CMIP6", sub("^CMIP6_outputs/1850-2099_runs/|^summaries/", "", f))
  if (file.exists(dest)) next
  message("Downloading ", basename(dest), " ...")
  # Download under a temporary name so an interrupted file is not skipped later
  download.file(paste(cmip6.repo, f, sep = "/"), paste0(dest, ".part"),
                mode = "wb", quiet = TRUE)
  file.rename(paste0(dest, ".part"), dest)
}
message("CMIP6 files present: ",
        length(list.files("./Data/CMIP6", pattern = "\\.nc$", recursive = TRUE)), " .nc")

message("All downloads complete.")
message("2026 output files:")
message("  SST: ./Data/era5_sst_NP_1950_2026.nc")
message("  SLP: ./Data/era5_slp_NP_1950_2026.nc")
message("  MLD: ./Data/oras5_mld_NP_1958_2026.nc")
message("  SST (boreal LMEs): ./Data/era5_sst_LME_1950_2026.nc")
message("Update file paths in analysis scripts to use 2026 versions.")
