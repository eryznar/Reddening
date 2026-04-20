# PURPOSE: Download all input data for the Reddening project
# Author: Mike Litzow
#
# Sections:
#   1. ERA5 SST        — 1950-2025, monthly, North Pacific (ecmwfr / CDS)
#   2. ERA5 SLP        — 1950-2025, monthly, North Pacific (ecmwfr / CDS)
#   3. ORAS5 MLD       — 1958-2025, monthly, North Pacific (cdsapi / CDS)
#   4. 2026 update     — Jan-Mar 2026 for ERA5 SST and SLP; appended to existing files
#                        ORAS5 MLD 2026 months appended via CDO
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

message("All downloads complete.")
message("2026 output files:")
message("  SST: ./Data/era5_sst_NP_1950_2026.nc")
message("  SLP: ./Data/era5_slp_NP_1950_2026.nc")
message("  MLD: ./Data/oras5_mld_NP_1958_2026.nc")
message("Update file paths in analysis scripts to use 2026 versions.")
