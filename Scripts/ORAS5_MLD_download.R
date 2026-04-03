# PURPOSE: Download ORAS5 monthly mean mixed layer depth
# Author: Mike Litzow
#
# Source: Copernicus Climate Data Store (CDS)
#   Dataset:  reanalysis-oras5
#   Variable: mixed_layer_depth_0_03 (monthly mean, m)
#   Coverage: 1958-2025 (monthly, global)
#
# CDS delivers files as ZIP archives of monthly NetCDF files. This script:
#   1. Downloads one year at a time as ZIP
#   2. Unzips the 12 monthly NetCDF files
#   3. Uses CDO to merge months and crop to North Pacific
#   4. Merges all years into one NetCDF
#
# Prerequisites (run once from terminal):
#   pip install cdsapi
#   brew install cdo
#   Create ~/.cdsapirc with your CDS API key.

py_path <- "/Users/MikeLitzow/Documents/github/Reddening/.venv/bin/python"
cdo_path <- "/usr/local/bin/cdo"

if (!file.exists(path.expand("~/.cdsapirc"))) {
  stop("~/.cdsapirc not found. Create it with your CDS API key.")
}

py_script <- sprintf('
import cdsapi, os, zipfile, subprocess, glob
from pathlib import Path

tmp_dir = Path("./Data/tmp_oras5")
tmp_dir.mkdir(parents=True, exist_ok=True)

cdo = "%s"
c = cdsapi.Client()

# Crop to 20-66N latitude band only.
# ORAS5 uses an irregular ORCA grid so antimeridian lon splits are unreliable.
# Longitude subsetting (110-250E) is done in R analysis scripts instead.
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

    # Unzip 12 monthly NetCDF files
    print(f"{year}: unzipping...")
    with zipfile.ZipFile(zip_file, "r") as z:
        z.extractall(tmp_dir)
    monthly_files = sorted(glob.glob(str(tmp_dir / f"*_{year}??_*.nc")))
    if not monthly_files:
        raise RuntimeError(f"No NC files found after unzipping {zip_file}")

    # Merge 12 months into one yearly file, then crop latitude
    merged_glo = str(tmp_dir / f"tmp_global_{year}.nc")
    subprocess.run([cdo, "mergetime"] + monthly_files + [merged_glo], check=True)
    subprocess.run([cdo, LAT_CROP, merged_glo, str(nc_file)], check=True)

    # Clean up intermediates
    for f in monthly_files + [merged_glo]:
        os.remove(f)
    nc_files.append(str(nc_file))
    print(f"{year}: done.")

print("Merging all years...")
out = "./Data/oras5_mld_NP_1958_2025.nc"
subprocess.run([cdo, "cat"] + nc_files + [out], check=True)
print(f"Done: {out}")
', cdo_path)

py_file <- tempfile(fileext = ".py")
writeLines(py_script, py_file)
system(paste(py_path, py_file))
