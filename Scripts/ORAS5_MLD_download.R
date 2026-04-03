# PURPOSE: Download ORAS5 monthly mean mixed layer depth
# Author: Mike Litzow
#
# Source: Copernicus Climate Data Store (CDS)
#   Dataset:  reanalysis-oras5
#   Variable: ocean_mixed_layer_thickness_defined_by_sigma_theta
#   Coverage: 1958-2024 (monthly, global)
#
# Output is global. Crop to North Pacific (20-66N, 110-250E) in analysis scripts.
# The North Pacific domain crosses the antimeridian, so spatial subsetting
# is not requested here — crop with ncdf4 or terra in R instead.
#
# Prerequisites (run once from terminal):
#   pip install cdsapi
#   Then create ~/.cdsapirc using your CDS profile page:
#     url: https://cds.climate.copernicus.eu/api
#     key: YOUR-API-KEY

py_path <- "/Users/MikeLitzow/Documents/github/Reddening/.venv/bin/python"

# Check for CDS credentials
if (!file.exists(path.expand("~/.cdsapirc"))) {
  stop("~/.cdsapirc not found. Create it with your CDS API key before running.")
}

py_script <- '
import cdsapi, os, xarray as xr
from pathlib import Path

tmp_dir = Path("./Data/tmp_oras5")
tmp_dir.mkdir(parents=True, exist_ok=True)

c = cdsapi.Client()

# Download year by year (global), crop to North Pacific, then merge.
# North Pacific domain: 20-66N, 110-250E (crosses antimeridian).
# CDS does not support spatial subsetting, so we crop after each download.

cropped = []

for year in range(1958, 2026):
    product_type = "consolidated" if year <= 2014 else "operational"
    tmp_file = tmp_dir / f"oras5_mld_{year}.nc"

    if not tmp_file.exists():
        print(f"Downloading {year}...")
        c.retrieve(
            "reanalysis-oras5",
            {
                "product_type": [product_type],
                "vertical_resolution": "single_level",
                "variable": "mixed_layer_depth_0_03",
                "year":  [str(year)],
                "month": [f"{m:02d}" for m in range(1, 13)],
            },
            str(tmp_file)
        )

    # Crop to North Pacific: 20-66N, 110E-180E and -180W to -110W (= 110-250E)
    with open(tmp_file, "rb") as f:
        magic = f.read(4)
    if magic == b"GRIB":
        engine = "cfgrib"
    elif magic[:3] == b"CDF":
        engine = "scipy"
    else:
        engine = "h5netcdf"
    ds = xr.open_dataset(tmp_file, engine=engine)
    lat_name = [c for c in ds.coords if "lat" in c.lower()][0]
    lon_name = [c for c in ds.coords if "lon" in c.lower()][0]
    ds = ds.sel({lat_name: slice(20, 66)})
    west = ds.sel({lon_name: slice(110, 180)})
    east = ds.sel({lon_name: slice(-180, -110)})
    ds_np = xr.concat([west, east], dim=lon_name)
    cropped.append(ds_np)
    ds.close()
    print(f"  {year} cropped.")

print("Merging all years...")
merged = xr.concat(cropped, dim="time")
merged.to_netcdf("./Data/oras5_mld_NP_1958_2025.nc")
print("Done: ./Data/oras5_mld_NP_1958_2025.nc")
'

py_file <- tempfile(fileext = ".py")
writeLines(py_script, py_file)
system(paste(py_path, py_file))
