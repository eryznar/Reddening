# PURPOSE: Download ORAS5 monthly mean mixed layer depth for the North Pacific domain
# Author: Mike Litzow
#
# ORAS5 MLD is available from the Copernicus Marine Service (CMEMS):
#   Dataset: global-reanalysis-phy-001-031
#   Variable: mlotst (ocean mixed layer thickness, m)
#   Coverage: 1958-present (monthly)
#
# Uses the official copernicusmarine Python CLI via system().
# Install once from terminal:  pip install copernicusmarine
#
# Log in once from terminal (stores credentials in ~/.copernicusmarine/):
#   copernicusmarine login
#
# North Pacific domain — same as ERA5 SST: 20-66N, 110-250E (250E = -110 in ±180)
# ORAS5 starts 1958; no data before that year.

cmd <- paste(
  "copernicusmarine subset",
  "--dataset-id global-reanalysis-phy-001-031-grepv2-monthly",
  "--variable mlotst",
  "--start-datetime 1958-01-01T00:00:00",
  "--end-datetime   2025-12-31T00:00:00",
  "--minimum-latitude  20",
  "--maximum-latitude  66",
  "--minimum-longitude  110",
  "--maximum-longitude -110",   # 250E in ±180 space
  "--minimum-depth 0",
  "--maximum-depth 1",
  "--output-filename oras5_mld_NP_1958_2025.nc",
  "--output-directory ./Data",
  "--force-download"
)

system(cmd)
