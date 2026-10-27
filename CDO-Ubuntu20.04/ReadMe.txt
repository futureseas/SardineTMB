===============================================================================
DATASET README
===============================================================================

TITLE: 
ROMS mixed layer depth (MLD) hindcast output for the California Current 
Ecosystem (2008–2020)

DATASET AUTHORS:
1) Stephanie C. Hopkins, University of California, Santa Cruz, Institute of 
Marine Sciences’ Fisheries Collaborative Program, 1156 High Street, Santa Cruz, 
California 95064, U.S.A. (ORCID: 0000-0002-9910-1148)
2) Mercedes Pozo Buil, University of South Carolina, School of the Earth, Ocean 
and Environment, 701 Sumter St., Columbia, SC 29208, U.S.A.
3) Jerome Fiechter, University of California, Santa Cruz, Institute of Marine 
Sciences, Ocean Sciences Department, 1156 High Street, Santa Cruz, California 
95064, U.S.A.

ASSOCIATED MANUSCRIPT:
Stephanie C. Hopkins, Desiree Tommasi, Mercedes Pozo Buil, Isaac C. Kaplan, 
Robert P. Wildermuth, Barbara Muhling, Jennifer L. Boldt, Alexander J. Jensen, 
Caitlin Allen Akselrud (2026). Harnessing transboundary survey integration for 
species distribution models: technical maneuvers for small pelagic fish. 
Canadian Journal of Fisheries and Aquatic Sciences. [In review]

ASSOCIATED CODE AND PROJECT PIPELINE:
Hopkins, S. C., & Tommasi, D. (2026). R Code & Simulation Pipeline corresponding 
to "Harnessing transboundary survey integration for species distribution models: 
technical maneuvers for small pelagic fish" (CJFAS publication). 
GitHub Repository: https://github.com/futureseas/SardineTMB	


LICENSING & TERMS OF USE:
- Dataset (NetCDF files): Creative Commons Attribution 4.0 International (CC BY 4.0)
  https://creativecommons.org/licenses/by/4.0/
- Code & Simulation Pipeline (GitHub): GNU General Public License v3.0 (GPL-3.0)
  https://www.gnu.org/licenses/gpl-3.0.html
- Manuscript & Text: Open Access subject to Canadian Journal of Fisheries and 
  Aquatic Sciences (CJFAS) publisher terms and copyright.

CITING THIS DATASET:
If you use these data, please cite both the dataset DOI and the associated 
manuscript:
  - Dataset Citation: Stephanie C. Hopkins, Mercedes Pozo Buil, and 
    Jerome Fiechter (2026). ROMS Mixed Layer Depth (MLD) Hindcast Output for the 
    California Current Ecosystem (2008–2020). SEANOE. 
    https://doi.org/10.17882/120759
  - Related Paper: Hopkins et al. (2026). Harnessing transboundary survey 
    integration for species distribution models: technical maneuvers for small 
    pelagic fish. Can. J. Fish. Aquat. Sci. [In review]
  - Code and Data Pipeline Citation: Hopkins, S. C., & Tommasi, D. (2026). R Code & 
    Simulation Pipeline corresponding to "Harnessing transboundary survey integration 
    for species distribution models: technical maneuvers for small pelagic fish" 
    (CJFAS publication). GitHub Repository: https://github.com/futureseas/SardineTMB

===============================================================================
1. OVERVIEW & SCOPE
===============================================================================
This dataset provides monthly mixed layer depth (MLD) hindcast outputs derived 
from a Regional Ocean Modeling System (ROMS) simulation of the California 
Current Ecosystem (2008–2020). 

The data are formatted to serve as direct inputs for fisheries spatial modeling. 
Specifically, each file contains elevation-masked MLD values along with appended 
dynamic bounding box parameters (`nbnds`) derived from model outputs. These 
bounding box definitions allow users to regrid the raw curvilinear model grid 
onto a standardized 10 x 10 km grid suitable for species distribution modeling 
(SDM) of Pacific Sardine (Sardinopsis sagax).

- Temporal Range: 2008-01-01 to 2020-12-31 (Monthly time steps)
- Spatial Range: California Current Ecosystem (WGS84)
- Nominal Average Spacing: ~0.106° Lon × ~0.075° Lat (~10 km resolution)
- Grid Features: Includes dynamic bounding box corner coordinates (`nbnds`) for 
  curvilinear-to-rectilinear regridding

===============================================================================
2. FILE INVENTORY & NAMING CONVENTIONS
===============================================================================
The dataset consists of 13 NetCDF (.nc) files (one per calendar year):

Filename Pattern: mld_hindcast_YYYY_masked.nc
  - Example: mld_hindcast_2020_masked.nc
  - Content: Monthly Mixed Layer Depth computed using potential density with a 
    fixed threshold criterion (0.03 kg/m³ difference from surface) with applied 
    elevation mask and dynamic bounding box (`nbnds`) definitions.

File Sizes: ~11 MB per file (~140 MB total)

===============================================================================
3. EMBEDDED METADATA & DATA STRUCTURE
===============================================================================
All files are saved in NetCDF-4 format adhering to Climate and Forecast (CF) 
Metadata Conventions. Full variable descriptions, standard names, units, 
missing/fill values (`NaNf`), and coordinate systems are embedded in each file 
header.

To inspect the embedded metadata, use:
  - Command Line (netCDF tools):  ncdump -h mld_hindcast_2020_masked.nc
  - Python (xarray): import xarray as xr; 
    ds = xr.open_dataset('mld_hindcast_2020_masked.nc'); print(ds)
  - R (ncdf4): library(ncdf4); nc <- nc_open('mld_hindcast_2020_masked.nc'); 
    print(nc)

Key Variables Included in File Headers:
  - Coordinate dimensions: xi_rho, eta_rho, time, nbnds
  - Variables: mld (Mixed Layer Depth, m), lon_rho / lat_rho (Grid center 
    coordinates), lon_bnds / lat_bnds (Dynamic cell bounding box corners)

===============================================================================
4. MODEL & PROVENANCE INFORMATION
===============================================================================
- Base Simulation: ROMS Northeast Pacific (ROMS-NEP) regional ocean model runs, 
  generated by Mercedes Pozo Buil.
- Model Code & Configuration Notice: While the underlying ROMS source code is 
  openly available, the specific numerical configuration, parameterizations, and 
  full raw simulation output remain proprietary to the original modelers. This 
  dataset contains processed, domain-masked MLD subset outputs with dynamic cell 
  boundaries designed specifically as input for downstream fisheries applications.
- Regridding & Processing Code: Open-source R scripts demonstrating how to 
  extract, regrid (10 x 10 km grid), and translate these dataset inputs into 
  species distribution models are archived under a GPL-3.0 license at:
  https://github.com/futureseas/SardineTMB.git. Following these procedures, the 
  spatial extent of the dataset spans roughly 162.45°W to 109.91°W and 18.10°N 
  to 66.53°N.

===============================================================================
5. CONTACT INFORMATION
===============================================================================
For questions regarding dataset generation, grid structure, or usage:
  - Contact Person: Stephanie C. Hopkins
  - Email: shopkin1@ucsc.edu
  - Institutions: University of California Santa Cruz / NOAA SWFSC / NOAA NWFSC