# SardineTMB

> **Peer Review Snapshot Note (`v0.2.0-alpha`):** This repository contains the R code, Linux CDO re-gridding routines, and spatiotemporal modeling pipeline accompanying the manuscript submission to the *Canadian Journal of Fisheries and Aquatic Sciences* (CJFAS). Code refinement, efficiency tuning, and script modularization remain ongoing on the `main` branch.

---

## Overview

This repository provides a reproducible workflow for standardizing, integrating, and modeling transboundary surface-trawl survey data for Pacific sardine across the California Current Ecosystem (CCE; 28.65°N–54.80°N, 134.21°W–114.79°W).

Haul-level biomass from U.S. and Canadian surveys (2008, 2012–2019) is standardized by swept-volume ($\text{m}^3$) and partitioned into functional maturity stages:

1. **Immature:** $< 15.0\text{ cm}$ ($< L_{50}$ threshold).
2. **Maturing:** $15.0\text{ cm}$ to $< 20.0\text{ cm}$ (young, partially mature fish, ages 1–2).
3. **Fully Mature:** $\ge 20.0\text{ cm}$ ($L_{100}$ threshold).

Spatial haul observations are coupled with $10\text{ km} \times 10\text{ km}$ Regional Ocean Modeling System (ROMS) hindcasts to evaluate catchability effects associated with Mixed Layer Depth (MLD). Using `sdmTMB`, spatiotemporal candidate models evaluate mesh geometries (standard, anisotropy, SPDE barrier), Tweedie/delta-hurdle distribution families, error structures, and spatial covariates via 5-fold multifactorial stratified cross-validation. Selected models are projected onto a standardized grid to generate regional abundance indices, Center of Gravity (COG), and Effective Volume Occupied (EVO).

---

## Manuscripts & Data Citations

### Primary Manuscript

> Hopkins, S. C., Tommasi, D., Pozo Buil, M., Kaplan, I. C., Wildermuth, R. P., Muhling, B., Boldt, J. L., Jensen, A. J., & Allen Akselrud, C. (2026). *Harnessing transboundary survey integration for species distribution models: technical maneuvers for small pelagic fish*. Canadian Journal of Fisheries and Aquatic Sciences [Under Review].

### Repository Citation

> Hopkins, S. C., & Tommasi, D. (2026). *R Code & Simulation Pipeline corresponding to "Harnessing transboundary survey integration for species distribution models: technical maneuvers for small pelagic fish"*. GitHub Repository: [https://github.com/futureseas/SardineTMB](https://github.com/futureseas/SardineTMB)

### Primary Data Sources

* **ROMS Hindcast MLD Layers:** Stephanie, C., Pozo Buil, M., & Fiechter, J. (2026). *ROMS mixed layer depth (MLD) hindcast output for the California Current Ecosystem (2008–2020)*. SEANOE. [https://doi.org/10.17882/120759](https://doi.org/10.17882/120759) (CC-BY 4.0).
* **U.S. ATM CPS Trawl Data:** NOAA Fisheries InPort Item 20693. [https://www.fisheries.noaa.gov/inport/item/20693](https://www.fisheries.noaa.gov/inport/item/20693)
* **Canadian IPES Trawl Data:** Fisheries and Oceans Canada (DFO) Open Data Portal. [Open Government Portal Link](https://open.canada.ca/data/en/dataset/aa13418a-12f3-40d1-bfa8-83c7b8c795b0) *(Night pelagic trawl survey data available upon request from DFO)*.
* **Coastline Geometry (Land Barrier Mesh):** Global Self-consistent, Hierarchical, High-resolution Geography Database (GSHHG v2.3.7). National Oceanic and Atmospheric Administration (NOAA). https://www.soest.hawaii.edu/wessel/gshhg/*.

---

## System Requirements & Dependencies

The complete pipeline spans R statistical computing environments and Linux bash environments (for CDO re-gridding).

### Hardware & Environment Testing

* **Operating Systems:** Ubuntu 22.04 / 24.04 LTS (WSL2 / Native Linux) and Windows 11.
* **R Environment:** `R (>= 4.5.2)`
* **Climate Data Operators (CDO):** `CDO v2.4.0` with `ecCodes v2.34.1` and `NetCDF v4.9.2` (located in `CDO-Ubuntu20.04/`).

### R Dependencies

```R
# Required R Libraries
install.packages(c(
  "sdmTMB", "TMB", "sf", "stars", 
  "tidyverse", "data.table", "DHARMa", 
  "Matrix", "ggplot2"
))

```

---

## Data Setup Before Running

1. HTML APIs have been enbedded in the code for direct extraction of NOAA data from https://www.fisheries.noaa.gov and requires no further action (step 1)
2. For Canadian data, you must both contact the DFO directly for access to the night trawl surveys dataset and download the IPES dataset from the Open Government Portal Link
   (https://open.canada.ca/data/en/dataset/aa13418a-12f3-40d1-bfa8-83c7b8c795b0) before running Step 4.
3. For Step 5, you must download the pre-processed ROMS MLD NetCDF files from SEANOE (https://doi.org/10.17882/120759).
4. Place the downloaded `.nc` files into the local directory `Data/ROMS2/Masked` before executing Step 6.
5. To create the barrier mesh in sdmTMB in Step 8, you will need the Coastline Geometry. The pre-processed coastline shapefile
   (GSHHG v2.3.7 cropped to 28.65°N–54.80°N, 134.21°W–114.79°W with a 100 m survey transect buffer) is provided in the `LandBarrier/` directory for immediate use in `sdmTMB` mesh construction.

---

---

## Repository Directory Tree

SardineTMB/
├── CDO-Ubuntu20.04/                 # Linux CDO shell scripts & ROMS re-gridding routines
├── Final.Model.Scripts/             # Executable Rmd scripts for baseline & length-bin models
├── LandBarrier/                     # Pre-processed land barrier shapefile (GSHHG full-res coastline cropped to study domain + 100m survey buffer)
├── Final.Figures/                   # [Generated by code] Publication figures and diagnostic plots
├── Final.Results/                   # [Generated by code] Decision tables, fitted model objects (.rds), and CIs
├── LICENSE                          # GNU General Public License v3.0
└── README.md                        # Repository documentation

---

---

## Analytical Workflow & Pipeline Execution

The pipeline is organized sequentially across 14 computational steps:

| Step | Objective | Primary Script | Core Output / Deliverable |
| --- | --- | --- | --- |
| **1** | Resolve U.S. survey precision, missing specimen entries, and sample quality | `ATM Survey Sardine Data Prep.First.Filter.and.Specimen.Data.Rmd` | Combined haul and specimen composition data (`2003.2024.ATM.Survey.with.specimen.data.csv`) |
| **2** | Allocate missing weight/length via spatio-temporal nearest-haul matrices | `ATM Survey Sardine Data Prep.No.Specimen.Data.Rmd` | Integrated specimen compositions (`2003.2024.ATM.Survey.with.integrated.specimen.data.csv`) |
| **3** | Calculate swept-volume effort, construct transects, and derive U.S. density | `ATM Survey Sardine Transect with Sweep Area Prep.Rmd` | Binned U.S. density csv(s) (`area.transects.final.*.csv`) & biomass audits |
| **4** | Clean, verify, and standardize Canadian IPES/DFO survey gear openings and catch | `IPES.and.Canadian.Data.V3.Rmd` | Standardized Canadian survey density data (`Canadian.set.1999.2024.csv`) |
| **5** | Merge U.S. & Canadian survey datasets across spatial coordinate projections | `Join.Canadian.American.Data.revised.Rmd` | Integrated transboundary summer dataset (`Canada.and.CPS.ATM.Survey.Data.w.latBin.csv`) |
| **6** | Mask and re-grid ROMS 3D hindcast fields to $10\text{ km} \times 10\text{ km}$ cells via Linux CDO | `CDO-Ubuntu20.04/CJFAS.MLD.work.R` | Masked & conservative/bilinear re-gridded NetCDF layers (`*_masked_regrid.nc`) |
| **7** | Intersect ROMS MLD and bathymetry grids with survey transect centroids | `test_regrid.function.plots.Rmd` | Transect dataset with environmental covariates (`mld.final_table_allHistoric.csv`) |
| **8** | Model selection routine for full dataset via 5-fold stratified cross-validation | `Integrated.Canadian.dat.cv.approach.to.model.selection.Rmd` | Model ranking metrics (`Model.decision.table.csv`) & top candidate model (`S2.fit.basic.model.rds`) |
| **9** | Spatiotemporal model selection across maturity length-bin candidate structures | `Integrated.Canadian.dat.w.LenBin.cv.approach.to.model.selection.V2.Rmd` | Top length-bin model fit with reml correction (`S2.fit.200plusLenBin.model_SptlVarying.w.sptl.reml.rds`) & cross-validation metrics |
| **11** | Optimized execution pipeline for global baseline model | `Final.Global.Model.Rmd` | Global model diagnostics, DHARMa residual maps, and parameter confidence intervals |
| **12** | Execution pipeline for final maturity length-bin model | `Final.LengthBin.Model.TEST.deltalognormal.2.Rmd` | Length-bin model diagnostics, DHARMa residuals, and spatial random fields |
| **13** | Derive regional/global abundance indices, COG, and EVO metrics | `Final.Model.Scripts/S2.fit.basic.model_10km2.grid.Rmd` | Projected abundance time-series and spatial distribution indices |
| **14** | Map spatially varying MLD and length-bin density fields and derive COG and EVO metrics | `Final.Model.Scripts/S2.DeltaLognormalLength.bins.model_10km2.grid.Rmd` | High-resolution publication figures for spatial density and management sub-areas |

---

## Auxiliary Tooling & Verification Statement

Gemini 3.6 Flash (Google LLC) was utilized as an auxiliary technical tool for code proofing, syntax error debugging, script refactoring, and documentation formatting. All algorithms, mathematical logic, statistical model choices, and numerical outputs were independently authored, verified, and validated by the primary authors.

---

---

## Contact & Support

For questions or issues regarding script execution, please open an Issue on this repository or contact the corresponding author: Stephanie C. Hopkins (shopkins@ucsc.edu).

---

## License

This repository is distributed under the **GNU General Public License v3.0**. See the [LICENSE](https://www.google.com/search?q=LICENSE) file for complete details.
