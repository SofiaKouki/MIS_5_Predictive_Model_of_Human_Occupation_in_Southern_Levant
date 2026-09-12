# MIS 5 Predictive Model of Human Occupation — Southern Levant

[![License: CC0 1.0](https://img.shields.io/badge/License-CC0%201.0-lightgrey.svg)](https://creativecommons.org/publicdomain/zero/1.0/)
[![R ≥ 4.2](https://img.shields.io/badge/R-%E2%89%A5%204.2-276DC3?logo=r)](https://www.r-project.org/)

---

## Overview

This repository contains the R pipeline for a **predictive model of MIS 5 human occupation** across the Southern Levant, as described in:

> Samawi, O., Kouki, S., Beller, J. A., Hallinan, E., Rose, J. I., Bicho, N., Nassr, A., Collard, M., & Al-Nahar, M. (in prep.). Hunting the Hunters: A Predictive Model of MIS 5 Human Occupation across the Southern Levant, with Application to Understudied Regions in Jordan. *Quaternary International*.

The model uses a **Multi-Criteria Decision Analysis (MCDA)** framework to identify landscape conditions associated with MIS 5 (~130–71 ka) occupation, integrating six topographic and hydrological variables derived from SRTM DEMs. Variable weights are calculated using Kullback–Leibler divergence between archaeological site and background landscape distributions. Sites from the Southern Levant, Northern Levant, and northwestern Arabia (n = 61) are used for model training; the predictive surface and all validation tests apply exclusively to the Southern Levant. A binary geological mask representing access to knappable raw material is applied as an independent exclusionary layer.

The model produces a continuous suitability surface — it does not predict presence or absence of sites, but identifies where the landscape conditions associated with known occupation are most concentrated. For larger study areas or where palaeoenvironmental reconstructions are available, the framework can be extended with past climatic variables to capture higher-resolution patterns of past human behaviour across landscapes.


---

## Repository Structure

```
MIS5-predictive-model/
├── README.md
├── .gitignore
├── MIS5_prediction_model_replicate.R   # Replicate published results (frozen settings)
├── MIS5_prediction_model.R             # Configurable version for new data or regions
└── Inputs/                             # NOT included — see Data Availability
    ├── DEM/
    │   ├── SouthernLevant.tif
    │   ├── Saudi.tif
    │   ├── Lebanon.tif
    │   └── Syria.tif
    ├── Sites/
    │   └── MIS5.csv
    └── Geo_data/
        └── raw_material_presence.tif   # Geology mask (optional)
```

Output folders are created automatically at runtime.

---

## Quick Start

### Replicate the published results

1. Clone the repository.
2. Place the input data in `Inputs/` (see Data Availability).
3. Open `MIS5_prediction_model_replicate.R` in RStudio.
4. Edit the `BASE_DIR` path at the top to match your machine:

```r
BASE_DIR <- "C:/path/to/your/project"
```

5. Run the entire script (Ctrl+Alt+R, or Source).

### Run with new data or a different study area

1. Open `MIS5_prediction_model.R` in RStudio.
2. Edit Section 0 at the top:

```r
# ── EDIT THESE to match your data ─────────────────────────────
BASE_DIR     <- "C:/path/to/your/project"
MODEL_REGION <- "YourRegion"          # must match a DEM filename
```

3. Place your DEM `.tif` files in `Inputs/DEM/`. The region name is taken from the filename:

```
Inputs/DEM/
├── YourRegion.tif      →  region "YourRegion"
├── OtherArea.tif       →  region "OtherArea"
```

4. Place your sites CSV in `Inputs/Sites/` and update `sites_csv` in Section 0 if the filename differs from `MIS5.csv`.
5. Run the entire script.

**Expected runtime:** 30–90 minutes depending on DEM size. If pre-computed Distance to Water rasters are available in `DistToWater/`, Section 3 loads them directly and skips the flow accumulation step.

---

## Methods Summary

### Variables

| Variable | Derivation | Direction |
|---|---|---|
| Elevation | SRTM DEM | Lower preferred |
| Slope | `terra::terrain()` | Higher preferred |
| Aspect | cos(aspect), rescaled [0–1] | North-facing preferred |
| TRI | `terra::terrain()` | Higher preferred |
| TWI | WhiteboxTools D8 flow accumulation | Higher preferred |
| Distance to drainage (DTD) | Euclidean distance from P95 flow accumulation channels | Closer preferred |


### Classification

The continuous suitability surface (0–100) is aggregated to **900 m** resolution and classified into **four quartile classes** (Low / Medium / High / Very High). A binary geological mask excludes areas lacking knappable raw material.

### Validation

Seven validation components are applied to the Southern Levant model region (n = 36 sites):

| Test | What it measures |
|---|---|
| Chi-square enrichment | Whether sites are non-randomly distributed across suitability classes |
| Mann–Whitney U | Whether site suitability scores exceed the landscape background |
| Cohen's d | Effect size of the site–landscape difference |
| Leave-one-out cross-validation (LOO-CV) | Predictive accuracy under weight recomputation |
| Location-based permutation test (n = 999) | Whether observed accuracy exceeds chance placement |
| Clopper–Pearson exact 95% CIs | Precision of classification accuracies given the small sample |
| LOO-CV weight stability | Whether any single site exerts leverage on the weighting scheme |

---

## Requirements

Install once before first use:

```r
install.packages(c("terra", "whitebox", "corrplot", "svglite"))
whitebox::install_whitebox()
```

| Package | Version used | Role |
|---|---|---|
| `terra` | 1.7-78 | Raster analysis |
| `whitebox` | 2.3.4 | Hydrological derivatives (TWI, flow accumulation) |
| `corrplot` | 0.92 | Correlation matrix visualisation |
| `svglite` | 2.1.3 | SVG figure output with editable text |

**R ≥ 4.2.0** required. WhiteboxTools binary is installed automatically.

---

## Data Availability

Input DEMs are 30 m SRTM tiles sourced from [NASA Earthdata](https://earthdata.nasa.gov/) and are not redistributed here. Archaeological site coordinates are available from the corresponding author upon reasonable request, subject to site confidentiality agreements.

The geology mask (`raw_material_presence.tif`) is derived from the [USGS Global Surficial Geology dataset](https://pubs.usgs.gov/of/1997/ofr-97-470/) and the [BGS World Geology map](https://www.bgs.ac.uk/datasets/world-geology/). It is optional; the pipeline runs without it.

---

## Reproducibility

Key reproducibility notes:

- Random seeds are fixed (`set.seed(123)`) for background sampling and permutation tests.
- LOO-CV scoring applies the same [0, 1] clamping as the full model normalisation.
- WhiteboxTools version: 2.3.4 (`whitebox::install_whitebox(version = "2.3.4")`).
- Package versions are listed in the Requirements table above.

---

## Citation

If you use this code or adapt it for your own study area, please cite:

> Samawi, O., Kouki, S., Beller, J. A., Hallinan, E., Rose, J. I., Bicho, N., Nassr, A., Collard, M., & Al-Nahar, M. (in prep.). Hunting the Hunters: A Predictive Model of MIS 5 Human Occupation across the Southern Levant, with Application to Understudied Regions in Jordan. *Quaternary International* (in press).

---

## License

Code: [CC0 1.0 Universal (Public Domain)](https://creativecommons.org/publicdomain/zero/1.0/) — no rights reserved. You may copy, modify, and distribute this code for any purpose without permission or attribution, though citation is appreciated.
Site data: available on request (see Data Availability).
