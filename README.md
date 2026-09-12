# MIS 5 Predictive Model of Human Occupation — Southern Levant

[![License: CC0 1.0](https://img.shields.io/badge/License-CC0%201.0-lightgrey.svg)](https://creativecommons.org/publicdomain/zero/1.0/)
[![R ≥ 4.2](https://img.shields.io/badge/R-%E2%89%A5%204.2-276DC3?logo=r)](https://www.r-project.org/)

---

## Overview

This repository contains the complete R pipeline for a **predictive model of MIS 5 human occupation** across the Southern Levant, as described in:

> Samawi, O., Kouki, S., Beller, J. A., Hallinan, E., Rose, J. I., Bicho, N., Nassr, A., Collard, M., & Al-Nahar, M. (in prep.). Hunting the Hunters: A Predictive Model of MIS 5 Human Occupation across the Southern Levant, with Application to Understudied Regions in Jordan. *Quaternary International*.

The model uses a **Multi-Criteria Decision Analysis (MCDA)** framework to identify landscape conditions associated with MIS 5 (~130–71 ka) occupation, integrating six topographic and hydrological variables derived from SRTM DEMs. Variable weights are calculated using Kullback–Leibler divergence between archaeological site and background landscape distributions. Sites from the Southern Levant, Northern Levant, and northwestern Arabia (n = 61) are used for model training; the predictive surface and all validation tests apply exclusively to the Southern Levant. A binary geological mask representing access to knappable raw material is applied as an independent exclusionary layer.

The model produces a continuous suitability surface — it does not predict presence or absence of sites, but identifies where the landscape conditions associated with known occupation are most concentrated. For larger study areas or where palaeoenvironmental reconstructions are available, the framework can be extended with past climatic variables to capture higher-resolution patterns of past human behaviour across landscapes.

Code development, pipeline implementation, and statistical validation by **S. Kouki**.

---

## Repository Structure

```
MIS5-predictive-model/
├── run_pipeline.R                        # Command-line entry point
├── MIS5_prediction_model_MERGED.R        # Full pipeline (all sections)
├── MIS5_prediction_model_replicate.R     # Original code to replicate published results
├── DESCRIPTION                           # R package metadata and dependencies
├── README.md
└── Inputs/                               # NOT included — see Data Availability
    ├── DEM/                              # Drop any number of regional .tif DEMs here
    │   ├── SouthernLevant.tif            # Region name = filename stem
    │   ├── Saudi.tif
    │   ├── Lebanon.tif
    │   └── Syria.tif
    ├── Sites/
    │   └── MIS5.csv                      # Site coordinates and attributes
    └── Geo_data/
        └── raw_material_presence.tif     # Geology mask (optional)
```

Output folders are created automatically at runtime and are gitignored:

```
<base_dir>/
└── outputs/                   # or whatever --output name you provide
    ├── Topo/                  # Slope, Aspect, TRI, TWI, DistToWater (per region)
    ├── Sites_data/            # Extracted site variables, background sample
    ├── Stats/                 # Descriptive stats, correlations, weights, validation metrics,
    │                          #   LOO-CV weight stability
    ├── MCDA/                  # Suitability rasters (30 m, 900 m, 4-class, masked)
    └── Figures/               # All plots (SVG and PNG)
```

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

### Weighting

Variable weights are derived from the **Kullback–Leibler divergence** between site and background distributions (Kullback & Leibler, 1951), following the information-theoretic approach to archaeological variable weighting. Weights are calculated both globally (full modelled area) and per bioclimatic zone (Mediterranean, Irano-Turanian, Saharo-Arabian).

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

Run once in R or RStudio before first use:

```r
install.packages(c("terra", "whitebox", "corrplot", "svglite"))
whitebox::install_whitebox()
```

| Package | Version tested | Role |
|---|---|---|
| `terra` | ≥ 1.7-0 | Raster analysis |
| `whitebox` | ≥ 2.3.0 | Hydrological derivatives (TWI, flow accumulation) |
| `corrplot` | ≥ 0.92 | Correlation matrix visualisation |
| `svglite` | ≥ 2.1.0 | SVG figure output with editable text |

**R ≥ 4.2.0** is required. WhiteboxTools is installed automatically via `whitebox::install_whitebox()`.

> **Windows note:** if `Rscript` is not recognised in your terminal, add R to your PATH or use the full path:
> `"C:/Program Files/R/R-4.x.x/bin/Rscript.exe" run_pipeline.R ...`

---

## Usage

There are two ways to use this repository depending on your goal.

---

### Option A — Replicate the published results

Use `MIS5_prediction_model_replicate.R`. This is the original, unmodified analysis script that produced the results reported in the paper. It is self-contained and requires no command-line setup.

1. Open `MIS5_prediction_model_replicate.R` in RStudio.
2. Edit the `BASE_DIR` path at the top of the script to point to your local copy of the input data:

```r
BASE_DIR <- "C:/path/to/your/data"
```

3. Run the script top-to-bottom (Ctrl+Alt+R in RStudio, or Source).

This script uses the exact same input files, thresholds, and classification as the published analysis. It does not expose command-line arguments or configurable regions — it is a fixed, reproducible record of the original run. If you need to change any of these configurations, use the `MIS5_prediction_model_MERGED.R` script and change the directory, folder, region names and configuration values.

---

### Option B — Run the full pipeline (new data or new regions)

Use `run_pipeline.R` + `MIS5_prediction_model_MERGED.R`. This is the generalised version of the pipeline: it discovers regions automatically from the DEM folder, accepts command-line arguments, and can be applied to any study area.

#### Step 1 — Prepare your input folder

Place your regional DEM `.tif` files in `Inputs/DEM/`. The region name is taken directly from the filename — no configuration file needed:

Our regions:
```
Inputs/DEM/
├── SouthernLevant.tif    →  region "SouthernLevant"
├── Saudi.tif             →  region "Saudi"
├── Lebanon.tif           →  region "Lebanon"
└── Syria.tif             →  region "Syria"
```

Any number of DEMs is supported. For a different study area, simply replace the files:

```
Inputs/DEM/
├── Region1.tif           →  region "Region1"
├── Region2.tif           →  region "Region2"
└── Region3.tif           →  region "Region3"
```

#### Step 2 — Open a terminal and run

**Minimal call** (uses default output folder name):

```bash
Rscript run_pipeline.R \
    --base_dir "your/project/directory" \
    --model_region "SouthernLevant"
```

**With a custom output folder:**

```bash
Rscript run_pipeline.R \
    --base_dir "your/project/directory" \
    --model_region "SouthernLevant" \
    --output "my_outputs"
```

**Different study area, different sites file:**

```bash
Rscript run_pipeline.R \
    --base_dir "/home/user/my_project" \
    --model_region "MyRegion" \
    --sites_file "my_sites.csv" \
    --output "outputs"
```

**Skip the geology mask** (even if the file exists):

```bash
Rscript run_pipeline.R \
    --base_dir "/home/user/my_project" \
    --model_region "MyRegion" \
    --no_geo_mask
```

#### All available arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--base_dir` | ✓ | — | Root project directory |
| `--model_region` | ✓ | — | DEM filename stem for the model region |
| `--sites_file` | | `MIS5.csv` | CSV filename inside `Inputs/Sites/` |
| `--output` | | `outputs` | Output folder name under `base_dir` |
| `--no_geo_mask` | | (mask used if found) | Skip geology mask |
| `--help` | | — | Print usage and exit |

#### What `--model_region` means

The pipeline computes topographic variables for **all** DEMs in `Inputs/DEM/`. The `--model_region` argument tells it which single region to use for the suitability model, validation, and figures. It must exactly match a DEM filename stem. If it doesn't, the script will print the available names and exit cleanly.

#### Running interactively in RStudio (no terminal)

You can also run `MIS5_prediction_model_MERGED.R` directly in RStudio. Edit the defaults at the top of Section 0:

```r
BASE_DIR     <- "C:/your/path"
MODEL_REGION <- "SouthernLevant"   # must match a DEM filename stem
```

Then Source the script. All CLI arguments are optional when running interactively — the script detects whether they have been set and falls back to the inline defaults if not.

**Expected runtime:** 30–90 minutes depending on DEM extent and available RAM. The WhiteboxTools flow accumulation step (Sections 2–3) is the most compute-intensive. If pre-computed Distance to Water rasters are available in `DistToWater/`, Section 3 loads them directly and skips the flow accumulation computation.

---

## Data Availability

Input DEMs are 30 m SRTM tiles sourced from [NASA Earthdata](https://earthdata.nasa.gov/) and are not redistributed here. Archaeological site coordinates are available from the corresponding author upon reasonable request, subject to site confidentiality agreements.

The geology mask (`raw_material_presence.tif`) is derived from the [USGS Global Surficial Geology dataset](https://pubs.usgs.gov/of/1997/ofr-97-470/) and the [BGS World Geology map](https://www.bgs.ac.uk/datasets/world-geology/). It is optional; the pipeline runs without it.

---

## Reproducibility

The original published results can be replicated exactly using `MIS5_prediction_model_replicate.R` with the archived input data (see Data Availability). Key reproducibility notes:

- Random seeds are fixed (`set.seed(123)`) for background sampling and permutation tests.
- LOO-CV scoring applies the same [0, 1] clamping as the full model normalisation.
- WhiteboxTools version should match: `whitebox::install_whitebox(version = "2.3.4")`.
- Full R environment is recorded in `renv.lock`. Restore with `renv::restore()`.

---

## Citation

If you use this code or adapt it for your own study area, please cite:

> Samawi, O., Kouki, S., Beller, J. A., Hallinan, E., Rose, J. I., Bicho, N., Nassr, A., Collard, M., & Al-Nahar, M. (in prep.). Hunting the Hunters: A Predictive Model of MIS 5 Human Occupation across the Southern Levant, with Application to Understudied Regions in Jordan. *Quaternary International*.

---

## License

Code: [CC0 1.0 Universal (Public Domain)](https://creativecommons.org/publicdomain/zero/1.0/) — no rights reserved. You may copy, modify, and distribute this code for any purpose without permission or attribution, though citation is appreciated.
Site data: available on request (see Data Availability).
