# turbulence_R_code
# Turbulence Dataset Analysis (Laplace Turbulence Index; LTI)

This repository provides an R script to reproduce the main analyses and figures for the study:

**“Proposal of a Laplace Turbulence Index (LTI) Based on the Non-Gaussianity of Temporal Wind Velocity Increments.”**

The workflow downloads the public dataset from Zenodo, preprocesses the time series, and generates figures that evaluate:
- the non-Gaussian (Laplace-type) behavior of wind-velocity increments,
- robustness to projection angle,
- model comparison (Normal vs Laplace) using AIC / log-likelihood,
- time-lag (Δt) dependence and scaling exponent,
- consistency between pole-mounted and drone-mounted ultrasonic anemometer measurements.

---

## Contents

- **R script**: 

---

## Requirements

### R packages
This script uses the following packages:

- archive
- httr2
- readr
- dplyr
- lubridate
- ggplot2
- data.table
- signal
- purrr
- stringr
- tibble
- tidyr

Install packages:

```r
install.packages(c(
  "archive","httr2","readr","dplyr","lubridate","ggplot2",
  "data.table","signal","purrr","stringr","tibble","tidyr"
))
```

Tested environment

R 4.x recommended

Windows / macOS / Linux (tested mainly on Windows)

Dataset

The script downloads:

Zenodo record: 17988010

ZIP file: turbulence_dataset.zip

Dataset will be extracted under:

~/zenodo_turbulence_dataset/

If you already downloaded the dataset manually, set:

dataset_root <- "PATH/TO/zenodo_turbulence_dataset"

You can see 24H observation data
http://ito-lab.mobi:3000/public-dashboards/23d519078ada48128335348fc2730135
