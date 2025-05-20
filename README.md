# ithir R Package

The `ithir` package provides data, models, and utility functions for digital soil informatics. It supports soil data analysis from both raster (gridded) and profile (point) formats, with tools for interpolation, prediction assessment, and general soil model diagnostics.

---

## 🔧 Key Features

- **Mass-preserving spline interpolation** for soil profile data (point and raster)
- **Fast raster spline implementation** using optimized matrix structures and `terra::app()`
- **Model evaluation tools**:
  - For **continuous predictions**: RMSE, concordance correlation, R², bias
  - For **categorical predictions**: confusion matrix, kappa, accuracy, entropy
- **Built-in example datasets** for quick testing and demonstration

---

## 📦 Installation

To install the development version from GitHub:

```r
install.packages("devtools")  # if needed
devtools::install_github("brendo1001/ithir_github/pkg")
