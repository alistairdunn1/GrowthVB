# growthVB: Von Bertalanffy Growth Curve Estimation for R

`growthVB` is an R package for estimating von Bertalanffy growth curves from age and length data. It provides both frequentist (nls) and Bayesian (brms) approaches to fitting, visualization tools for diagnostics, and methods for summarizing parameter estimates.

## Installation

```r
# Install from a local directory
install.packages("path/to/growthVB", repos = NULL, type = "source")

# Or install directly from GitHub (if hosted there)
# devtools::install_github("username/growthVB")
```

## Features

- Fit von Bertalanffy growth curves using non-linear least squares (NLS)
- Optionally fit Bayesian von Bertalanffy curves using brms
- Model CV as a function of mean length
- Sex-specific growth modeling with optional interaction terms
- Length bin sampling corrections
- Parameter estimation with 95% confidence/credible intervals
- Visualize growth curves and underlying data
- Produce diagnostic plots
- Generate age-length heatmaps and age frequency summaries
- Posterior predictive checks for Bayesian models

## Core Functions

### Model Fitting

#### `fit_vb_nls()`: Frequentist Von Bertalanffy Model Fitting

```r
fit_vb_nls(
  age,                 # Numeric vector of ages
  length,              # Numeric vector of lengths
  sex = NULL,          # Optional factor for sex-specific models (e.g., "M", "F")
  length_bins = NULL,  # Optional numeric vector of length bin midpoints
  bin_counts = NULL,   # Optional count of samples in each length bin
  start = NULL,        # Optional starting values for parameters (list: Linf, k, t0)
  cv_model = TRUE,     # Whether to model CV as function of mean length
  control = NULL,      # Control parameters for nls fitting (list)
  verbose = TRUE,      # Whether to print progress and warnings
  ...                  # Additional arguments passed to nls
)
```

**Features:**
- **Sex-specific modeling**: When `sex` is provided, fits separate growth curves for each sex
- **Length bin sampling**: Corrects for length-stratified sampling through `length_bins` and `bin_counts`
- **CV modeling**: Models heteroscedasticity where variance increases with fish size
- **Robust fitting**: Includes intelligent starting value estimation and convergence handling

#### `fit_vb_brms()`: Bayesian Von Bertalanffy Model Fitting

```r
fit_vb_brms(
  age,                  # Numeric vector of ages
  length,               # Numeric vector of lengths
  sex = NULL,           # Optional factor for sex-specific models (e.g., "M", "F")
  length_bins = NULL,   # Optional numeric vector of length bin midpoints
  bin_counts = NULL,    # Optional count of samples in each length bin
  priors = NULL,        # Optional list of prior specifications
  iter = 2000,          # MCMC iterations
  warmup = 1000,        # Warmup iterations
  chains = 4,           # Number of MCMC chains
  cores = getOption("mc.cores", 1), # Number of cores for parallel processing
  cv_model = TRUE,      # Whether to model CV as function of mean length
  verbose = TRUE,       # Whether to print progress and diagnostics
  ...                   # Additional arguments passed to brms
)
```

**Features:**
- **Full Bayesian inference**: Posterior distributions for all parameters
- **Sex-specific modeling**: When `sex` is provided, fits separate growth parameters for each sex
- **Length bin sampling**: Corrects for length-stratified sampling through weighting
- **CV modeling**: Models heteroscedasticity where variance increases with fish size
- **Prior specification**: Flexible prior definition for all parameters
- **MCMC control**: Fine-grained control of sampling behavior

### Analysis
- `summarize_vb()`: Summarize parameter estimates with confidence/credible intervals
- `compare_vb()`: Statistical comparison of multiple growth models

### Visualization
- `plot_vb()`: Plot growth curves with data points
- `plot_vb_diagnostics()`: Produce standard regression diagnostics plots
- `plot_vb_posteriors()`: Generate posterior predictive plots for Bayesian models
- `plot_age_length_heatmap()`: Create heatmap visualizations of age-length observations
- `plot_vb_age_counts()`: Plot frequency distribution of age samples

## Example Usage

### Basic Usage

```r
library(growthVB)

# Simple example with simulated data
age <- 1:15
length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5)

# Fit model
fit <- fit_vb_nls(age = age, length = length)

# Summarize parameters
summarize_vb(fit)

# Plot growth curve
plot_vb(fit)
```

### Sex-specific Growth Models

```r
# Create data with sex differences
age <- rep(1:15, 2)
sex <- rep(c("F", "M"), each = 15)
length_f <- 120 * (1 - exp(-0.15 * (1:15 - (-0.5)))) + rnorm(15, 0, 5)
length_m <- 90 * (1 - exp(-0.25 * (1:15 - (-0.5)))) + rnorm(15, 0, 5)
length <- c(length_f, length_m)

# Fit sex-specific model
fit_sex <- fit_vb_nls(age = age, length = length, sex = sex)

# Plot sex-specific growth curves
plot_vb(fit_sex)
```

### Length Bin Sampling Correction

```r
# Example with length bin sampling
age <- sample(1:15, 100, replace = TRUE)
true_length <- 100 * (1 - exp(-0.2 * (age - (-0.5))))
length <- true_length + rnorm(100, 0, 5)

# Define length bins and counts for sampling correction
length_bins <- seq(10, 110, by = 10)
bin_counts <- table(cut(length, breaks = seq(5, 115, by = 10), labels = length_bins))

# Fit with length bin correction
fit_binned <- fit_vb_nls(
  age = age, 
  length = length,
  length_bins = as.numeric(names(bin_counts)),
  bin_counts = as.numeric(bin_counts)
)

# Compare with uncorrected model
fit_uncorrected <- fit_vb_nls(age = age, length = length)
compare_vb(list(Corrected = fit_binned, Uncorrected = fit_uncorrected))
```

### Bayesian Model

```r
library(growthVB)

# Only runs if brms is available
if (requireNamespace("brms", quietly = TRUE)) {
  # Fit Bayesian model
  fit_bayes <- fit_vb_brms(
    age = age, 
    length = length,
    priors = list(
      Linf = c(100, 20),  # prior mean and SD for Linf
      k = c(0.2, 0.1),    # prior mean and SD for k
      t0 = c(0, 1)        # prior mean and SD for t0
    ),
    chains = 2,           # using fewer chains for example speed
    iter = 1000
  )
  
  # Plot posterior predictions
  plot_vb_posteriors(fit_bayes)
}
```

For a complete example, see the example script in the package:
```r
file.path(system.file(package = "growthVB"), "examples", "example_usage.R")
```

## References

Von Bertalanffy, L. (1957). Quantitative laws in metabolism and growth. *The Quarterly Review of Biology*, 32(3), 217-231.

Kimura, D. K. (2008). Extending the von Bertalanffy growth model using explanatory variables. *Canadian Journal of Fisheries and Aquatic Sciences*, 65(9), 1879-1891.

Okuda, T., Somhlaba, S., Sarralde, R., Mori, M., Rojo, V., & Dunn, A. (2025). Characterisation of the toothfish fishery in Subarea 48.6 through the 2024/25 season. CCAMLR Working Group Paper, WG-FSA-2025/34, Hobart.

## License

MIT
