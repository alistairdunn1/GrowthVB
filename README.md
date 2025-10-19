# growthVB: Von Bertalanffy Growth Curve Estimation for R

`growthVB` is an R package for estimating von Bertalanffy growth curves from age and length data. It provides both maximum likelihood estimation (MLE) and Bayesian (brms) approaches to fitting, visualization tools for diagnostics, and methods for summarizing parameter estimates.

## Installation

```r
# Install directly from GitHub (recommended)
devtools::install_github("alistairdunn1/growthVB", subdir = "growthVB")

# Or install from a local directory (from the root directory containing growthVB folder)
# install.packages("growthVB", repos = NULL, type = "source")

# Or if you're in a different directory, specify the full path
# install.packages("path/to/GrowthVB/growthVB", repos = NULL, type = "source")

# To build the package from source, run:
# Rscript build_package.R
```

## Features

- Fit von Bertalanffy growth curves using maximum likelihood estimation (MLE)
- Optionally fit Bayesian von Bertalanffy curves using brms
- Explicitly model heteroscedasticity with CV as a function of predicted length
- Sex-specific growth modeling with separate parameters by sex
- Length bin sampling corrections
- Parameter estimation with 95% confidence/credible intervals
- Predictions with confidence and prediction intervals
- Visualize growth curves and underlying data
- Produce diagnostic plots including empirical CV analysis
- Generate age-length heatmaps with optional smoothers and age frequency summaries
- Posterior predictive checks for Bayesian models

## Core Functions

### Model Fitting

#### `fit_vb_mle()`: Frequentist Von Bertalanffy Model Fitting

```r
fit_vb_mle(
  age,                 # Numeric vector of ages
  length,              # Numeric vector of lengths
  sex = NULL,          # Optional factor for sex-specific models (e.g., "M", "F")
  length_bins = NULL,  # Optional numeric vector of length bin midpoints
  sampling_prob = 1,   # Optional vector of sampling probabilities (defaults to 1)
  ci_level = 0.95,     # Confidence interval level (default 0.95)
  optim_method = "L-BFGS-B", # Optimization method for stats::optim (default "L-BFGS-B")
  maxit = 1000         # Maximum number of iterations for optimization
)
```

**Features:**

- **Maximum likelihood estimation**: Uses direct optimization for flexible parameter estimation
- **Sex-specific modeling**: When `sex` is provided, fits separate growth curves for each sex
- **CV modeling**: Explicitly models coefficient of variation as a function of mean length
- **Heteroscedastic variance**: Accounts for increasing variance with fish length
- **Parameter uncertainty**: Provides standard errors and confidence intervals for all parameters
- **Robust fitting**: Uses bounded optimization (L-BFGS-B) with intelligent starting values for reliable parameter estimation and confidence intervals

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

### Visualization

- `plot_vb()`: Plot growth curves with data points
- `plot_vb_diagnostics()`: Produce standard regression diagnostics plots
- `plot_vb_posteriors()`: Generate posterior predictive plots for Bayesian models
- `plot_age_length_heatmap()`: Create heatmap visualizations of age-length observations with optional smoothers
- `plot_vb_age_counts()`: Plot frequency distribution of age samples by group variables
- `plot_empirical_cv()`: Plot empirical coefficient of variation by age for data exploration

## Example Usage

### Basic Usage

```r
library(growthVB)

# Simple example with simulated data (heteroscedastic errors)
set.seed(123)
true_Linf <- 120
true_k <- 0.25
true_t0 <- -0.5
true_cv <- 0.1

# Generate ages (multiple observations per age)
age <- rep(1:15, each = 6)
age <- age + rnorm(length(age), 0, 0.1)  # Add some noise

# Generate lengths with heteroscedastic error (SD = CV * predicted length)
mean_length <- true_Linf * (1 - exp(-true_k * (age - true_t0)))
sd_length <- true_cv * mean_length
length <- rnorm(length(age), mean = mean_length, sd = sd_length)

# Fit model with MLE and heteroscedastic errors
fit <- fit_vb_mle(age = age, length = length)

# Print model results
print(fit)

# Summarize parameters
summarize_vb(fit)

# Plot growth curve
plot_vb(fit)

# Generate diagnostic plots
diagnostics <- plot_vb_diagnostics(fit)

# Create age-length heatmap with smoother
heatmap_plot <- plot_age_length_heatmap(age = age, length = length, add_smoother = TRUE)

# Plot empirical CV by age
cv_plot <- plot_empirical_cv(age = age, length = length)

# Generate predictions with confidence intervals
new_ages <- data.frame(age = seq(1, 15, by = 1))
predictions <- predict(fit, newdata = new_ages, interval = "confidence")
print(predictions)
```

### Sex-specific Growth Models

```r
# Create data with sex differences and heteroscedastic errors
set.seed(456)
n_per_sex <- 50
ages <- seq(1, 15, length.out = n_per_sex)

# Create sex variable
sex <- rep(c("F", "M"), each = n_per_sex)

# Different parameters by sex
# Females: larger maximum size (Linf), slower growth rate (k)
Linf_f <- 120
k_f <- 0.15
t0_f <- -0.5

# Males: smaller maximum size, faster growth rate
Linf_m <- 90
k_m <- 0.25
t0_m <- -0.3

# Common CV parameter
cv <- 0.08

# Calculate mean lengths
mean_length_f <- Linf_f * (1 - exp(-k_f * (ages - t0_f)))
mean_length_m <- Linf_m * (1 - exp(-k_m * (ages - t0_m)))
mean_length <- c(mean_length_f, mean_length_m)

# Calculate standard deviations (heteroscedastic)
sd_length <- cv * mean_length

# Generate observed lengths
age <- rep(ages, 2)
length <- rnorm(length(age), mean = mean_length, sd = sd_length)

# Fit sex-specific model
fit_sex <- fit_vb_mle(age = age, length = length, sex = sex)

# Print results
print(fit_sex)

# Summarize parameters
summarize_vb(fit_sex)

# Plot sex-specific growth curves
plot_vb(fit_sex)

# Plot age counts by sex
age_counts_plot <- plot_vb_age_counts(age = age, sex = sex)

# Compare CV patterns between sexes
cv_plot_sex <- plot_empirical_cv(age = age, length = length, sex = sex)

# Make predictions for specific ages by sex
new_data <- data.frame(
  age = rep(c(1, 5, 10, 15), 2),
  sex = rep(c("F", "M"), each = 4)
)
predictions <- predict(fit_sex, newdata = new_data, interval = "confidence")
print(predictions)
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
fit_binned <- fit_vb_mle(
  age = age, 
  length = length,
  length_bins = as.numeric(names(bin_counts)),
  bin_counts = as.numeric(bin_counts)
)

# Compare with uncorrected model
fit_uncorrected <- fit_vb_mle(age = age, length = length)

# Compare parameter estimates
cat("Binned model parameters:\n")
print(fit_binned$parameters)
cat("Uncorrected model parameters:\n") 
print(fit_uncorrected$parameters)
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
