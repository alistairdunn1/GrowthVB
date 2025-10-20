# growthVB: Von Bertalanffy Growth Curve Estimation for R

`growthVB` is an R package for estimating von Bertalanffy growth curves from age and length data. It provides maximum likelihood estimation (MLE), Bayesian (brms), and spatial modelling approaches, with visualisation tools for diagnostics and methods for summarising parameter estimates.

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

- **Standard von Bertalanffy models**: Fit growth curves using maximum likelihood estimation (MLE) or Bayesian methods (brms)
- **Heteroscedasticity modelling**: Explicitly model CV as a function of predicted length
- **Sex-specific growth modelling**: Separate parameters by sex with proper statistical treatment
- **Length bin sampling corrections**: Account for stratified sampling designs
- **Parameter estimation**: 95% confidence/credible intervals for all parameters
- **Diagnostics**: Residual analysis, posterior predictive checks, empirical CV analysis

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
  optim_method = "L-BFGS-B", # Optimisation method for stats::optim (default "L-BFGS-B")
  maxit = 1000         # Maximum number of iterations for optimisation
)
```

**Features:**

- **Maximum likelihood estimation**: Uses direct optimisation for flexible parameter estimation
- **Sex-specific modelling**: When `sex` is provided, fits separate growth curves for each sex
- **CV modelling**: Explicitly models coefficient of variation as a function of mean length
- **Heteroscedastic variance**: Accounts for increasing variance with fish length
- **Parameter uncertainty**: Provides standard errors and confidence intervals for all parameters
- **Robust fitting**: Uses bounded optimisation (L-BFGS-B) with intelligent starting values for reliable parameter estimation and confidence intervals

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
- **Sex-specific modelling**: When `sex` is provided, fits separate growth parameters for each sex
- **Length bin sampling**: Corrects for length-stratified sampling through weighting
- **CV modelling**: Models heteroscedasticity where variance increases with fish size
- **Prior specification**: Flexible prior definition for all parameters
- **MCMC control**: Fine-grained control of sampling behaviour

#### Analysis

- `summarise_vb()`: Summarise parameter estimates with confidence/credible intervals

### Visualisation

#### Standard Plotting

- `plot_vb()`: Plot growth curves with data points
- `plot_vb_mle_diagnostics()`: Produce standard regression diagnostics plots for MLE models
- `plot_vb_posteriors()`: Generate posterior predictive plots for Bayesian models
- `plot_vb_bayes_diagnostics()`: Comprehensive MCMC convergence and Bayesian model diagnostics
- `plot_vb_growth_pp_checks()`: Growth-specific posterior predictive checks for Bayesian models
- `plot_vb_predictions()`: Plot von Bertalanffy growth curves from prediction data with confidence/prediction intervals
- `plot_age_length_heatmap()`: Create heatmap visualisations of age-length observations with optional smoothers
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

# Summarise parameters
summarise_vb(fit)

# Plot growth curve
plot_vb(fit)

# Generate diagnostic plots
diagnostics <- plot_vb_mle_diagnostics(fit)

# Create age-length heatmap with smoother
heatmap_plot <- plot_age_length_heatmap(age = age, length = length, add_smoother = TRUE)

# Plot empirical CV by age
cv_plot <- plot_empirical_cv(age = age, length = length)

# Generate predictions with confidence intervals
new_ages <- data.frame(age = seq(1, 15, by = 1))
predictions <- predict(fit, newdata = new_ages, interval = "confidence")
print(predictions)

# Plot predictions with original data
plot_vb_predictions(predictions, original_data = fit$data)
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

# Summarise parameters
summarise_vb(fit_sex)

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

# Create detailed predictions over age range for plotting
new_ages <- data.frame(
  age = rep(seq(0, 20, by = 0.5), each = 2), 
  sex = rep(c("F", "M"), times = 41)
)
predictions_detailed <- predict(fit_sex, newdata = new_ages, interval = "prediction")

# Plot growth curves with prediction intervals and original data
plot_vb_predictions(predictions_detailed, original_data = fit_sex$data)
```

### Advanced Prediction Plotting

The `plot_vb_predictions()` function provides comprehensive visualisation of von Bertalanffy growth predictions:

```r
# Fit a model
fit <- fit_vb_mle(age = age, length = length, sex = sex)

# Generate predictions over desired age range
new_ages <- data.frame(
  age = rep(seq(0, 25, by = 1), each = 2), 
  sex = rep(c("M", "F"), times = 26)
)
predictions <- predict(fit, newdata = new_ages, interval = "prediction")

# Create comprehensive plot with all features
plot_vb_predictions(
  predictions = predictions,
  original_data = fit$data,        # Show original data points
  show_points = TRUE,              # Display data points
  show_intervals = TRUE,           # Show prediction intervals
  facet_by_sex = TRUE,            # Separate panels for each sex
  alpha_ribbon = 0.3,             # Transparency for intervals
  alpha_points = 0.5,             # Transparency for points
  point_size = 0.9,               # Size of data points
  line_size = 1,                  # Thickness of fitted curves
  title = "Von Bertalanffy Growth Curves with Prediction Intervals"
)

# Single model without sex
fit_single <- fit_vb_mle(age = age, length = length)
new_ages_single <- data.frame(age = seq(0, 25, by = 1))
predictions_single <- predict(fit_single, newdata = new_ages_single, interval = "confidence")

# Plot without faceting
plot_vb_predictions(predictions_single, original_data = fit_single$data, facet_by_sex = FALSE)
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
      t0 = c(0, 1),       # prior mean and SD for t0
      CV = c(0.1, 0.05)   # prior mean and SD for CV
    ),
    chains = 2,           # using fewer chains for example speed
    iter = 1000
  )
  
  # Basic posterior plots
  plot_vb_posteriors(fit_bayes)
  
  # Comprehensive Bayesian diagnostics
  full_diagnostics <- plot_vb_bayes_diagnostics(fit_bayes)
  
  # Access specific diagnostic categories
  full_diagnostics$convergence$trace    # MCMC trace plots
  full_diagnostics$convergence$rhat     # R-hat convergence diagnostic
  full_diagnostics$pp_checks$density    # Posterior predictive density check
  full_diagnostics$residuals$residual_scatter  # Residual diagnostics
  
  # Growth-specific posterior predictive checks
  growth_checks <- plot_vb_growth_pp_checks(fit_bayes)
  growth_checks$growth_patterns         # Growth trajectory validation
  growth_checks$age_effects            # Age-specific CV patterns
  growth_checks$length_distribution$qq_plot  # Q-Q plot comparison
}
```

### Citation

To cite growthVB in publications, please use:

```r
citation("growthVB")
```

This will provide the appropriate citation format. For manual citation:

> Dunn, A. (2025). growthVB: Von Bertalanffy Growth Curve Estimation with Heteroscedastic Errors. R package version 0.2.0. [https://github.com/alistairdunn1/growthVB](https://github.com/alistairdunn1/growthVB)

When using specific functionality, please also cite:

- **Bayesian models**: Bürkner (2017) for brms, Carpenter et al. (2017) for Stan
- **von Bertalanffy equation**: von Bertalanffy (1957)
- **Heteroscedastic growth models**: Kimura (2008)

## References

Bürkner, P. C. (2017). brms: An R package for Bayesian multilevel models using Stan. *Journal of Statistical Software*, 80(1), 1-28.

Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt, M., ... & Riddell, A. (2017). Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1).

Kimura, D. K. (2008). Extending the von Bertalanffy growth model using explanatory variables. *Canadian Journal of Fisheries and Aquatic Sciences*, 65(9), 1879-1891.

Okuda, T., Somhlaba, S., Sarralde, R., Mori, M., Rojo, V., & Dunn, A. (2025). Characterisation of the toothfish fishery in Subarea 48.6 through the 2024/25 season. CCAMLR Working Group Paper, WG-FSA-2025/34, Hobart.

Von Bertalanffy, L. (1957). Quantitative laws in metabolism and growth. *The Quarterly Review of Biology*, 32(3), 217-231.

## License

MIT
