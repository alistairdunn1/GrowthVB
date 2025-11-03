# growthVB: Von Bertalanffy Growth Curve Estimation for R

`growthVB` is an R package for estimating von Bertalanffy growth curves from age and length data. It provides maximum likelihood estimation (MLE) and Bayesian (brms) approaches, with statistical testing, visualisation tools, and output methods.

**Current Version:** 0.2.0

## Installation

```r
# Install from GitHub
devtools::install_github("alistairdunn1/growthVB", subdir = "growthVB")

# Or install from local directory
install.packages("growthVB", repos = NULL, type = "source")

# Build from source
Rscript build_package.R
```

## Features

- **Standard von Bertalanffy models**: Fit growth curves using MLE or Bayesian methods (brms)
- **Heteroscedasticity modelling**: Model CV as a function of predicted length
- **Sex-specific growth modelling**: Separate parameters by sex
- **Length bin sampling corrections**: Account for stratified sampling designs
- **Parameter estimation**: 95% confidence/credible intervals
- **Diagnostics**: Residual analysis, posterior predictive checks, empirical CV analysis
- **Growth curve comparison**: Test differences in growth trajectories between groups
- **Visualisation**: Density plots showing bootstrap permutation distributions
- **Summary and print methods**: Formatted output for results
- **Documentation**: Help with examples for all functionality

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

- **Maximum likelihood estimation**: Uses optimisation for parameter estimation
- **Sex-specific modelling**: When `sex` is provided, fits separate growth curves for each sex
- **CV modelling**: Models coefficient of variation as a function of mean length
- **Heteroscedastic variance**: Accounts for increasing variance with fish length
- **Parameter uncertainty**: Provides standard errors and confidence intervals
- **Fitting**: Uses bounded optimisation (L-BFGS-B) with starting values

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

- **Bayesian inference**: Posterior distributions for all parameters
- **Sex-specific modelling**: When `sex` is provided, fits separate growth parameters for each sex
- **Length bin sampling**: Corrects for length-stratified sampling through weighting
- **CV modelling**: Models heteroscedasticity where variance increases with fish size
- **Prior specification**: Prior definition for all parameters
- **MCMC control**: Control of sampling behaviour

#### Analysis

- `summarise_vb()`: Summarise parameter estimates with confidence/credible intervals

### Model Output Structure

Both `fit_vb_mle()` and `fit_vb_brms()` return list objects containing:

- **`parameters`**: Parameter estimates with standard errors/credible intervals
- **`data`**: Original data with analysis results, including:
  - `ID`: Unique identifier preserving original input order (1, 2, 3, ...)
  - `age`: Original age values
  - `length`: Original length values
  - `sex`: Sex classification (if provided)
  - `fitted`: Fitted/predicted length values from the von Bertalanffy model
  - `residual`: Residuals calculated as observed - fitted values
- **`model`**: The fitted model object (nls for MLE, brmsfit for Bayesian)
- **`call`**: The function call that created the model

### Statistical Testing

#### `compare_vb_mle()`: Bootstrap Permutation Tests for Group Comparisons

```r
compare_vb_mle(
  age,                    # Numeric vector of ages
  length,                 # Numeric vector of lengths  
  group,                  # Factor or character vector specifying groups
  sex = NULL,             # Optional factor for sex-specific models
  n_bootstrap = 1000,     # Number of bootstrap permutations (default 1000)
  parameters = c("Linf", "k", "t0", "CV"), # Parameters to test (default all)
  alpha = 0.05,           # Significance level (default 0.05)
  verbose = TRUE,         # Print progress messages (default TRUE)
  seed = NULL,            # Optional seed for reproducible results
  min_obs = 50,           # Minimum observations per group/sex combination
  age_stratified = TRUE,  # Use age-stratified permutation (default TRUE)
  age_bin_width = 2,      # Width of age bins in years (default 2)
  test_curves = TRUE,     # Test growth curve differences (default TRUE)
  curve_ages = NULL       # Ages for curve evaluation (default: seq from min to max)
)
```

**Key Features:**

- **Bootstrap permutation testing**: Non-parametric statistical tests for parameter differences
- **Age-stratified permutation**: Maintains age distributions for valid tests
- **Sex-specific model support**: Handles complex sex-specific parameter comparisons
- **Growth curve testing**: Tests overall trajectory differences beyond individual parameters
- **Data filtering**: Automatic removal of groups/combinations with insufficient data
- **Comprehensive output**: P-values, effect sizes, and detailed method information

### Visualisation

#### Standard Plotting

- `plot_vb()`: Plot growth curves with data points
- `plot_vb_mle_diagnostics()`: Produce standard regression diagnostics plots for MLE models
- `plot_vb_posteriors()`: Generate posterior predictive plots for Bayesian models
- `plot_vb_bayes_diagnostics()`: Comprehensive MCMC convergence and Bayesian model diagnostics
- `plot_vb_growth_pp_checks()`: Growth-specific posterior predictive checks for Bayesian models
- `plot_vb_predictions()`: Plot von Bertalanffy growth curves from prediction data with confidence/prediction intervals
- `plot_age_length_heatmap()`: Create heatmap visualisations of age-length observations with optional smoothers and grouping variables
- `plot_vb_age_counts()`: Plot frequency distribution of age samples by group variables
- `plot_empirical_cv()`: Plot empirical coefficient of variation by age for data exploration

## Key Enhancements in Version 0.2.0

### Statistical Testing

- **Bootstrap Permutation Tests**: Non-parametric testing framework for comparing growth parameters between groups
- **Age-Stratified Permutation**: Maintains age distributions during permutation for statistical tests
- **Growth Curve Comparison**: Tests differences in growth trajectories using maximum absolute deviance
- **Sex-Specific Model Support**: Handles comparisons of sex-specific growth parameters
- **Filtering**: Removal of groups/sex combinations with insufficient data

### Output and Visualisation

- **Summary Methods**: Formatted reports with statistical test results and method information
- **Density Plot Visualisation**: Shows bootstrap null distributions with observed statistics overlaid
- **Plotting**: Control over plot appearance, colours, and layout
- **Documentation**: Help with examples for functions

### Data Quality and Validation

- **Error Handling**: Input validation and error messages
- **Convergence Diagnostics**: Detection and handling of model fitting issues
- **Starting Value Optimisation**: Uses observed parameters as starting values for bootstrap convergence

## Example Usage

### Basic Usage

```r
library(growthVB)

# Example with simulated data (heteroscedastic errors)
set.seed(123)
true_Linf <- 120
true_k <- 0.25
true_t0 <- -0.5
true_cv <- 0.1

# Generate ages (multiple observations per age)
age <- rep(1:15, each = 6)
age <- age + rnorm(length(age), 0, 0.1)  # Add noise

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

# Access enhanced data with ID, fitted values, and residuals
head(fit$data)
# Columns: ID, age, length, fitted, residual

# Check that original order is preserved by ID
all(fit$data$ID == seq_along(age))  # Should be TRUE

# Examine residuals
summary(fit$data$residual)

# Plot growth curve
plot_vb(fit)

# Generate diagnostic plots
diagnostics <- plot_vb_mle_diagnostics(fit)

# Create age-length heatmap with smoother
heatmap_plot <- plot_age_length_heatmap(age = age, length = length, add_smoother = TRUE)

# Create heatmap with grouping variable (e.g., by year, study, site)
# group_var <- sample(c("2020", "2021", "2022"), length(age), replace = TRUE)
# heatmap_grouped <- plot_age_length_heatmap(age = age, length = length, 
#                                           group = group_var, add_smoother = TRUE)

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

# CV parameter
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

### Prediction Plotting

The `plot_vb_predictions()` function provides visualisation of von Bertalanffy growth curve predictions:

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

### Statistical Comparisons

### `compare_vb_mle()`

Compare von Bertalanffy growth parameters between groups using bootstrap permutation testing. This function evaluates whether growth parameters differ between two or more groups by comparing observed differences against a null distribution generated through random permutation of group labels.

**Features:**

- **Sex model support**: Works with both single and sex-specific von Bertalanffy models
- **Data quality filtering**: Removes groups/sex combinations with insufficient observations (default minimum 50)
- **Age-stratified permutation**: Preserves age distributions between groups for statistical testing (default 2-year age bins)
- **Bootstrap optimisation**: Uses observed parameters as starting values for convergence

```r
# Simulate data for two populations with different growth parameters
set.seed(123)
n_per_group <- 50
age1 <- runif(n_per_group, 1, 15)
age2 <- runif(n_per_group, 1, 15)

# Population A: Linf=100, k=0.2
length1 <- 100 * (1 - exp(-0.2 * (age1 - (-0.5)))) + rnorm(n_per_group, 0, 8)
# Population B: Linf=120, k=0.15 (different parameters)  
length2 <- 120 * (1 - exp(-0.15 * (age2 - (-0.5)))) + rnorm(n_per_group, 0, 10)

age <- c(age1, age2)
length <- c(length1, length2)
population <- rep(c("Pop_A", "Pop_B"), each = n_per_group)

# Perform bootstrap permutation test
comparison <- compare_vb_mle(
  age = age, 
  length = length, 
  group = population,
  n_bootstrap = 500,
  min_obs = 40,          # Minimum observations per group
  age_stratified = TRUE, # Use age-stratified permutation
  age_bin_width = 2,     # 2-year age bins for stratification
  seed = 123
)

# View results
print(comparison$p_values)
print(comparison$significant)
print(comparison$group_parameters)

# Use summary method for output
summary(comparison)

# Or simply print the object (uses print method)
comparison

# Visualise the permutation test results
plot(comparison)

# Plot specific parameters only
plot(comparison, parameters = c("Linf", "k"))

# Test with sex-specific models
sex <- rep(c("M", "F"), length.out = length(age))
comparison_sex <- compare_vb_mle(
  age = age, 
  length = length, 
  group = population,
  sex = sex,             # Sex-specific models
  parameters = c("Linf", "k"), # Test specific parameters
  n_bootstrap = 500,
  min_obs = 20,          # Lower threshold for sex combinations
  age_bin_width = 3,     # 3-year age bins
  seed = 123
)
```

The function returns:

- **observed_diffs**: Observed differences between groups for each parameter
- **p_values**: P-values for each parameter comparison
- **significant**: Logical indicating which parameters show differences
- **group_parameters**: Parameter estimates for each group
- **null_distributions**: Bootstrap null distributions for visualisation
- **curve_comparison**: Results from growth curve comparison test (if enabled)
- **method_info**: Test settings and sample information

**Summary and Visualisation Methods:**

The `summary()` method provides a formatted report including:

- Test method and configuration details
- Sample sizes and group information
- Parameter estimates for each group
- Statistical test results with p-values and significance indicators
- Growth curve comparison results (if performed)

Use `print(result)` or simply `result` for the same output, or `summary(result, digits = 3)` to control decimal precision.

**Visualisation:**

The `plot()` method creates density plots showing null distributions from the bootstrap permutation test:

```r
# Plot all parameters and curve comparison
plot(result)

# Plot specific parameters only
plot(result, parameters = c("Linf", "k"), include_curves = FALSE)

# Customise appearance
plot(result, observed_colour = "blue", alpha = 0.5, ncol = 3)
```

The plots show:

- Density distributions of the null hypothesis (no difference between groups)
- Red vertical lines indicating observed test statistics
- P-values and significance indicators in titles
- Growth curve comparison plots (if enabled)

**Function Parameters:**

- `age`, `length`, `group`: Vectors of age, length, and group classifications
- `sex`: Factor for sex-specific models (expands parameters to sex-specific versions)
- `n_bootstrap`: Number of bootstrap permutations (default 1000)
- `parameters`: Character vector of parameters to test (default: "Linf", "k", "t0", "CV")
- `min_obs`: Minimum observations per group/sex combination (default 50)
- `age_stratified`: Whether to use age-stratified permutation (default TRUE)
- `age_bin_width`: Width of age bins in years for stratification (default 2)
- `alpha`: Significance level (default 0.05)
- `verbose`: Whether to print progress messages (default TRUE)
- `seed`: Seed for reproducible results

Parameters can be selectively tested using the `parameters` argument. For sex models, parameters are expanded to sex-specific versions (e.g., "Linf_M", "Linf_F").

### Bayesian Models

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
  
  # Access Bayesian model data with fitted values and residuals
  head(fit_bayes$data)
  # Columns: ID, age, length, [sex], fitted, residual
  
  # Summary of Bayesian residuals (based on posterior mean fitted values)
  summary(fit_bayes$data$residual)
  
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

## Changelog

### Version 0.2.0 (November 2025)

**Enhancements:**

- **Statistical Testing Framework**: Added `compare_vb_mle()` for bootstrap permutation testing
- **Age-Stratified Permutation**: Preserves age distributions for statistical tests
- **Growth Curve Comparison**: Tests growth trajectory differences using maximum absolute deviance
- **Sex-Specific Model Support**: Support for comparing sex-specific growth parameters
- **Output Methods**: Added `summary()` and `print()` methods for formatted results
- **Visualisation Framework**: Added `plot()` method for bootstrap permutation test results
- **Documentation**: Help files and mathematical methodology documentation

**Technical Improvements:**

- Bootstrap optimisation using observed parameters as starting values
- Error handling and input validation
- Filtering of insufficient data combinations
- British English spelling throughout (colour, visualisation, etc.)
- Test coverage and validation

**New Functions:**

- `compare_vb_mle()`: Bootstrap permutation testing for group comparisons
- `summary.vb_comparison()`: Formatted output
- `print.vb_comparison()`: Simplified output method
- `plot.vb_comparison()`: Density plot visualisation of test results

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
