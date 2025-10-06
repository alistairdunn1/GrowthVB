#' Example Usage of growthVB Package
#'
#' This script demonstrates how to use the growthVB package to fit and
#' analyze von Bertalanffy growth curves.

# Load required packages
library(growthVB)
library(ggplot2)

# Generate simulated data
set.seed(123)
n_samples <- 200

# Parameters for males
Linf_M <- 120
k_M <- 0.25
t0_M <- -0.5

# Parameters for females
Linf_F <- 150
k_F <- 0.18
t0_F <- -0.3

# Generate balanced dataset
age_M <- sample(1:20, n_samples / 2, replace = TRUE)
age_F <- sample(1:20, n_samples / 2, replace = TRUE)
age <- c(age_M, age_F)

# Generate lengths with CV dependent on mean length
mean_length_M <- Linf_M * (1 - exp(-k_M * (age_M - t0_M)))
mean_length_F <- Linf_F * (1 - exp(-k_F * (age_F - t0_F)))

# CV increases with length
cv_M <- 0.05 + 0.01 * sqrt(mean_length_M)
cv_F <- 0.05 + 0.01 * sqrt(mean_length_F)

length_M <- rnorm(n_samples / 2, mean_length_M, mean_length_M * cv_M)
length_F <- rnorm(n_samples / 2, mean_length_F, mean_length_F * cv_F)
length <- c(length_M, length_F)

sex <- factor(c(rep("M", n_samples / 2), rep("F", n_samples / 2)))
year <- sample(2020:2023, n_samples, replace = TRUE)

# Create a data frame
data <- data.frame(age = age, length = length, sex = sex, year = year)

# 1. Fit von Bertalanffy model using NLS
vb_fit <- fit_vb_nls(age = data$age, length = data$length, sex = data$sex)

# 2. Summarize parameters
params <- summarize_vb(vb_fit)
print(params)

# 3. Plot the growth curves
p1 <- plot_vb(vb_fit)
print(p1)

# 4. Plot diagnostics
diagnostics <- plot_vb_diagnostics(vb_fit)
print(diagnostics$M$residuals_vs_fitted)
print(diagnostics$F$qq_plot)

# 5. Plot age-length heatmap
p2 <- plot_age_length_heatmap(age = data$age, length = data$length, sex = data$sex)
print(p2)

# 6. Plot age counts by year
p3 <- plot_vb_age_counts(age = data$age, year = data$year, sex = data$sex)
print(p3)

# 7. Fit Bayesian model (only run if brms is installed)
if (requireNamespace("brms", quietly = TRUE)) {
  cat("Fitting Bayesian model with brms. This may take some time...\n")

  # Use a subset of data to speed up the example
  sample_indices <- sample(1:nrow(data), 50)
  sample_data <- data[sample_indices, ]

  # Fit Bayesian model
  vb_brms <- fit_vb_brms(
    age = sample_data$age,
    length = sample_data$length,
    sex = sample_data$sex,
    iter = 2000 # Smaller number of iterations for the example
  )

  # Summarize parameters
  brms_params <- summarize_vb(vb_brms)
  print(brms_params)

  # Plot growth curves
  p4 <- plot_vb(vb_brms)
  print(p4)

  # Plot posterior predictive checks
  if (requireNamespace("bayesplot", quietly = TRUE)) {
    post_plots <- plot_vb_posteriors(vb_brms)
    print(post_plots$M$posterior_predictive)
    print(post_plots$F$parameter_densities)
  }
}
