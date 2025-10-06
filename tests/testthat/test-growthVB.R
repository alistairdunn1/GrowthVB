library(growthVB)

test_that("fit_vb_nls works correctly", {
  # Create a very clean dataset for testing
  set.seed(123)
  age <- seq(1, 20, by = 1)
  # Generate perfect von Bertalanffy data with no noise
  true_Linf <- 100
  true_k <- 0.2
  true_t0 <- -0.5
  # Very small noise to ensure convergence
  length <- true_Linf * (1 - exp(-true_k * (age - true_t0))) + rnorm(length(age), 0, 0.1)

  # Fit the model
  fit <- fit_vb_nls(age = age, length = length) # Check that the fit has the correct structure
  expect_s3_class(fit, "vb_nls")
  expect_true(all(c("parameters", "data", "fits", "model") %in% names(fit)))

  # Check that parameters are reasonably close to the true values
  expect_true(abs(fit$parameters$Linf - 100) < 20)
  expect_true(abs(fit$parameters$k - 0.2) < 0.1)
  expect_true(abs(fit$parameters$t0 - (-0.5)) < 0.5)
})

test_that("summarize_vb works correctly", {
  # Create a very clean dataset for testing
  set.seed(123)
  age <- seq(1, 20, by = 1)
  # Generate perfect von Bertalanffy data with no noise
  true_Linf <- 100
  true_k <- 0.2
  true_t0 <- -0.5
  # Very small noise to ensure convergence
  length <- true_Linf * (1 - exp(-true_k * (age - true_t0))) + rnorm(length(age), 0, 0.1)

  # Fit the model
  fit <- fit_vb_nls(age = age, length = length) # Summarize
  summary <- summarize_vb(fit)

  # Check structure
  expect_true(all(c("Parameter", "Estimate", "Lower_CI", "Upper_CI") %in% names(summary)))
  expect_equal(nrow(summary), 3) # Linf, k, t0
})

test_that("plotting functions don't error", {
  # Create a very clean dataset for testing
  set.seed(123)
  # Generate data for both sexes
  age_M <- seq(1, 20, by = 1)
  age_F <- seq(1, 20, by = 1)
  age <- c(age_M, age_F)

  # Generate perfect von Bertalanffy data with different params by sex
  true_Linf_M <- 80
  true_k_M <- 0.3
  true_t0_M <- -0.2

  true_Linf_F <- 100
  true_k_F <- 0.2
  true_t0_F <- -0.5

  # Very small noise to ensure convergence
  length_M <- true_Linf_M * (1 - exp(-true_k_M * (age_M - true_t0_M))) + rnorm(length(age_M), 0, 0.1)
  length_F <- true_Linf_F * (1 - exp(-true_k_F * (age_F - true_t0_F))) + rnorm(length(age_F), 0, 0.1)
  length <- c(length_M, length_F)

  sex <- c(rep("M", length(age_M)), rep("F", length(age_F)))

  # Fit the model
  fit <- fit_vb_nls(age = age, length = length, sex = sex) # Check that plots can be created without error
  expect_error(plot_vb(fit), NA)
  expect_error(plot_vb_diagnostics(fit), NA)
  expect_error(plot_age_length_heatmap(age, length, sex), NA)
  expect_error(plot_vb_age_counts(age, sex = sex), NA)
})
