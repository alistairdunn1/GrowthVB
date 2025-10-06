library(growthVB)
library(testthat)

# Test basic functionality without fitting
test_that("Basic package functions exist", {
  expect_true(exists("fit_vb_nls"))
  expect_true(exists("summarize_vb"))
  expect_true(exists("plot_vb"))
  expect_true(exists("plot_vb_diagnostics"))
  expect_true(exists("plot_age_length_heatmap"))
  expect_true(exists("plot_vb_age_counts"))
})

# Test age-length heatmap function directly
test_that("plot_age_length_heatmap works", {
  # Create test data
  set.seed(123)
  ages <- rep(1:10, each = 5)
  lengths <- 50 + ages * 5 + rnorm(length(ages), 0, 2)

  # Check that the function runs without error
  expect_error(p <- plot_age_length_heatmap(age = ages, length = lengths), NA)
  expect_s3_class(p, "gg")
})

# Test age count summary
test_that("plot_vb_age_counts works", {
  # Create test data
  set.seed(123)
  ages <- rep(1:10, each = 5)
  years <- sample(2020:2023, length(ages), replace = TRUE)
  sexes <- sample(c("M", "F"), length(ages), replace = TRUE)

  # Test with age only
  expect_error(p1 <- plot_vb_age_counts(age = ages), NA)
  expect_s3_class(p1, "gg")

  # Test with age and year
  expect_error(p2 <- plot_vb_age_counts(age = ages, year = years), NA)
  expect_s3_class(p2, "gg")

  # Test with age and sex
  expect_error(p3 <- plot_vb_age_counts(age = ages, sex = sexes), NA)
  expect_s3_class(p3, "gg")

  # Test with age, year, and sex
  expect_error(p4 <- plot_vb_age_counts(age = ages, year = years, sex = sexes), NA)
  expect_s3_class(p4, "gg")
})
