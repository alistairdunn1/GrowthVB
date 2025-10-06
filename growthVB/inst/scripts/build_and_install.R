# This script helps with building and installing the growthVB package

# First, set the working directory to the package root
# If running from the GrowthVB project root, navigate to the growthVB subdirectory:
# setwd("growthVB")
# If running from elsewhere, set the full path:
# setwd("C:/path/to/GrowthVB/growthVB")

# Build the documentation with roxygen2
if (requireNamespace("roxygen2", quietly = TRUE)) {
  roxygen2::roxygenize()
} else {
  warning("Package roxygen2 is required to build documentation")
}

# Build and install the package
install.packages(".", repos = NULL, type = "source")

# Run tests
if (requireNamespace("testthat", quietly = TRUE)) {
  testthat::test_package("growthVB")
} else {
  warning("Package testthat is required to run tests")
}

# Check the package
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::check()
} else {
  warning("Package devtools is required to run R CMD check")
}
