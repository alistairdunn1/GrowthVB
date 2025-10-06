# Load the growthVB package
library(growthVB)

# Create a simple testing function for plot_vb_age_counts
test_plot_vb_age_counts <- function() {
  # Generate test data
  set.seed(123)
  ages <- c(1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5)
  years <- c(2020, 2020, 2021, 2020, 2021, 2022, 2020, 2021, 2021, 2022, 2020, 2021, 2021, 2022, 2022)
  sexes <- c("M", "M", "F", "M", "F", "M", "F", "F", "M", "M", "F", "F", "M", "M", "F")
  
  # Test 1: Basic plot
  cat("Test 1: Basic age counts plot\n")
  p1 <- plot_vb_age_counts(age = ages)
  print(p1)
  
  # Test 2: With year
  cat("\nTest 2: Age counts by year\n")
  p2 <- plot_vb_age_counts(age = ages, year = years)
  print(p2)
  
  # Test 3: With sex
  cat("\nTest 3: Age counts by sex\n")
  p3 <- plot_vb_age_counts(age = ages, sex = sexes)
  print(p3)
  
  # Test 4: With year and sex
  cat("\nTest 4: Age counts by year and sex\n")
  p4 <- plot_vb_age_counts(age = ages, year = years, sex = sexes)
  print(p4)
  
  cat("\nAll tests passed!\n")
}

# Run the tests
test_plot_vb_age_counts()
