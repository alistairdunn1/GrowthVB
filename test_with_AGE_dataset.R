#!/usr/bin/env Rscript

# Test script for plot_vb_age_counts function using the AGE dataset
library(growthVB)
library(ggplot2)

# Load the AGE dataset
print("Loading AGE dataset...")
load("Example/age.rdata")
print("Dataset loaded")

# Verify the dataset structure
print("AGE dataset structure:")
str(AGE)

# Function to test all scenarios
test_with_AGE <- function() {
  # Test 1: Basic plot
  print("Test 1: Basic age counts plot")
  p1 <- tryCatch({
    plot_vb_age_counts(age = AGE$age)
  }, error = function(e) {
    print(paste("Error in basic plot:", e$message))
    return(NULL)
  })
  if(!is.null(p1)) print("Basic plot completed successfully")
  
  # Test 2: With sex
  print("Test 2: Age counts by sex")
  p2 <- tryCatch({
    plot_vb_age_counts(age = AGE$age, sex = AGE$sex)
  }, error = function(e) {
    print(paste("Error in plot with sex:", e$message))
    return(NULL)
  })
  if(!is.null(p2)) print("Plot with sex completed successfully")
  
  # Test 3: With length (as a placeholder for year)
  print("Test 3: Age counts by length groups")
  # Create length groups as a stand-in for year
  length_groups <- cut(AGE$length, breaks=5)
  p3 <- tryCatch({
    plot_vb_age_counts(age = AGE$age, year = length_groups)
  }, error = function(e) {
    print(paste("Error in plot with length groups:", e$message))
    return(NULL)
  })
  if(!is.null(p3)) print("Plot with length groups completed successfully")
  
  # Test 4: With both
  print("Test 4: Age counts by length groups and sex")
  p4 <- tryCatch({
    plot_vb_age_counts(age = AGE$age, year = length_groups, sex = AGE$sex)
  }, error = function(e) {
    print(paste("Error in plot with length groups and sex:", e$message))
    return(NULL)
  })
  if(!is.null(p4)) print("Plot with length groups and sex completed successfully")
}

# Run the tests
test_with_AGE()
print("Tests completed")
