#!/usr/bin/env Rscript

# Script to build and reinstall the package
print("Building and installing the package...")

# Set working directory to the package
setwd("growthVB")

# Build and install
system("R CMD INSTALL .")

# Return to original directory
setwd("..")

print("Package installed successfully!")

# Test if we can use the updated function
library(growthVB)

# Check if the plot_vb_age_counts function is defined
if(exists("plot_vb_age_counts", envir = asNamespace("growthVB"))) {
  print("Function plot_vb_age_counts exists in the package")
} else {
  print("Function plot_vb_age_counts does NOT exist in the package")
}

# Test with the simplest case
print("Testing with simple data...")
simple_ages <- 1:10
tryCatch({
  p <- plot_vb_age_counts(simple_ages)
  print("Simple test successful!")
}, error = function(e) {
  print(paste("Error in simple test:", e$message))
})

# Test with AGE data
print("Testing with AGE data...")
load("Example/age.rdata")
tryCatch({
  p <- plot_vb_age_counts(AGE$age)
  print("AGE test successful!")
}, error = function(e) {
  print(paste("Error in AGE test:", e$message))
})
