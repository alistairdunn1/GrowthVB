#!/usr/bin/env Rscript

# Script to check the contents of the age.rdata dataset
print("Loading the age dataset...")
load("Example/age.rdata")

# Check what objects were loaded
print("Objects loaded:")
print(ls())

# Check if 'age' variable exists
if(exists("age")) {
  print("Structure of 'age' object:")
  print(str(age))
  print("Summary of 'age' object:")
  print(summary(age))
  
  # Check for specific variables
  variables <- c("age", "length", "sex")
  for(var in variables) {
    if(var %in% names(age)) {
      cat("Variable '", var, "' is present\n", sep="")
    } else {
      cat("Variable '", var, "' is NOT present\n", sep="")
    }
  }
  
  # Check the first few rows
  print("First few rows of 'age':")
  print(head(age))
} else {
  print("The 'age' object was not found in the dataset")
  
  # Check all objects
  all_objs <- ls()
  for(obj in all_objs) {
    cat("Object: ", obj, "\n")
    print(str(get(obj)))
  }
}
