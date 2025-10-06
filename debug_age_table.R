#!/usr/bin/env Rscript

# Debug test script
print("Loading AGE dataset...")
load("Example/age.rdata")
print("Dataset loaded")

# Check age variable
print("AGE$age class:")
print(class(AGE$age))

print("First few age values:")
print(head(AGE$age))

# Try direct table call
print("Testing table function directly:")
age_table <- table(AGE$age)
print("Table created successfully")
print("Table structure:")
print(str(age_table))
print("Table head:")
print(head(age_table))

# Test converting names to numeric
print("Testing conversion of table names to numeric:")
names_vector <- as.numeric(names(age_table))
print("First few converted names:")
print(head(names_vector))

# Create counts data frame directly
print("Creating counts data frame:")
counts <- data.frame(
  age = as.numeric(names(age_table)),
  count = as.numeric(age_table)
)
print("Counts data frame created")
print("Structure:")
print(str(counts))
print("Head:")
print(head(counts))

# Try plotting
print("Testing basic ggplot:")
library(ggplot2)
p <- tryCatch({
  ggplot(counts, aes(x = age, y = count)) +
    geom_col(fill = "steelblue") +
    labs(x = "Age", y = "Count") +
    theme_minimal()
  print("Plot created successfully")
}, error = function(e) {
  print(paste("Error creating plot:", e$message))
  return(NULL)
})

if(!is.null(p)) print("Plot creation successful")
