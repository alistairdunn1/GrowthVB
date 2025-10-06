#!/usr/bin/env Rscript

# Debug script to trace execution of plot_vb_age_counts line by line
print("Loading AGE dataset...")
load("Example/age.rdata")
print("Dataset loaded")

print("AGE dataset structure:")
str(AGE)

print("Running debug version of plot_vb_age_counts...")

# Debug by showing each step
debug_plot <- function() {
  # Get the age data
  age_data <- AGE$age
  print("Age data extracted:")
  print(class(age_data))
  print(head(age_data))
  
  # Convert to numeric
  age_numeric <- as.numeric(age_data)
  print("Converted to numeric:")
  print(class(age_numeric))
  print(head(age_numeric))
  
  # Create table
  print("Creating table...")
  age_table <- table(age_numeric)
  print("Table created:")
  print(head(age_table))
  print(class(age_table))
  
  # Extract names and counts
  print("Extracting table components...")
  table_names <- names(age_table)
  print("Table names:")
  print(head(table_names))
  
  age_values <- as.numeric(table_names)
  print("Age values:")
  print(head(age_values))
  
  count_values <- as.numeric(age_table)
  print("Count values:")
  print(head(count_values))
  
  # Create data frame
  print("Creating data frame...")
  counts_df <- data.frame(
    age = age_values,
    count = count_values
  )
  print("Data frame created:")
  print(head(counts_df))
  
  # Create plot
  print("Creating plot...")
  p <- ggplot2::ggplot(counts_df, ggplot2::aes(x = age, y = count)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::labs(x = "Age", y = "Count") +
    ggplot2::theme_minimal()
  
  print("Plot created successfully!")
  return(p)
}

# Run the debug function
tryCatch({
  debug_plot()
  print("Success!")
}, error = function(e) {
  print(paste("Error at:", e$message))
  print(paste("In:", e$call))
})
