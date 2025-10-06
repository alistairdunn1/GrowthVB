# Test installation instructions for growthVB package
# Run this script from the GrowthVB project root directory

# Check current working directory
cat("Current working directory:", getwd(), "\n")

# Check if growthVB subdirectory exists
if (dir.exists("growthVB")) {
  cat("✓ growthVB subdirectory found\n")

  # Check if DESCRIPTION file exists
  if (file.exists("growthVB/DESCRIPTION")) {
    cat("✓ DESCRIPTION file found in growthVB/\n")

    # Try to install the package
    tryCatch(
      {
        install.packages("growthVB", repos = NULL, type = "source")
        cat("✓ Package installation completed\n")
      },
      error = function(e) {
        cat("✗ Installation failed:", e$message, "\n")
      }
    )
  } else {
    cat("✗ DESCRIPTION file not found in growthVB/\n")
  }
} else {
  cat("✗ growthVB subdirectory not found\n")
  cat("Available directories:\n")
  print(list.dirs(".", recursive = FALSE))
}
