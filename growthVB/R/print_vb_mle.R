#' Print method for vb_mle objects
#'
#' @param x A fitted model object of class vb_mle
#' @param digits Number of digits to display (default 4)
#' @param ... Additional arguments (not used)
#'
#' @return The object x, invisibly
#'
#' @export
print.vb_mle <- function(x, digits = 4, ...) {
  cat("Von Bertalanffy Growth Model (Maximum Likelihood Estimation)\n\n")

  # Check if we have multiple models (sex-specific)
  if (is.list(x$model) && !inherits(x$model, "vb_optim")) {
    cat("Sex-specific models:\n")
    for (s in names(x$model)) {
      cat("\n", s, ":", "\n", sep = "")
      print_single_model(x$model[[s]], digits)
    }
  } else {
    # Single model
    print_single_model(x$model, digits)
  }

  invisible(x)
}

#' Helper function for printing model results
#'
#' @param model A single vb_optim model
#' @param digits Number of digits to display
#'
#' @keywords internal
print_single_model <- function(model, digits) {
  # Extract parameters and standard errors
  params <- model$parameters
  se <- model$se

  # Create parameter table
  param_table <- data.frame(
    Parameter = names(params),
    Estimate = params,
    row.names = NULL
  )

  # Add standard errors if available
  if (!is.null(se)) {
    param_table$Std.Error <- se

    # Calculate z/t values
    param_table$z.value <- param_table$Estimate / param_table$Std.Error

    # Calculate p-values (using normal approximation)
    param_table$p.value <- 2 * stats::pnorm(-abs(param_table$z.value))

    # Add significance stars
    param_table$Signif <- symnum(param_table$p.value,
      corr = FALSE,
      na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols = c("***", "**", "*", ".", " ")
    )
  }

  # Print parameter table
  print(format(param_table, digits = digits), row.names = FALSE)

  # Print significance codes if we have p-values
  if (!is.null(se)) {
    cat("---\n")
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }

  # Print convergence info
  cat("\nConvergence info:\n")
  cat("  Log-likelihood:", -model$value, "\n")
  cat(
    "  Convergence code:", model$convergence,
    if (model$convergence == 0) "(successful)" else "(check convergence!)", "\n"
  )
  cat("  Iterations:", model$counts[1], "\n")

  if (!is.null(model$message)) {
    cat("  Message:", model$message, "\n")
  }
}
