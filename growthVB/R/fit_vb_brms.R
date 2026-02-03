#' Fit von Bertalanffy Growth Function using Bayesian Methods (brms)
#'
#' This function fits a von Bertalanffy growth function to age and length data using
#' Bayesian methods via the brms package. The CV of the length is modeled as a function
#' of the mean length.
#'
#' @param age A numeric vector of ages
#' @param length A numeric vector of lengths
#' @param sex An optional factor or character vector specifying the sex for each observation
#' @param priors Optional brms prior object/list for the model parameters. If provided, this replaces defaults entirely.
#' @param prior_overrides Optional named character vector or list to selectively override default priors
#'   for specific non-linear parameters. Names must be in c("Linf","k","t0","CV"). Values are prior
#'   specification strings accepted by brms (e.g. "normal(120, 30)"). Ignored if `priors` is supplied.
#' @param fixed_params Optional named vector specifying parameters to fix at given values.
#'   Names should be from c("Linf", "k", "t0", "CV"). This sets a very tight prior (effectively fixing the parameter).
#'   Example: c(Linf = 100) to fix Linf at 100 while estimating k, t0, and CV.
#' @param chains Number of MCMC chains (default 4)
#' @param iter Number of iterations for each chain (default 4000)
#' @param cores Number of CPU cores to use for parallel processing (default: number of chains)
#' @param parallel_sex Logical. If TRUE and sex-specific models are being fit, run models for different sexes in parallel (default TRUE)
#' @param ... Additional arguments passed to brms::brm()
#'
#' @return A list containing:
#'   \item{models}{The brms model object(s)}
#'   \item{parameters}{Summary of parameter estimates}
#'   \item{predictions}{Predictions for plotting}
#'   \item{data}{Original data with fitted values and residuals, including ID column preserving input order}
#'
#' @examples
#' \dontrun{
#' # Simple example with simulated data
#' age <- 1:15
#' length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5)
#' fit <- fit_vb_brms(age = age, length = length)
#'
#' # Sex-specific models run in parallel
#' sex <- rep(c("M", "F"), length.out = 15)
#' fit_sex <- fit_vb_brms(age = age, length = length, sex = sex, parallel_sex = TRUE)
#'
#' # Example with selective prior overrides (modify specific parameters)
#' # Override only Linf and k priors, keep default t0 and CV priors
#' fit_override <- fit_vb_brms(
#'   age = age, length = length,
#'   prior_overrides = list(
#'     Linf = "normal(120, 20)", # More informative prior for Linf
#'     k = "gamma(2, 10)" # Gamma prior for k (ensures positive)
#'   )
#' )
#'
#' # Example with completely custom priors (replaces all defaults)
#' library(brms)
#' custom_priors <- c(
#'   prior_string("normal(110, 15)", nlpar = "Linf", lb = 0),
#'   prior_string("gamma(3, 15)", nlpar = "k", lb = 0),
#'   prior_string("normal(-0.5, 0.5)", nlpar = "t0"),
#'   prior_string("gamma(2, 20)", nlpar = "CV", lb = 0)
#' )
#' fit_custom <- fit_vb_brms(
#'   age = age, length = length,
#'   priors = custom_priors
#' )
#'
#' # Example combining sex-specific models with custom priors
#' fit_sex_custom <- fit_vb_brms(
#'   age = age, length = length, sex = sex,
#'   prior_overrides = list(
#'     Linf = "normal(130, 25)", # Expect larger maximum length
#'     CV = "gamma(1.5, 15)" # More constrained CV
#'   ),
#'   parallel_sex = TRUE
#' )
#' }
#'
#' @importFrom stats predict
#' @importFrom brms brm fixef
#' @export
fit_vb_brms <- function(age, length, sex = NULL,
                        priors = NULL, prior_overrides = NULL,
                        fixed_params = NULL,
                        chains = 4, iter = 4000, cores = chains,
                        parallel_sex = TRUE, ...) {
  # Check if brms is installed
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required for this function. Please install it.",
      call. = FALSE
    )
  }

  # Create data frame for analysis with original order index
  data <- data.frame(
    ID = seq_along(age),  # Add unique identifier preserving original order
    age = age, 
    length = length
  )

  # Add sex if provided
  if (!is.null(sex)) {
    data$sex <- sex
    # Filter out any NA values but preserve ID column
    data <- data[!is.na(data$age) & !is.na(data$length) & !is.na(data$sex), ]
  } else {
    # Filter out any NA values but preserve ID column
    data <- data[!is.na(data$age) & !is.na(data$length), ]
  }

  # Build priors: explicit 'priors' takes precedence; otherwise use defaults with optional overrides
  if (is.null(priors)) {
    make_default_priors <- function() {
      list(
        Linf = brms::prior_string("normal(150, 100)", nlpar = "Linf", lb = 0),
        k    = brms::prior_string("normal(0.1, 100)", nlpar = "k", lb = 0),
        t0   = brms::prior_string("normal(0, 1)", nlpar = "t0"),
        CV   = brms::prior_string("normal(0.1, 10)", nlpar = "CV", lb = 0)
      )
    }

    # Start from defaults
    prior_list <- make_default_priors()

    # Handle fixed parameters by setting very tight priors
    if (!is.null(fixed_params)) {
      # Validate fixed parameter names
      allowed_names <- c("Linf", "k", "t0", "CV")
      if (!all(names(fixed_params) %in% allowed_names)) {
        invalid_names <- setdiff(names(fixed_params), allowed_names)
        stop("Invalid parameter names in fixed_params: ", paste(invalid_names, collapse = ", "),
             ". Must be from: Linf, k, t0, CV")
      }
      
      # Set very tight priors (sd = 0.001) for fixed parameters
      for (param_name in names(fixed_params)) {
        param_value <- fixed_params[[param_name]]
        # Create a very tight normal prior centered at the fixed value
        lb <- if (param_name %in% c("Linf", "k", "CV")) 0 else NULL
        prior_string <- sprintf("normal(%s, 0.001)", param_value)
        prior_list[[param_name]] <- brms::prior_string(prior_string, nlpar = param_name, lb = lb)
      }
    }

    # Apply selective overrides if provided (these override defaults but not fixed_params)
    if (!is.null(prior_overrides)) {
      if (is.list(prior_overrides)) {
        # convert to character vector preserving names
        prior_overrides <- unlist(prior_overrides)
      }
      if (!is.null(names(prior_overrides))) {
        allowed <- intersect(names(prior_overrides), c("Linf", "k", "t0", "CV"))
        for (nm in allowed) {          # Skip if this parameter is fixed (fixed_params take precedence)
          if (!is.null(fixed_params) && nm %in% names(fixed_params)) {
            warning(sprintf("Parameter %s is fixed. Ignoring prior_override for this parameter.", nm))
            next
          }          # lower bounds for parameters that must be >= 0
          lb <- if (nm %in% c("Linf", "k", "CV")) 0 else NULL
          # Construct a prior for this parameter from the string
          prior_list[[nm]] <- brms::prior_string(as.character(prior_overrides[[nm]]), nlpar = nm, lb = lb)
        }
      }
    }
    # Collapse to a vector/list that brms::brm accepts
    priors <- do.call(c, unname(prior_list))
  }

  # Function to fit a Bayesian model
  fit_brms_model <- function(subset_data, cores_per_model = cores) {
    # von Bertalanffy model with CV dependent on mean length
    model <- brms::brm(
      brms::bf(length ~ eta, nl = TRUE) +
        brms::nlf(eta ~ 1 + Linf * (1.0 - exp(-k * (age - t0)))) +
        brms::nlf(sigma ~ eta * CV) +
        brms::lf(Linf ~ 1, k ~ 1, t0 ~ 1, CV ~ 1),
      data = subset_data,
      chains = chains,
      cores = cores_per_model,
      prior = priors,
      family = brms::brmsfamily("gaussian", link_sigma = "identity"),
      iter = iter,
      ...
    )

    return(model)
  }

  # Results list
  results <- list()

  # Split data by sex if provided, otherwise fit single model
  if (!is.null(sex) && length(unique(sex[!is.na(sex)])) > 1) {
    sex_levels <- unique(sex[!is.na(sex)])
    models_list <- list()
    preds_list <- list()
    data_list <- list()

    if (parallel_sex && length(sex_levels) > 1) {
      # Parallel processing for sex-specific models

      # Check if parallel package is available
      if (!requireNamespace("parallel", quietly = TRUE)) {
        warning("Package 'parallel' not available. Running sex-specific models sequentially.")
        parallel_sex <- FALSE
      } else {
        # Determine cores allocation
        n_sexes <- length(sex_levels)

        # If we have enough cores, allocate cores per sex model
        # Otherwise, use 1 core per sex and reduce chains per model
        if (cores >= n_sexes * chains) {
          cores_per_sex <- cores %/% n_sexes
        } else if (cores >= n_sexes) {
          cores_per_sex <- 1
          chains <- min(chains, cores %/% n_sexes)
          if (chains < 2) {
            warning("Insufficient cores for parallel sex fitting with multiple chains. Consider increasing cores or setting parallel_sex=FALSE.")
            parallel_sex <- FALSE
          }
        } else {
          warning("Insufficient cores for parallel sex fitting. Running sequentially.")
          parallel_sex <- FALSE
        }

        if (parallel_sex) {
          message(sprintf(
            "Fitting %d sex-specific models in parallel using %d total cores (%d cores per model, %d chains each)",
            n_sexes, cores, cores_per_sex, chains
          ))

          # Create cluster
          cl <- parallel::makeCluster(n_sexes)

          # Ensure cluster is closed on exit
          on.exit(parallel::stopCluster(cl), add = TRUE)

          # Export necessary objects to cluster
          parallel::clusterEvalQ(cl, {
            library(brms)
          })

          parallel::clusterExport(cl, c("fit_brms_model", "priors", "chains", "iter", "cores_per_sex"),
            envir = environment()
          )

          # Prepare data for each sex
          sex_data_list <- lapply(sex_levels, function(s) {
            list(sex = s, data = data[data$sex == s, ])
          })
          names(sex_data_list) <- sex_levels

          # Fit models in parallel
          fit_one_sex <- function(sex_info) {
            model <- fit_brms_model(sex_info$data, cores_per_model = cores_per_sex)

            # Generate predictions for plotting
            age_seq <- sort(unique(sex_info$data$age))
            new_data <- data.frame(age = age_seq)
            preds <- stats::predict(model, newdata = new_data, probs = c(0.025, 0.975))

            preds_df <- cbind(as.data.frame(preds), new_data)
            preds_df$Sex <- sex_info$sex
            preds_df$Model <- "von Bertalanffy"

            return(list(model = model, predictions = preds_df))
          }

          # Run parallel fitting
          parallel_results <- parallel::parLapply(cl, sex_data_list, fit_one_sex)

          # Extract results
          models_list <- lapply(parallel_results, function(x) x$model)
          names(models_list) <- sex_levels

          preds_list <- lapply(parallel_results, function(x) x$predictions)
          names(preds_list) <- sex_levels
        }
      }
    }

    # Sequential processing (fallback or if parallel_sex = FALSE)
    if (!parallel_sex) {
      message(sprintf("Fitting %d sex-specific models sequentially", length(sex_levels)))

      for (s in sex_levels) {
        subset_data <- data[data$sex == s, ]
        model <- fit_brms_model(subset_data)
        models_list[[s]] <- model

        # Calculate fitted values and residuals for this sex
        fitted_values <- brms::fitted(model)[, "Estimate"]  # Extract fitted values
        residuals <- subset_data$length - fitted_values
        
        # Add fitted values and residuals to subset data
        subset_data$fitted <- fitted_values
        subset_data$residual <- residuals
        data_list[[s]] <- subset_data

        # Generate predictions for plotting
        age_seq <- sort(unique(subset_data$age))
        new_data <- data.frame(age = age_seq)
        preds <- stats::predict(model, newdata = new_data, probs = c(0.025, 0.975))

        preds_df <- cbind(as.data.frame(preds), new_data)
        preds_df$Sex <- s
        preds_df$Model <- "von Bertalanffy"

        preds_list[[s]] <- preds_df
      }
    } else {
      # Parallel branch: compute fitted values/residuals in main process
      for (s in sex_levels) {
        subset_data <- data[data$sex == s, ]
        model <- models_list[[s]]
        fitted_values <- brms::fitted(model)[, "Estimate"]
        residuals <- subset_data$length - fitted_values

        subset_data$fitted <- fitted_values
        subset_data$residual <- residuals
        data_list[[s]] <- subset_data
      }
    }

    # Get parameter summaries
    params <- lapply(models_list, function(m) brms::fixef(m))

    # Combine data with fitted values and residuals, preserving order by ID
    combined_data <- do.call(rbind, data_list)

    # Store results with data ordered by ID
    results$models <- models_list
    results$parameters <- params
    results$predictions <- preds_list
    results$data <- combined_data[order(combined_data$ID), ]  # Include data with fitted/residuals ordered by ID
  } else {
    # Fit a single model to all data
    message("Fitting single combined model")
    model <- fit_brms_model(data)

    # Generate predictions for plotting
    age_seq <- sort(unique(data$age))
    new_data <- data.frame(age = age_seq)
    preds <- stats::predict(model, newdata = new_data, probs = c(0.025, 0.975))

    preds_df <- cbind(as.data.frame(preds), new_data)
    preds_df$Model <- "von Bertalanffy"

    # Calculate fitted values and residuals for the original data
    fitted_values <- brms::fitted(model)[, "Estimate"]  # Extract fitted values
    residuals <- data$length - fitted_values
    
    # Add fitted values and residuals to data
    data$fitted <- fitted_values
    data$residual <- residuals

    # Store results with data ordered by ID
    results$models <- model
    results$parameters <- brms::fixef(model)
    results$predictions <- preds_df
    results$data <- data[order(data$ID), ]  # Include original data ordered by ID
  }

  class(results) <- c("vb_brms", "list")
  return(results)
}
