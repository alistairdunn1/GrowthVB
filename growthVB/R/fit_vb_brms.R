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
#'   for specific non-linear parameters. Names must be in c("Linf","k","t0","tau"). Values are prior
#'   specification strings accepted by brms (e.g. "normal(120, 30)"). Ignored if `priors` is supplied.
#' @param chains Number of MCMC chains (default 4)
#' @param iter Number of iterations for each chain (default 4000)
#' @param ... Additional arguments passed to brms::brm()
#'
#' @return A list containing:
#'   \item{models}{The brms model object(s)}
#'   \item{parameters}{Summary of parameter estimates}
#'   \item{predictions}{Predictions for plotting}
#'
#' @examples
#' \dontrun{
#' # Simple example with simulated data
#' age <- 1:15
#' length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5)
#' fit <- fit_vb_brms(age = age, length = length)
#' }
#'
#' @importFrom stats predict
#' @importFrom brms brm fixef
#' @export
fit_vb_brms <- function(age, length, sex = NULL,
                        priors = NULL, prior_overrides = NULL,
                        chains = 4, iter = 4000, ...) {
  # Check if brms is installed
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required for this function. Please install it.",
      call. = FALSE
    )
  }

  # Create data frame for analysis
  data <- data.frame(age = age, length = length)

  # Add sex if provided
  if (!is.null(sex)) {
    data$sex <- sex
    # Filter out any NA values
    data <- data[!is.na(data$age) & !is.na(data$length) & !is.na(data$sex), ]
  } else {
    # Filter out any NA values
    data <- data[!is.na(data$age) & !is.na(data$length), ]
  }

  # Build priors: explicit 'priors' takes precedence; otherwise use defaults with optional overrides
  if (is.null(priors)) {
    make_default_priors <- function() {
      list(
        Linf = brms::prior_string("normal(150, 100)", nlpar = "Linf", lb = 0),
        k    = brms::prior_string("normal(0.3, 100)", nlpar = "k", lb = 0),
        t0   = brms::prior_string("normal(0, 1)", nlpar = "t0"),
        tau  = brms::prior_string("normal(0, 100)", nlpar = "tau", lb = 0)
      )
    }

    # Start from defaults
    prior_list <- make_default_priors()

    # Apply selective overrides if provided
    if (!is.null(prior_overrides)) {
      if (is.list(prior_overrides)) {
        # convert to character vector preserving names
        prior_overrides <- unlist(prior_overrides)
      }
      if (!is.null(names(prior_overrides))) {
        allowed <- intersect(names(prior_overrides), c("Linf", "k", "t0", "tau"))
        for (nm in allowed) {
          # lower bounds for parameters that must be >= 0
          lb <- if (nm %in% c("Linf", "k", "tau")) 0 else NULL
          # Construct a prior for this parameter from the string
          prior_list[[nm]] <- brms::prior_string(as.character(prior_overrides[[nm]]), nlpar = nm, lb = lb)
        }
      }
    }
    # Collapse to a vector/list that brms::brm accepts
    priors <- do.call(c, unname(prior_list))
  }

  # Function to fit a Bayesian model
  fit_brms_model <- function(subset_data) {
    # von Bertalanffy model with CV dependent on mean length
    model <- brms::brm(
      brms::bf(length ~ eta, nl = TRUE) +
        brms::nlf(eta ~ 1 + Linf * (1.0 - exp(-k * (age - t0)))) +
        brms::nlf(sigma ~ eta * tau) +
        brms::lf(Linf ~ 1, k ~ 1, t0 ~ 1, tau ~ 1),
      data = subset_data,
      chains = chains,
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
    # Process each sex separately
    sex_levels <- unique(sex[!is.na(sex)])
    models_list <- list()
    preds_list <- list()

    for (s in sex_levels) {
      subset_data <- data[data$sex == s, ]
      model <- fit_brms_model(subset_data)
      models_list[[s]] <- model

      # Generate predictions for plotting
      age_seq <- sort(unique(subset_data$age))
      new_data <- data.frame(age = age_seq)
  preds <- stats::predict(model, newdata = new_data, probs = c(0.025, 0.975))

      preds_df <- cbind(as.data.frame(preds), new_data)
      preds_df$Sex <- s
      preds_df$Model <- "von Bertalanffy"

      preds_list[[s]] <- preds_df
    }

    # Get parameter summaries
    params <- lapply(models_list, function(m) brms::fixef(m))

    # Store results
    results$models <- models_list
    results$parameters <- params
    results$predictions <- preds_list
  } else {
    # Fit a single model to all data
    model <- fit_brms_model(data)

    # Generate predictions for plotting
    age_seq <- sort(unique(data$age))
    new_data <- data.frame(age = age_seq)
  preds <- stats::predict(model, newdata = new_data, probs = c(0.025, 0.975))

    preds_df <- cbind(as.data.frame(preds), new_data)
    preds_df$Model <- "von Bertalanffy"

    # Store results
    results$models <- model
    results$parameters <- brms::fixef(model)
    results$predictions <- preds_df
  }

  class(results) <- c("vb_brms", "list")
  return(results)
}
