#' Fit von Bertalanffy Growth Function using Bayesian Methods (brms)
#'
#' This function fits a von Bertalanffy growth function to age and length data using
#' Bayesian methods via the brms package. The CV of the length is modeled as a function
#' of the mean length.
#'
#' @param age A numeric vector of ages
#' @param length A numeric vector of lengths
#' @param sex An optional factor or character vector specifying the sex for each observation
#' @param priors Optional list of priors for the model parameters
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
#' @export
fit_vb_brms <- function(age, length, sex = NULL,
                        priors = NULL, chains = 4, iter = 4000, ...) {
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

  # Set default priors if not provided
  if (is.null(priors)) {
    priors <- c(
      brms::prior(normal(150, 100), nlpar = "Linf", lb = 0),
      brms::prior(normal(0.3, 100), nlpar = "k", lb = 0),
      brms::prior(normal(0, 1), nlpar = "t0"),
      brms::prior(normal(0, 100), nlpar = "tau", lb = 0)
    )
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
      preds <- brms::predict(model, newdata = new_data, probs = c(0.025, 0.975))

      preds_df <- data.frame(preds) %>%
        dplyr::bind_cols(new_data) %>%
        dplyr::mutate(Sex = s, Model = "von Bertalanffy")

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
    preds <- brms::predict(model, newdata = new_data, probs = c(0.025, 0.975))

    preds_df <- data.frame(preds) %>%
      dplyr::bind_cols(new_data) %>%
      dplyr::mutate(Model = "von Bertalanffy")

    # Store results
    results$models <- model
    results$parameters <- brms::fixef(model)
    results$predictions <- preds_df
  }

  class(results) <- c("vb_brms", "list")
  return(results)
}
