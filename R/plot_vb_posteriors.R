#' Plot Posterior Predictive Checks for von Bertalanffy Bayesian Models
#'
#' This function creates posterior predictive check plots for Bayesian von Bertalanffy growth models.
#'
#' @param model A model object returned by fit_vb_brms()
#' @param ndraws Number of posterior draws to include (default 50)
#' @param theme_fn Optional ggplot2 theme function to apply (default theme_minimal)
#'
#' @return A list of ggplot2 objects
#'
#' @examples
#' \dontrun{
#' # Simple example with simulated data
#' age <- 1:15
#' length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5)
#' fit <- fit_vb_brms(age = age, length = length)
#' plots <- plot_vb_posteriors(fit)
#' plots$posterior_predictive
#' }
#'
#' @export
plot_vb_posteriors <- function(model, ndraws = 50, theme_fn = ggplot2::theme_minimal()) {
  # Check if brms is installed
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required for this function. Please install it.",
      call. = FALSE
    )
  }

  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop("Package 'bayesplot' is required for this function. Please install it.",
      call. = FALSE
    )
  }

  if (!inherits(model, "vb_brms")) {
    stop("This function only works with brms models from fit_vb_brms()")
  }

  plots <- list()

  # Function to create posterior plots for a specific model
  create_posterior_plots <- function(brms_model, sex_label = NULL) {
    # Prepare plot title suffix
    title_suffix <- if (is.null(sex_label)) "" else paste(" -", sex_label)

    # Posterior predictive check
    pp_check <- bayesplot::pp_check(brms_model, ndraws = ndraws) +
      ggplot2::labs(title = paste("Posterior Predictive Check", title_suffix)) +
      theme_fn

    # Parameter posterior densities
    post_samples <- brms::posterior_samples(brms_model)
    param_names <- c("b_Linf_Intercept", "b_k_Intercept", "b_t0_Intercept")
    param_labels <- c("Linf", "k", "t0")

    # Reshape for plotting
    post_df <- data.frame(
      parameter = rep(param_labels, each = nrow(post_samples)),
      value = c(
        post_samples$b_Linf_Intercept,
        post_samples$b_k_Intercept,
        post_samples$b_t0_Intercept
      )
    )

    # Create parameter density plots
    param_plot <- ggplot2::ggplot(post_df, ggplot2::aes(x = value)) +
      ggplot2::geom_density(fill = "blue", alpha = 0.5) +
      ggplot2::facet_wrap(~parameter, scales = "free") +
      ggplot2::labs(
        title = paste("Parameter Posterior Distributions", title_suffix),
        x = "Value", y = "Density"
      ) +
      theme_fn

    # Parameter pair plots
    pairs_plot <- NULL
    if (requireNamespace("GGally", quietly = TRUE)) {
      pairs_df <- data.frame(
        Linf = post_samples$b_Linf_Intercept,
        k = post_samples$b_k_Intercept,
        t0 = post_samples$b_t0_Intercept
      )

      pairs_plot <- GGally::ggpairs(pairs_df) +
        ggplot2::labs(title = paste("Parameter Correlations", title_suffix)) +
        theme_fn
    }

    # Conditional effects
    ce_plot <- ggplot2::ggplot() +
      ggplot2::geom_text(ggplot2::aes(
        x = 0.5, y = 0.5,
        label = "Use brms::conditional_effects() directly for more options"
      )) +
      ggplot2::theme_void() +
      ggplot2::labs(title = paste("Conditional Effects", title_suffix))

    if (requireNamespace("brms", quietly = TRUE)) {
      ce <- try(brms::conditional_effects(brms_model, spaghetti = TRUE, ndraws = 10), silent = TRUE)
      if (!inherits(ce, "try-error")) {
        ce_plot <- brms::plot(ce)[[1]] +
          ggplot2::labs(title = paste("Conditional Effects", title_suffix)) +
          theme_fn
      }
    }

    return(list(
      posterior_predictive = pp_check,
      parameter_densities = param_plot,
      parameter_pairs = pairs_plot,
      conditional_effects = ce_plot
    ))
  }

  # Check if we have multiple models by sex
  if (is.list(model$models) && !inherits(model$models, "brmsfit")) {
    # Create plots for each sex
    for (s in names(model$models)) {
      plots[[s]] <- create_posterior_plots(model$models[[s]], s)
    }
  } else {
    # Single model
    plots <- create_posterior_plots(model$models)
  }

  return(plots)
}
