#' Growth-Specific Posterior Predictive Checks
#'
#' This function creates posterior predictive checks specifically designed for
#' von Bertalanffy growth models, focusing on growth-relevant patterns and diagnostics.
#'
#' @param model A model object returned by fit_vb_brms()
#' @param ndraws Number of posterior draws to include (default 100)
#' @param check_types Character vector specifying which checks to perform.
#'   Options: "growth_patterns", "age_effects", "length_distribution", "residual_patterns", "all"
#'
#' @return A list of ggplot2 objects with growth-specific posterior predictive checks
#'
#' @examples
#' \dontrun{
#' # Growth model with posterior predictive checks
#' age <- 1:15
#' length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5)
#' fit <- fit_vb_brms(age = age, length = length)
#'
#' # Growth-specific checks
#' growth_checks <- plot_vb_growth_pp_checks(fit)
#' }
#'
#' @export
plot_vb_growth_pp_checks <- function(model,
                                     ndraws = 100,
                                     check_types = "all") {
  # Check required packages
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required for this function. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop("Package 'bayesplot' is required for this function. Please install it.", call. = FALSE)
  }

  if (!inherits(model, "vb_brms")) {
    stop("This function only works with brms models from fit_vb_brms()")
  }

  # Set default check types
  if ("all" %in% check_types) {
    check_types <- c("growth_patterns", "age_effects", "length_distribution", "residual_patterns")
  }

  pp_checks <- list()

  # Function to create growth-specific PP checks
  create_growth_pp_checks <- function(brms_model, sex_label = NULL) {
    title_suffix <- if (is.null(sex_label)) "" else paste(" -", sex_label)

    checks <- list()

    # Get observed data
    y_obs <- brms::get_y(brms_model)
    original_data <- brms_model$data

    # Generate posterior predictions
    y_rep <- try(posterior_predict(brms_model, ndraws = ndraws), silent = TRUE)
    if (inherits(y_rep, "try-error")) {
      stop("Could not generate posterior predictions from model")
    }

    # 1. GROWTH PATTERNS
    if ("growth_patterns" %in% check_types) {
      # Growth trajectory comparisons
      if ("age" %in% names(original_data)) {
        # Create age bins for trajectory comparison
        age_breaks <- quantile(original_data$age, probs = seq(0, 1, 0.25))
        age_groups <- cut(original_data$age, breaks = age_breaks, include.lowest = TRUE)

        # Mean length by age group - observed vs predicted
        obs_by_age <- tapply(y_obs, age_groups, mean, na.rm = TRUE)
        pred_by_age <- apply(y_rep, 1, function(x) tapply(x, age_groups, mean, na.rm = TRUE))

        age_comparison_data <- data.frame(
          age_group = rep(names(obs_by_age), ndraws + 1),
          mean_length = c(obs_by_age, as.vector(pred_by_age)),
          type = rep(c("Observed", rep("Predicted", ndraws)), each = length(obs_by_age))
        )

        growth_pattern_plot <- ggplot2::ggplot(age_comparison_data, ggplot2::aes(x = age_group, y = mean_length)) +
          ggplot2::geom_boxplot(
            data = subset(age_comparison_data, type == "Predicted"),
            ggplot2::aes(group = age_group), alpha = 0.3, fill = "lightblue"
          ) +
          ggplot2::geom_point(
            data = subset(age_comparison_data, type == "Observed"),
            colour = "red", size = 3
          ) +
          ggplot2::labs(
            title = paste("Growth Patterns: Mean Length by Age Group", title_suffix),
            x = "Age Group", y = "Mean Length",
            caption = "Red points: observed, Boxes: posterior predictive distribution"
          ) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

        checks$growth_patterns <- growth_pattern_plot
      }
    }

    # 2. AGE EFFECTS
    if ("age_effects" %in% check_types && "age" %in% names(original_data)) {
      # Length-at-age variability
      age_unique <- sort(unique(original_data$age))

      if (length(age_unique) > 3) {
        # Calculate CV by age for observed and predicted
        cv_obs <- tapply(seq_along(y_obs), original_data$age, function(idx) {
          vals <- y_obs[idx]
          if (length(vals) > 1) sd(vals) / mean(vals) else NA
        })

        cv_pred <- apply(y_rep, 1, function(x) {
          tapply(seq_along(x), original_data$age, function(idx) {
            vals <- x[idx]
            if (length(vals) > 1) sd(vals) / mean(vals) else NA
          })
        })

        cv_data <- data.frame(
          age = rep(as.numeric(names(cv_obs)), ndraws + 1),
          cv = c(cv_obs, as.vector(unlist(cv_pred))),
          type = rep(c("Observed", rep("Predicted", ndraws)), each = length(cv_obs))
        )
        cv_data <- cv_data[!is.na(cv_data$cv), ]

        age_cv_plot <- ggplot2::ggplot(cv_data, ggplot2::aes(x = age, y = cv)) +
          ggplot2::geom_line(
            data = subset(cv_data, type == "Predicted"),
            ggplot2::aes(group = interaction(type, rep(1:ndraws, each = length(cv_obs)))),
            alpha = 0.1, colour = "lightblue"
          ) +
          ggplot2::geom_point(
            data = subset(cv_data, type == "Observed"),
            colour = "red", size = 2
          ) +
          ggplot2::geom_smooth(
            data = subset(cv_data, type == "Observed"),
            se = FALSE, colour = "red", method = "loess"
          ) +
          ggplot2::labs(
            title = paste("Coefficient of Variation by Age", title_suffix),
            x = "Age", y = "CV (SD/Mean)",
            caption = "Red: observed, Blue lines: posterior predictions"
          )

        checks$age_effects <- age_cv_plot
      }
    }

    # 3. LENGTH DISTRIBUTION
    if ("length_distribution" %in% check_types) {
      # Overall length distribution comparison
      length_dist_plot <- bayesplot::pp_check(brms_model, ndraws = ndraws, type = "dens_overlay") +
        ggplot2::labs(title = paste("Length Distribution PP Check", title_suffix))

      # Quantile-quantile plot
      qq_data <- data.frame(
        theoretical = sort(apply(y_rep, 2, median)),
        observed = sort(y_obs)
      )

      qq_plot <- ggplot2::ggplot(qq_data, ggplot2::aes(x = theoretical, y = observed)) +
        ggplot2::geom_point(alpha = 0.6, colour = "royalblue") +
        ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
        ggplot2::labs(
          title = paste("Q-Q Plot: Observed vs Predicted", title_suffix),
          x = "Predicted Quantiles", y = "Observed Quantiles"
        )

      checks$length_distribution <- list(
        density = length_dist_plot,
        qq_plot = qq_plot
      )
    }

    # 4. RESIDUAL PATTERNS
    if ("residual_patterns" %in% check_types) {
      # Prediction intervals vs observed
      fitted_summary <- try(fitted(brms_model, probs = c(0.025, 0.975)), silent = TRUE)
      if (inherits(fitted_summary, "try-error")) {
        # Skip this check if fitted method fails
        checks$residual_patterns <- list(
          error_message = ggplot2::ggplot() +
            ggplot2::geom_text(ggplot2::aes(
              x = 0.5, y = 0.5,
              label = "Residual pattern checks not available"
            )) +
            ggplot2::theme_void()
        )
      } else {
        # Calculate coverage
        coverage <- mean(y_obs >= fitted_summary[, "Q2.5"] & y_obs <= fitted_summary[, "Q97.5"])

        coverage_data <- data.frame(
          fitted = fitted_summary[, "Estimate"],
          observed = y_obs,
          lower = fitted_summary[, "Q2.5"],
          upper = fitted_summary[, "Q97.5"],
          in_interval = y_obs >= fitted_summary[, "Q2.5"] & y_obs <= fitted_summary[, "Q97.5"]
        )

        coverage_plot <- ggplot2::ggplot(coverage_data, ggplot2::aes(x = fitted, y = observed)) +
          ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
            alpha = 0.3, colour = "lightblue"
          ) +
          ggplot2::geom_point(ggplot2::aes(colour = in_interval), alpha = 0.6) +
          ggplot2::geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
          ggplot2::scale_colour_manual(
            values = c("FALSE" = "red", "TRUE" = "royalblue"),
            name = "In 95% PI"
          ) +
          ggplot2::labs(
            title = paste("Prediction Interval Coverage", title_suffix),
            subtitle = paste("Coverage:", round(coverage * 100, 1), "%"),
            x = "Fitted Values", y = "Observed Values"
          )

        # Residual autocorrelation (if ordered by age)
        if ("age" %in% names(original_data)) {
          residuals <- y_obs - fitted_summary[, "Estimate"]
          age_order <- order(original_data$age)

          residual_acf_data <- data.frame(
            lag1_residual = residuals[age_order[-length(age_order)]],
            residual = residuals[age_order[-1]]
          )

          acf_plot <- ggplot2::ggplot(residual_acf_data, ggplot2::aes(x = lag1_residual, y = residual)) +
            ggplot2::geom_point(alpha = 0.6, colour = "royalblue") +
            ggplot2::geom_smooth(se = TRUE, colour = "black", method = "lm") +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
            ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
            ggplot2::labs(
              title = paste("Residual Autocorrelation (Age-ordered)", title_suffix),
              x = "Residual at Age t", y = "Residual at Age t+1"
            )

          checks$residual_patterns <- list(
            coverage = coverage_plot,
            autocorrelation = acf_plot
          )
        } else {
          checks$residual_patterns <- list(coverage = coverage_plot)
        }
      }
    }

    return(checks)
  }

  # Create checks based on model structure
  if (is.list(model$models) && !inherits(model$models, "brmsfit")) {
    # Multiple models by sex
    for (s in names(model$models)) {
      pp_checks[[s]] <- create_growth_pp_checks(model$models[[s]], s)
    }
    message("Growth PP checks organised by sex: ", paste(names(pp_checks), collapse = ", "))
  } else {
    # Single model
    pp_checks <- create_growth_pp_checks(model$models)
    message("Growth PP checks for combined model")
  }

  return(pp_checks)
}
