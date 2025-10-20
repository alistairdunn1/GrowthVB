#' Comprehensive Bayesian Diagnostics for von Bertalanffy Models
#'
#' This function creates comprehensive MCMC and posterior predictive diagnostic plots
#' for Bayesian von Bertalanffy growth models, including convergence diagnostics,
#' enhanced posterior predictive checks with violin plots showing full posterior
#' distributions of residuals, and model-specific growth diagnostics.
#'
#' @param model A model object returned by fit_vb_brms()
#' @param ndraws Number of posterior draws to include in plots (default 50)
#' @param diagnostic_types Character vector specifying which diagnostics to include.
#'   Options: "convergence", "pp_checks", "residuals", "parameters", "growth_specific", "all"
#' @param pars Optional character vector of parameter names to focus on
#'
#' @return A list of ggplot2 objects organised by diagnostic type
#'
#' @examples
#' \dontrun{
#' # Simple example with simulated data
#' age <- 1:15
#' length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5)
#' fit <- fit_vb_brms(age = age, length = length)
#'
#' # Full diagnostics
#' diagnostics <- plot_vb_bayes_diagnostics(fit)
#'
#' # Specific diagnostic types
#' convergence <- plot_vb_bayes_diagnostics(fit, diagnostic_types = "convergence")
#' pp_checks <- plot_vb_bayes_diagnostics(fit, diagnostic_types = "pp_checks")
#' }
#'
#' @export
plot_vb_bayes_diagnostics <- function(model,
                                      ndraws = 100,
                                      diagnostic_types = "all",
                                      pars = NULL) {
  # Check required packages
  required_packages <- c("brms", "bayesplot", "posterior", "tidyr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required for this function. Please install it."),
        call. = FALSE
      )
    }
  }

  if (!inherits(model, "vb_brms")) {
    stop("This function only works with brms models from fit_vb_brms()")
  }

  # Set default diagnostic types
  if ("all" %in% diagnostic_types) {
    diagnostic_types <- c("convergence", "pp_checks", "residuals", "parameters", "growth_specific")
  }

  diagnostics <- list()

  # Function to create diagnostics for a specific model
  create_diagnostics <- function(brms_model, sex_label = NULL) {
    title_suffix <- if (is.null(sex_label)) "" else paste0("(", sex_label, ")")

    model_diagnostics <- list()

    # 1. CONVERGENCE DIAGNOSTICS
    if ("convergence" %in% diagnostic_types) {
      # Extract posterior samples for convergence checks
      posterior_samples <- brms::as_draws_array(brms_model)

      # Default parameters if not specified
      if (is.null(pars)) {
        pars <- c("b_Linf_Intercept", "b_k_Intercept", "b_t0_Intercept", "b_CV_Intercept")
        # Filter to available parameters
        available_pars <- dimnames(posterior_samples)$variable
        pars <- pars[pars %in% available_pars]
      }

      # Trace plots
      trace_plot <- bayesplot::mcmc_trace(posterior_samples, pars = pars) +
        ggplot2::labs(title = paste("MCMC Trace Plots", title_suffix))

      # Rank plots for convergence
      rank_plot <- bayesplot::mcmc_rank_overlay(posterior_samples, pars = pars) +
        ggplot2::labs(title = paste("MCMC Rank Plots", title_suffix))

      # R-hat diagnostic
      rhats <- brms::rhat(brms_model)
      rhat_plot <- bayesplot::mcmc_rhat(rhats) +
        ggplot2::labs(title = paste("R-hat Diagnostics", title_suffix))

      # Effective sample size
      neff_ratios <- brms::neff_ratio(brms_model)
      neff_plot <- bayesplot::mcmc_neff(neff_ratios) +
        ggplot2::labs(title = paste("Effective Sample Size Ratios", title_suffix))

      # Autocorrelation
      acf_plot <- bayesplot::mcmc_acf(posterior_samples, pars = pars) +
        ggplot2::labs(title = paste("Autocorrelation Functions", title_suffix))

      model_diagnostics$convergence <- list(
        trace = trace_plot,
        rank = rank_plot,
        rhat = rhat_plot,
        neff = neff_plot,
        autocorr = acf_plot
      )
    }

    # 2. ENHANCED POSTERIOR PREDICTIVE CHECKS
    if ("pp_checks" %in% diagnostic_types) {
      # Basic posterior predictive check
      pp_basic <- bayesplot::pp_check(brms_model, ndraws = ndraws) +
        ggplot2::labs(title = paste("Posterior Predictive Check - Density", title_suffix))

      # PP check for mean
      pp_mean <- bayesplot::pp_check(brms_model, type = "stat", stat = "mean", ndraws = ndraws) +
        ggplot2::labs(title = paste("PP Check - Mean", title_suffix))

      # PP check for standard deviation
      pp_sd <- bayesplot::pp_check(brms_model, type = "stat", stat = "sd", ndraws = ndraws) +
        ggplot2::labs(title = paste("PP Check - Standard Deviation", title_suffix))

      # PP check for quantiles
      pp_quantiles <- bayesplot::pp_check(brms_model, type = "stat_2d", stat = c("mean", "sd"), ndraws = ndraws) +
        ggplot2::labs(title = paste("PP Check - Mean vs SD", title_suffix))

      # Rootogram
      pp_rootogram <- try(
        {
          bayesplot::pp_check(brms_model, type = "rootogram", ndraws = ndraws) +
            ggplot2::labs(title = paste("PP Check - Rootogram", title_suffix))
        },
        silent = TRUE
      )
      if (inherits(pp_rootogram, "try-error")) {
        pp_rootogram <- ggplot2::ggplot() +
          ggplot2::geom_text(ggplot2::aes(x = 0.5, y = 0.5, label = "Rootogram not available")) +
          ggplot2::theme_void()
      }

      model_diagnostics$pp_checks <- list(
        density = pp_basic,
        mean = pp_mean,
        sd = pp_sd,
        quantiles = pp_quantiles,
        rootogram = pp_rootogram
      )
    }

    # 3. RESIDUAL DIAGNOSTICS
    if ("residuals" %in% diagnostic_types) {
      # Get residuals - these are S3 methods that should work with brms models
      residuals_result <- try(residuals(brms_model, summary = FALSE), silent = TRUE)
      fitted_result <- try(fitted(brms_model, summary = FALSE), silent = TRUE)

      if (inherits(residuals_result, "try-error") || inherits(fitted_result, "try-error")) {
        # Skip residual diagnostics if methods fail
        model_diagnostics$residuals <- list(
          error_message = ggplot2::ggplot() +
            ggplot2::geom_text(ggplot2::aes(
              x = 0.5, y = 0.5,
              label = "Residual diagnostics not available for this model"
            )) +
            ggplot2::theme_void()
        )
      } else {
        residuals <- residuals_result
        fitted_values <- fitted_result

        # Prepare data for violin plots
        # Create bins for fitted values to group residuals
        fitted_mean <- apply(fitted_values, 2, mean)
        n_bins <- min(10, max(5, length(fitted_mean) %/% 10)) # Adaptive binning
        fitted_bins <- cut(fitted_mean, breaks = n_bins, labels = FALSE)

        # Create long-form data for violin plots
        residual_long <- data.frame()
        for (i in 1:ncol(residuals)) {
          temp_df <- data.frame(
            observation = i,
            fitted_mean = fitted_mean[i],
            fitted_bin = fitted_bins[i],
            residual = residuals[, i]
          )
          residual_long <- rbind(residual_long, temp_df)
        }

        # Calculate bin centres for x-axis positioning
        bin_centers <- aggregate(fitted_mean ~ fitted_bin, residual_long, mean, na.rm = TRUE)
        residual_long$fitted_center <- bin_centers$fitted_mean[match(residual_long$fitted_bin, bin_centers$fitted_bin)]

        # Violin plot of residual distributions
        residual_violin <- ggplot2::ggplot(residual_long, ggplot2::aes(x = factor(fitted_bin), y = residual)) +
          ggplot2::geom_violin(fill = "royalblue", alpha = 0.6, colour = "darkblue") +
          ggplot2::geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.size = 0.5) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "red", size = 1) +
          ggplot2::scale_x_discrete(labels = round(bin_centers$fitted_mean, 1)) +
          ggplot2::labs(
            title = paste("Posterior Residual Distributions vs Fitted", title_suffix),
            x = "Fitted Values (binned)", y = "Residuals"
          )

        # Alternative: point + ribbon plot showing quantiles
        residual_summary <- aggregate(residual ~ fitted_bin, residual_long, function(x) {
          quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
        })
        residual_quantiles <- data.frame(
          fitted_center = bin_centers$fitted_mean,
          q05 = residual_summary$residual[, 1],
          q25 = residual_summary$residual[, 2],
          median = residual_summary$residual[, 3],
          q75 = residual_summary$residual[, 4],
          q95 = residual_summary$residual[, 5]
        )

        residual_quantile_plot <- ggplot2::ggplot(residual_quantiles, ggplot2::aes(x = fitted_center)) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = q05, ymax = q95), fill = "royalblue", alpha = 0.3) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75), fill = "royalblue", alpha = 0.5) +
          ggplot2::geom_line(ggplot2::aes(y = median), colour = "darkblue", size = 1) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
          ggplot2::labs(
            title = paste("Residual Quantiles vs Fitted", title_suffix),
            x = "Fitted Values", y = "Residuals",
            subtitle = "Dark band: 50% interval, Light band: 90% interval, Line: median"
          )

        model_diagnostics$residuals <- list(
          violin_plot = residual_violin,
          quantile_bands = residual_quantile_plot
        )
      }
    }

    # 4. PARAMETER DIAGNOSTICS
    if ("parameters" %in% diagnostic_types) {
      # Parameter densities with prior overlays (if available)
      posterior_samples <- brms::posterior_samples(brms_model)
      param_cols <- grep("^b_.*_Intercept$", names(posterior_samples), value = TRUE)

      if (length(param_cols) > 0) {
        # Reshape for plotting
        param_data <- posterior_samples[param_cols]
        names(param_data) <- gsub("b_(.*)_Intercept", "\\1", names(param_data))

        param_long <- tidyr::pivot_longer(param_data, everything(), names_to = "parameter", values_to = "value")

        param_density <- ggplot2::ggplot(param_long, ggplot2::aes(x = value)) +
          ggplot2::geom_density(fill = "royalblue", alpha = 0.6) +
          ggplot2::facet_wrap(~parameter, scales = "free") +
          ggplot2::labs(
            title = paste("Parameter Posterior Densities", title_suffix),
            x = "Value", y = "Density"
          )

        # Parameter intervals
        param_intervals <- bayesplot::mcmc_areas(
          as.matrix(param_data),
          prob = 0.8, # 80% intervals
          prob_outer = 0.95 # 95% intervals
        ) +
          ggplot2::labs(title = paste("Parameter Credible Intervals", title_suffix))

        # Parameter scatter plots (correlations) - key for von Bertalanffy parameters
        param_scatter <- NULL
        if (ncol(param_data) >= 2) {
          # Create pairwise scatter plots to show parameter correlations
          param_pairs <- bayesplot::mcmc_pairs(
            as.matrix(param_data),
            pars = names(param_data)[1:min(4, ncol(param_data))], # Limit to first 4 parameters
            off_diag_args = list(size = 0.5, alpha = 0.6),
            diag_fun = "dens"
          ) +
            ggplot2::labs(title = paste("Parameter Correlations", title_suffix))
          param_scatter <- param_pairs
        }

        model_diagnostics$parameters <- list(
          densities = param_density,
          intervals = param_intervals,
          correlations = param_scatter
        )
      }
    }

    # 5. GROWTH-SPECIFIC DIAGNOSTICS
    if ("growth_specific" %in% diagnostic_types) {
      # Get original data
      original_data <- brms_model$data

      if ("age" %in% names(original_data) && "length" %in% names(original_data)) {
        # Age-specific residuals
        residuals <- try(residuals(brms_model, summary = TRUE), silent = TRUE)
        if (inherits(residuals, "try-error")) {
          # Skip if residuals method fails
          model_diagnostics$growth_specific <- list(
            error_message = ggplot2::ggplot() +
              ggplot2::geom_text(ggplot2::aes(
                x = 0.5, y = 0.5,
                label = "Growth-specific diagnostics not available"
              )) +
              ggplot2::theme_void()
          )
        } else {
          # Get full residuals matrix (not summary) for violin plots
          residuals_full <- try(residuals(brms_model, summary = FALSE), silent = TRUE)
          if (!inherits(residuals_full, "try-error")) {
            # Create bins for age to group residuals
            age_values <- original_data$age
            n_age_bins <- min(8, max(4, length(unique(age_values)))) # Adaptive binning
            age_bins <- cut(age_values, breaks = n_age_bins, labels = FALSE)

            # Create long-form data for violin plots by age
            age_residual_long <- data.frame()
            for (i in 1:ncol(residuals_full)) {
              temp_df <- data.frame(
                observation = i,
                age = age_values[i],
                age_bin = age_bins[i],
                residual = residuals_full[, i]
              )
              age_residual_long <- rbind(age_residual_long, temp_df)
            }

            # Calculate bin centres for x-axis positioning
            age_bin_centers <- aggregate(age ~ age_bin, age_residual_long, mean, na.rm = TRUE)
            age_residual_long$age_center <- age_bin_centers$age[match(age_residual_long$age_bin, age_bin_centers$age_bin)]

            # Violin plot of residual distributions by age
            age_residual_violin <- ggplot2::ggplot(age_residual_long, ggplot2::aes(x = factor(age_bin), y = residual)) +
              ggplot2::geom_violin(fill = "royalblue", alpha = 0.6, colour = "darkblue") +
              ggplot2::geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.size = 0.5) +
              ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "red", size = 1) +
              ggplot2::scale_x_discrete(labels = round(age_bin_centers$age, 1)) +
              ggplot2::labs(
                title = paste("Posterior Residual Distributions by Age", title_suffix),
                x = "Age (binned)", y = "Residuals"
              )

            age_residual_plot <- age_residual_violin
          } else {
            # Fallback to summary residuals if full residuals fail
            residuals_summary <- try(residuals(brms_model, summary = TRUE), silent = TRUE)
            if (!inherits(residuals_summary, "try-error")) {
              age_residual_data <- data.frame(
                age = original_data$age,
                residual = residuals_summary[, "Estimate"],
                residual_lower = residuals_summary[, "Q2.5"],
                residual_upper = residuals_summary[, "Q97.5"]
              )

              age_residual_plot <- ggplot2::ggplot(age_residual_data, ggplot2::aes(x = age, y = residual)) +
                ggplot2::geom_point(alpha = 0.6, colour = "royalblue") +
                ggplot2::geom_errorbar(ggplot2::aes(ymin = residual_lower, ymax = residual_upper),
                  alpha = 0.3, colour = "royalblue"
                ) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
                ggplot2::geom_smooth(se = FALSE, colour = "black", method = "loess") +
                ggplot2::labs(
                  title = paste("Residuals by Age (Summary)", title_suffix),
                  x = "Age", y = "Residuals"
                )
            } else {
              age_residual_plot <- ggplot2::ggplot() +
                ggplot2::geom_text(ggplot2::aes(x = 0.5, y = 0.5, label = "Age residuals not available")) +
                ggplot2::theme_void()
            }
          }

          # Growth curve with prediction intervals
          age_seq <- seq(min(original_data$age), max(original_data$age), length.out = 50)
          new_data <- data.frame(age = age_seq)
          predictions <- try(posterior_predict(brms_model, newdata = new_data, summary = TRUE), silent = TRUE)
          if (inherits(predictions, "try-error")) {
            predictions <- fitted(brms_model, newdata = new_data, probs = c(0.025, 0.975))
          }

          pred_data <- data.frame(
            age = age_seq,
            fit = predictions[, "Estimate"],
            lower = predictions[, "Q2.5"],
            upper = predictions[, "Q97.5"]
          )

          growth_curve_plot <- ggplot2::ggplot() +
            ggplot2::geom_ribbon(
              data = pred_data,
              ggplot2::aes(x = age, ymin = lower, ymax = upper),
              fill = "royalblue", alpha = 0.3
            ) +
            ggplot2::geom_line(
              data = pred_data,
              ggplot2::aes(x = age, y = fit),
              colour = "royalblue", linewidth = 1
            ) +
            ggplot2::geom_point(
              data = original_data,
              ggplot2::aes(x = age, y = length),
              alpha = 0.6
            ) +
            ggplot2::labs(
              title = paste("Growth Curve with Prediction Intervals", title_suffix),
              x = "Age", y = "Length"
            )

          # Enhanced prediction uncertainty visualization with ribbons
          pred_uncertainty_plot <- NULL
          age_range <- range(original_data$age)
          detailed_ages <- seq(age_range[1], age_range[2], length.out = 50)
          detailed_pred_data <- data.frame(age = detailed_ages)

          # Get full posterior predictions for ribbon plots
          pred_samples <- try(brms::posterior_predict(brms_model, newdata = detailed_pred_data), silent = TRUE)
          if (!inherits(pred_samples, "try-error")) {
            # Calculate prediction quantiles for multiple ribbon layers
            pred_quantiles <- data.frame(
              age = detailed_ages,
              median = apply(pred_samples, 2, quantile, 0.5),
              q10 = apply(pred_samples, 2, quantile, 0.1),
              q25 = apply(pred_samples, 2, quantile, 0.25),
              q75 = apply(pred_samples, 2, quantile, 0.75),
              q90 = apply(pred_samples, 2, quantile, 0.9),
              q05 = apply(pred_samples, 2, quantile, 0.05),
              q95 = apply(pred_samples, 2, quantile, 0.95)
            )

            pred_uncertainty_plot <- ggplot2::ggplot(pred_quantiles, ggplot2::aes(x = age)) +
              # Outer ribbon (90% prediction interval)
              ggplot2::geom_ribbon(ggplot2::aes(ymin = q05, ymax = q95),
                fill = "lightblue", alpha = 0.3
              ) +
              # Middle ribbon (80% prediction interval)
              ggplot2::geom_ribbon(ggplot2::aes(ymin = q10, ymax = q90),
                fill = "royalblue", alpha = 0.4
              ) +
              # Inner ribbon (50% prediction interval)
              ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75),
                fill = "darkblue", alpha = 0.5
              ) +
              # Median line
              ggplot2::geom_line(ggplot2::aes(y = median), colour = "darkblue", linewidth = 1.2) +
              # Observed data points
              ggplot2::geom_point(
                data = original_data,
                ggplot2::aes(x = age, y = length),
                alpha = 0.7, colour = "red", size = 1.5
              ) +
              ggplot2::labs(
                title = paste("Enhanced Prediction Uncertainty", title_suffix),
                x = "Age", y = "Predicted Length",
                subtitle = "Ribbons: 50% (dark), 80% (medium), 90% (light) prediction intervals"
              )
          }

          model_diagnostics$growth_specific <- list(
            age_residuals = age_residual_plot,
            growth_curve = growth_curve_plot,
            prediction_uncertainty = pred_uncertainty_plot
          )
        }
      }
    }

    return(model_diagnostics)
  }

  # Create diagnostics based on model structure
  if (is.list(model$models) && !inherits(model$models, "brmsfit")) {
    # Multiple models by sex
    for (s in names(model$models)) {
      diagnostics[[s]] <- create_diagnostics(model$models[[s]], s)
    }
    message("Bayesian diagnostics organised by sex: ", paste(names(diagnostics), collapse = ", "))
  } else {
    # Single model
    diagnostics <- create_diagnostics(model$models)
    message("Bayesian diagnostics for combined model")
  }

  return(diagnostics)
}
