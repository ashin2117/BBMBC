#' Beta-Binomial Likelihood Ratio Test for Group Differences
#'
#' Performs a Likelihood Ratio Test (LRT) using Beta-Binomial models to evaluate
#' whether the underlying proportion of events differs significantly between two groups.
#'
#' @param data A data frame containing the experimental data.
#' @param success_col Character string specifying the column name for successes.
#' @param total_col Character string specifying the column name for total trials.
#' @param group_col Character string specifying the column name with group labels (must have 2 levels).
#' @param unequal_dispersion Logical. If TRUE, models group-specific overdispersion. Default: FALSE.
#'
#' @return A list of class `bb_test_result` containing p-value, chi-square statistic,
#'   degrees of freedom, log-odds ratio, and its standard error.
#'
#' @importFrom glmmTMB glmmTMB betabinomial fixef
#' @importFrom stats anova vcov
#' @export
#'
#' @examples
#' set.seed(2025)
#' trial_data <- data.frame(
#'   arm = rep(c("Control", "Vaccine"), each = 15),
#'   total_bcells = as.integer(stats::runif(30, 100000, 500000))
#' )
#' trial_data$target_cells <- ifelse(
#'   trial_data$arm == "Control",
#'   stats::rbinom(15, trial_data$total_bcells, stats::rbeta(15, 5, 995)),
#'   stats::rbinom(15, trial_data$total_bcells, stats::rbeta(15, 10, 990))
#' )
#'
#' result <- test_diff_bb(
#'   data = trial_data,
#'   success_col = "target_cells",
#'   total_col = "total_bcells",
#'   group_col = "arm",
#'   unequal_dispersion = FALSE
#' )
#' print(result$p_value)
test_diff_bb <- function(data, success_col, total_col, group_col, unequal_dispersion = FALSE) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }

  missing_cols <- setdiff(c(success_col, total_col, group_col), names(data))
  if (length(missing_cols) > 0) {
    stop("The following columns were not found in 'data': ", paste(missing_cols, collapse = ", "))
  }

  successes <- data[[success_col]]
  totals <- data[[total_col]]
  group <- data[[group_col]]

  if (any(successes > totals, na.rm = TRUE)) {
    stop("Success counts cannot exceed total counts.")
  }

  group_factor <- as.factor(group)
  if (length(levels(group_factor)) != 2) {
    stop("The 'group_col' variable must have exactly two distinct levels.")
  }

  df_model <- data.frame(
    group = group_factor,
    successes = successes,
    failures = totals - successes
  )

  disp_form <- if (unequal_dispersion) ~ group else ~ 1

  fit_null <- suppressWarnings(try(
    glmmTMB::glmmTMB(cbind(successes, failures) ~ 1,
                     dispformula = disp_form,
                     family = glmmTMB::betabinomial(link = "logit"),
                     data = df_model),
    silent = TRUE
  ))

  fit_full <- suppressWarnings(try(
    glmmTMB::glmmTMB(cbind(successes, failures) ~ group,
                     dispformula = disp_form,
                     family = glmmTMB::betabinomial(link = "logit"),
                     data = df_model),
    silent = TRUE
  ))

  res <- list(
    p_value = NA_real_,
    chisq_stat = NA_real_,
    df = NA_integer_,
    log_odds_ratio = NA_real_,
    se = NA_real_
  )

  if (!inherits(fit_null, "try-error") && !inherits(fit_full, "try-error")) {
    anova_res <- try(stats::anova(fit_null, fit_full), silent = TRUE)

    if (!inherits(anova_res, "try-error")) {
      res$p_value <- anova_res$"Pr(>Chisq)"[2]
      res$chisq_stat <- anova_res$Chisq[2]
      res$df <- anova_res$"Chi Df"[2]
    }

    coef_name <- paste0("group", levels(group_factor)[2])
    cond_fixef <- glmmTMB::fixef(fit_full)$cond

    if (coef_name %in% names(cond_fixef)) {
      res$log_odds_ratio <- cond_fixef[coef_name]
      res$se <- sqrt(diag(stats::vcov(fit_full)$cond)[coef_name])
    }
  }

  structure(res, class = "bb_test_result")
}
