#' Calculate Statistical Power for Beta-Binomial Models (Parallelized)
#'
#' Simulates overdispersed cell frequency data under the alternative hypothesis
#' to estimate the statistical power of a Beta-Binomial generalized linear model.
#'
#' @param n_0 Integer. Biological replicates in control arm (default: 20).
#' @param n_1 Integer. Biological replicates in treatment arm (default: 20).
#' @param p_base Numeric. Baseline proportion in control group (default: 0.005).
#' @param fold_change Numeric. Ratio of treatment mean to control mean (default: 2.0).
#' @param rho_0 Numeric. Overdispersion for control group (default: 0.006).
#' @param rho_1 Numeric. Overdispersion for treatment group (default: 0.006).
#' @param alpha Numeric. Significance threshold for likelihood ratio test (default: 0.05).
#' @param min_cells Integer. Minimum total cells per sample (default: 1e5).
#' @param max_cells Integer. Maximum total cells per sample (default: 1e6).
#' @param seed Integer. Random seed for reproducibility (default: 2025).
#' @param n_sims Integer. Number of simulations to run (default: 1000).
#' @param n_cores Integer. Number of cores for parallel processing. Defaults to 2 (CRAN safe default).
#'
#' @return A list of class `power_result` containing the estimated statistical power,
#'   convergence metrics, and the input parameters used.
#'
#' @importFrom glmmTMB glmmTMB betabinomial
#' @importFrom stats rbeta rbinom anova runif
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster clusterSetRNGStream detectCores
#' @export
#'
#' @examples
#' \donttest{
#' # Evaluate statistical power across a range of sample sizes
#' sample_sizes <- c(10, 15, 20, 25)
#' power_results <- numeric(length(sample_sizes))
#' names(power_results) <- paste0("N=", sample_sizes, "/arm")
#'
#' for (i in seq_along(sample_sizes)) {
#'   n <- sample_sizes[i]
#'
#'   # Run simulation for the current sample size
#'   # Note: n_sims = 100 is used here for demonstration speed.
#'   # For formal analysis, use n_sims = 1000 or higher.
#'   pwr_res <- calc_power_bb(
#'     n_0 = n, n_1 = n,
#'     p_base = 0.005,
#'     fold_change = 2.0,
#'     rho_0 = 0.006, rho_1 = 0.006,
#'     n_sims = 100,
#'     n_cores = 1
#'   )
#'
#'   power_results[i] <- pwr_res$power
#' }
#'
#' # Print the resulting power curve estimates
#' print(power_results)
#' }
calc_power_bb <- function(n_0 = 20,
                          n_1 = 20,
                          p_base = 0.005,
                          fold_change = 2.0,
                          rho_0 = 0.006,
                          rho_1 = 0.006,
                          alpha = 0.05,
                          min_cells = 1e5,
                          max_cells = 1e6,
                          seed = 2025,
                          n_sims = 1000,
                          n_cores = min(2, max(1, parallel::detectCores() - 1))) {

  p_treat <- p_base * fold_change
  if (p_treat >= 1) {
    stop("p_base * fold_change must be less than 1.")
  }

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  parallel::clusterSetRNGStream(cl, iseed = seed)

  shape1_0 <- p_base * (1 - rho_0) / rho_0
  shape2_0 <- (1 - p_base) * (1 - rho_0) / rho_0

  shape1_1 <- p_treat * (1 - rho_1) / rho_1
  shape2_1 <- (1 - p_treat) * (1 - rho_1) / rho_1

  group <- factor(c(rep(0, n_0), rep(1, n_1)))
  total_samples <- n_0 + n_1
  equal_dispersion <- isTRUE(all.equal(rho_0, rho_1))

  p_values <- foreach::foreach(
    i = seq_len(n_sims),
    .combine = "c",
    .packages = c("glmmTMB", "stats")
  ) %dopar% {
    n_cells <- as.integer(stats::runif(total_samples, min = min_cells, max = max_cells))
    beta_props_0 <- stats::rbeta(n_0, shape1_0, shape2_0)
    beta_props_1 <- stats::rbeta(n_1, shape1_1, shape2_1)
    beta_props <- c(beta_props_0, beta_props_1)

    successes <- stats::rbinom(total_samples, size = n_cells, prob = beta_props)
    failures <- n_cells - successes
    df <- data.frame(group = group, successes = successes, failures = failures)

    if (equal_dispersion) {
      fit_null <- suppressWarnings(try(
        glmmTMB::glmmTMB(cbind(successes, failures) ~ 1, dispformula = ~ 1, family = glmmTMB::betabinomial(link = "logit"), data = df),
        silent = TRUE
      ))
      fit_full <- suppressWarnings(try(
        glmmTMB::glmmTMB(cbind(successes, failures) ~ group, dispformula = ~ 1, family = glmmTMB::betabinomial(link = "logit"), data = df),
        silent = TRUE
      ))
    } else {
      fit_null <- suppressWarnings(try(
        glmmTMB::glmmTMB(cbind(successes, failures) ~ 1, dispformula = ~ group, family = glmmTMB::betabinomial(link = "logit"), data = df),
        silent = TRUE
      ))
      fit_full <- suppressWarnings(try(
        glmmTMB::glmmTMB(cbind(successes, failures) ~ group, dispformula = ~ group, family = glmmTMB::betabinomial(link = "logit"), data = df),
        silent = TRUE
      ))
    }

    iter_pval <- NA_real_
    if (!inherits(fit_null, "try-error") && !inherits(fit_full, "try-error")) {
      anova_res <- try(stats::anova(fit_null, fit_full), silent = TRUE)
      if (!inherits(anova_res, "try-error")) {
        iter_pval <- anova_res$"Pr(>Chisq)"[2]
      }
    }

    iter_pval
  }

  valid_pvals <- p_values[!is.na(p_values)]
  power_est <- if (length(valid_pvals) > 0) mean(valid_pvals <= alpha) else NA_real_

  structure(
    list(
      power = power_est,
      convergences = length(valid_pvals),
      failed_fits = n_sims - length(valid_pvals),
      parameters = list(
        n_0 = n_0,
        n_1 = n_1,
        p_base = p_base,
        fold_change = fold_change,
        p_treat = p_treat,
        rho_0 = rho_0,
        rho_1 = rho_1,
        alpha = alpha,
        model_dispformula = if (equal_dispersion) "~ 1" else "~ group",
        n_cores = n_cores
      )
    ),
    class = "power_result"
  )
}
