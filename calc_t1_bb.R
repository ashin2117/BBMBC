#' Calculate Type I Error for Beta-Binomial Models
#'
#' Simulates immunology data under the null hypothesis of equal means
#' (no biological difference in cell frequency) using parallel processing.
#' Handles balanced/unbalanced sample sizes and dynamic overdispersion.
#'
#' @param n_0 Integer. The number of biological replicates in the control arm (default: 20).
#' @param n_1 Integer. The number of biological replicates in the treatment arm (default: 20).
#' @param p_base Numeric. Baseline proportion of the target cell population for BOTH groups (default: 0.005).
#' @param rho_0 Numeric. Intra-class correlation (overdispersion) for the control group (default: 0.006).
#' @param rho_1 Numeric. Intra-class correlation (overdispersion) for the treatment group (default: 0.006).
#' @param alpha Numeric. Significance threshold for the likelihood ratio test (default: 0.05).
#' @param min_cells Integer. Minimum number of total cells analyzed per sample (default: 1e5).
#' @param max_cells Integer. Maximum number of total cells analyzed per sample (default: 1e6).
#' @param seed Integer. Random seed for reproducibility (default: 2026).
#' @param n_sims Integer. The number of simulations to run (default: 1000).
#' @param n_cores Integer. Number of cores to use for parallel processing. Defaults to 2 (CRAN safe default).
#'
#' @return A list of class `t1_error_result` containing the estimated Type I error rate,
#'   convergence metrics, and the input parameters.
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
#' # Run small simulation example
#' res <- calc_t1_bb(
#'   n_0 = 10, n_1 = 10,
#'   rho_0 = 0.006, rho_1 = 0.006,
#'   n_sims = 10,
#'   n_cores = 1
#' )
#' print(res$type_1_error_rate)
#' }
calc_t1_bb <- function(n_0 = 20,
                       n_1 = 20,
                       p_base = 0.005,
                       rho_0 = 0.006,
                       rho_1 = 0.006,
                       alpha = 0.05,
                       min_cells = 1e5,
                       max_cells = 1e6,
                       seed = 2026,
                       n_sims = 1000,
                       n_cores = min(2, max(1, parallel::detectCores() - 1))) {

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)
  parallel::clusterSetRNGStream(cl, iseed = seed)

  shape1_0 <- p_base * (1 - rho_0) / rho_0
  shape2_0 <- (1 - p_base) * (1 - rho_0) / rho_0

  shape1_1 <- p_base * (1 - rho_1) / rho_1
  shape2_1 <- (1 - p_base) * (1 - rho_1) / rho_1

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
  t1_error <- if (length(valid_pvals) > 0) mean(valid_pvals <= alpha) else NA_real_

  structure(
    list(
      type_1_error_rate = t1_error,
      convergences = length(valid_pvals),
      failed_fits = n_sims - length(valid_pvals),
      parameters = list(
        n_0 = n_0,
        n_1 = n_1,
        p_base = p_base,
        rho_0 = rho_0,
        rho_1 = rho_1,
        alpha = alpha,
        model_dispformula = if (equal_dispersion) "~ 1" else "~ group",
        n_cores = n_cores
      )
    ),
    class = "t1_error_result"
  )
}
