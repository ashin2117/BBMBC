#' Calculate Type I Error Rate for Wilcoxon Rank-Sum Test
#'
#' @description
#' Simulates overdispersed count data based on a Beta-Binomial distribution to evaluate 
#' the empirical Type I error rate of the non-parametric Wilcoxon rank-sum test. This function 
#' is utilized to benchmark rank-based testing against parametric frameworks when analyzing 
#' the relative abundance of rare memory B-cell responses in flow cytometry data.
#'
#' @param n_0 Integer. Sample size for the control group. Default is 20.
#' @param n_1 Integer. Sample size for the treatment group. Default is 20.
#' @param p_base Numeric. Expected baseline proportion of the target cell population. Default is 0.005.
#' @param rho_0 Numeric. Intra-class correlation (dispersion) for the control group. Default is 0.006.
#' @param rho_1 Numeric. Intra-class correlation (dispersion) for the treatment group. Default is 0.006.
#' @param alpha Numeric. Nominal significance level for hypothesis testing. Default is 0.05.
#' @param min_cells Numeric. Minimum bound for the uniform distribution of the parent cell pool. Default is 1e5.
#' @param max_cells Numeric. Maximum bound for the uniform distribution of the parent cell pool. Default is 1e6.
#' @param seed Integer. Random seed for reproducible simulations. Default is 2026.
#' @param n_sims Integer. Number of simulation iterations. Default is 1000.
#' @param n_cores Integer. Number of cores to use for parallel processing. Default is 1. 
#'
#' @return An S3 object of class \code{t1_error_result_wilcox} containing the empirical Type I error rate, 
#' the number of successful convergences, failed fits, and the input parameters.
#'
#' @importFrom stats runif rbeta rbinom wilcox.test
#' @importFrom parallel makeCluster stopCluster clusterSetRNGStream
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar% %do%
#' @export
#'
#' @examples
#' # Simulate Type I error for comparing memory B-cell proportions 
#' # with equal biological overdispersion between arms.
#' res <- calc_t1_wilcox(
#'   n_0 = 15,
#'   n_1 = 15,
#'   p_base = 0.01,
#'   rho_0 = 0.01,
#'   rho_1 = 0.01,
#'   n_sims = 100,
#'   n_cores = 1
#' )
#' print(res$type_1_error_rate)
calc_t1_wilcox <- function(n_0 = 20,
                           n_1 = 20,
                           p_base = 0.005,
                           rho_0 = 0.006,
                           rho_1 = 0.006,
                           alpha = 0.05,
                           min_cells = 1e5,
                           max_cells = 1e6,
                           seed = 2026,
                           n_sims = 1000,
                           n_cores = 1) {
  
  # Input Validation
  stopifnot(
    "n_0 and n_1 must be strictly positive integers" = n_0 > 0 && n_1 > 0,
    "p_base must be between 0 and 1" = p_base > 0 && p_base < 1,
    "rho_0 and rho_1 must be between 0 and 1" = rho_0 > 0 && rho_0 < 1 && rho_1 > 0 && rho_1 < 1,
    "alpha must be between 0 and 1" = alpha > 0 && alpha < 1,
    "max_cells must be greater than or equal to min_cells" = max_cells >= min_cells,
    "n_cores must be a positive integer" = n_cores >= 1
  )
  
  shape1_0 <- p_base * (1 - rho_0) / rho_0
  shape2_0 <- (1 - p_base) * (1 - rho_0) / rho_0
  
  shape1_1 <- p_base * (1 - rho_1) / rho_1
  shape2_1 <- (1 - p_base) * (1 - rho_1) / rho_1
  
  group <- factor(c(rep(0, n_0), rep(1, n_1)))
  total_samples <- n_0 + n_1
  
  # Safe Parallel/Sequential Setup
  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
    parallel::clusterSetRNGStream(cl, iseed = seed)
    `%op%` <- foreach::`%dopar%`
  } else {
    # Isolate global environment modifications
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
    `%op%` <- foreach::`%do%`
  }
  
  i <- NULL 
  
  p_values <- foreach::foreach(
    i = seq_len(n_sims),
    .combine = "c",
    .packages = c("stats")
  ) %op% {
    n_cells <- as.integer(stats::runif(total_samples, min = min_cells, max = max_cells))
    beta_props_0 <- stats::rbeta(n_0, shape1_0, shape2_0)
    beta_props_1 <- stats::rbeta(n_1, shape1_1, shape2_1)
    beta_props <- c(beta_props_0, beta_props_1)
    
    successes <- stats::rbinom(total_samples, size = n_cells, prob = beta_props)
    prop <- successes / n_cells
    
    wilcox_res <- suppressWarnings(try(stats::wilcox.test(prop ~ group), silent = TRUE))
    
    if (!inherits(wilcox_res, "try-error")) {
      wilcox_res$p.value
    } else {
      NA_real_
    }
  }
  
  valid_pvals <- p_values[!is.na(p_values)]
  t1_error <- if (length(valid_pvals) > 0) mean(valid_pvals <= alpha) else NA_real_
  
  structure(
    list(
      type_1_error_rate = t1_error,
      convergences = length(valid_pvals),
      failed_fits = n_sims - length(valid_pvals),
      parameters = list(
        n_0 = n_0, n_1 = n_1, p_base = p_base, rho_0 = rho_0, rho_1 = rho_1,
        alpha = alpha, min_cells = min_cells, max_cells = max_cells, n_cores = n_cores
      )
    ),
    class = "t1_error_result_wilcox"
  )
}