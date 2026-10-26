#' Estimate Overdispersion (Rho) via Method of Moments
#'
#' Calculates a rapid empirical estimate of the intra-class correlation
#' (overdispersion parameter, rho) for proportion data. This method of moments
#' estimator uses the ratio of the sample variance of the proportions to the
#' theoretical maximum variance of a true binomial distribution with the same mean.
#'
#' @param props Numeric vector of observed proportions (must be between 0 and 1).
#' @param na.rm Logical. Should missing values (including NaN) be removed before
#'   computation? (default: TRUE).
#'
#' @return A numeric value representing the estimated overdispersion parameter (rho).
#'   Returns `NA` if the sample mean is exactly 0 or 1, as the denominator becomes 0.
#'
#' @importFrom stats var na.omit
#' @export
#'
#' @examples
#' # Simulate 100 proportions with a known mean of 0.2 and overdispersion of 0.05
#' set.seed(123)
#' p_base <- 0.2
#' rho_true <- 0.05
#' shape1 <- p_base * (1 - rho_true) / rho_true
#' shape2 <- (1 - p_base) * (1 - rho_true) / rho_true
#'
#' sim_props <- stats::rbeta(100, shape1, shape2)
#'
#' # Estimate rho empirically
#' estimated_rho <- calc_rho_mom(sim_props)
#' print(estimated_rho)
calc_rho_mom <- function(props, na.rm = TRUE) {
  if (na.rm) {
    props <- stats::na.omit(props)
  }

  if (anyNA(props)) {
    return(NA_real_)
  }

  if (length(props) < 2) {
    stop("At least two observations are required to calculate sample variance.")
  }

  if (any(props < 0 | props > 1)) {
    stop("All values in 'props' must be proportions between 0 and 1.")
  }

  p_bar <- mean(props)

  # Prevent division by zero if all proportions are exactly 0 or exactly 1
  if (p_bar == 0 || p_bar == 1) {
    warning("Sample mean is exactly 0 or 1. Rho is undefined (division by zero). Returning NA.")
    return(NA_real_)
  }

  s2 <- stats::var(props)
  rho <- s2 / (p_bar * (1 - p_bar))

  if (rho < 0) {
    warning("Estimated rho is negative. Data may be underdispersed or sample size is too small.")
  }

  return(rho)
}
