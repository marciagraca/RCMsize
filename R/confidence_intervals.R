#' Confidence Interval for Seroprevalence
#'
#' Calculates the confidence interval for a seroprevalence estimate with a specified confidence level.
#'
#' @param SP Seroprevalence estimate.
#' @param n Sample size.
#' @param conf.level Confidence level (default is 0.95).
#' @param method Method for calculating the confidence interval (default is "asymptotic"). Available methods: c("asymptotic","exact","ac","wilson","logit","cloglog")
#' @return A vector with the lower and upper limits of the confidence interval.
#' @examples
#' IC_SP(0.25, 100, conf.level = 0.95, method = "asymptotic")
#' @references
#' The methods available in this function are some of the available in the binom package. For more information, see \url{https://CRAN.R-project.org/package=binom }
#' @importFrom binom binom.confint
#' @export

IC_SP <- function(SP, n, conf.level = 0.95, method = "asymptotic"){
  x <- SP*n
  info <- binom::binom.confint(x,n,conf.level = conf.level,methods = method)
  ll <- info$lower
  ul <- info$upper
  return(c(ll,ul))
}
#' Confidence Interval for Seroprevalence with Continuity Correction (Wald Method)
#'
#' Calculates the confidence interval for seroprevalence using the Wald method with continuity correction.
#'
#' @param SP Seroprevalence estimate.
#' @param n Sample size.
#' @param conf.level Confidence level (default is 0.95).
#' @return A vector with the lower and upper limits of the confidence interval.
#' @examples
#' IC_SP_Waldcc(0.25, 100, conf.level = 0.95)
#' @importFrom stats qnorm
#' @export

IC_SP_Waldcc <- function(SP,n, conf.level=0.95){
  alpha <- 1 - conf.level
  z <- stats::qnorm(1 - alpha / 2)
  ME <- z*sqrt(SP*(1-SP)/n)
  ll <- SP-ME-1/(2*n)
  ul <- SP+ME+1/(2*n)
  return(c(ll,ul))
}

#' Confidence Interval for the Seroconversion Rate (SCR)
#'
#' Calculates the confidence interval for the seroconversion rate (SCR) using
#' the confidence interval of seroprevalence.
#'
#' @param SP_interval A vector with the lower and upper limits of seroprevalence.
#' @param SRR Seroreversion rate.
#' @param ages Vector with the proportions of different ages in the population (age structure).
#' @param A_max Maximum age considered in the population.
#' @param limits Lower and upper limits for the calculation of \code{SCR}.
#' @return A vector with the lower and upper limits for the seroconversion rate \code{SCR}.
#' @examples
#' A_max <- 80
#' age_distribution <- rep(1 / A_max, A_max)
#' IC_SCR(c(0.1, 0.2), 0.01, age_distribution, A_max, limits = c(0, 1))
#' @export

IC_SCR <- function(SP_interval,SRR,ages,A_max,limits=c(0,1)){
  SCR_ll <-stats::uniroot(function(SCR) seroprevalence(ages,A_max,SCR,SRR) - SP_interval[1], lower = limits[1], upper = limits[2], tol = 1e-5)$root
  SCR_ul <- stats::uniroot(function(SCR) seroprevalence(ages,A_max,SCR,SRR) - SP_interval[2], lower = limits[1], upper = limits[2], tol = 1e-5)$root
  return(c(SCR_ll,SCR_ul))
}
