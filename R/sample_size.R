#' Sample Size Calculation
#'
#' Estimates the required sample size so that the confidence interval width
#' for SCR does not exceed a specified limit.
#'
#' @param SCR Seroconversion rate.
#' @param RL Desired relative width.
#' @param SRR Seroreversion rate.
#' @param ages Vector with the proportions of different ages in the population (age structure).
#' @param A_max Maximum age considered in the population.
#' @param limits Lower and upper limits for the calculation of \code{SCR}.
#' @param max_iter Maximum number of iterations.
#' @param conf.level Confidence level (default is 0.95).
#' @param method Method for calculating the confidence interval. Available methods: "waldcc" and the methods in IC_SP documentation.
#' @return A list with the required sample size, the confidence interval for
#' seroprevalence, and the confidence interval for \code{SCR}.
#' @examples
#' A_max <- 80
#' age_distribution <- rep(1 / A_max, A_max)
#' sample_s(0.03, 1, 0.01, age_distribution, A_max, limits = c(0, 1))
#' @details
#' **Disclaimer**: The sample size function may not produce accurate values for
#' scenarios involving extremely low SCR (e.g.,
#' elimination scenarios). Users are advised to exercise caution and consider
#' the results critically when applying this function to such cases.
#'
#' @importFrom stats uniroot
#' @export

sample_s <- function(SCR, RL, SRR, ages, A_max, limits,max_iter = 10000,conf.level=0.95,method="asymptotic") {
  if (RL <= 0) {
    stop("RL should be a positive value")
  }
  ampl <- RL*SCR
  n <- 50
  sp <- seroprevalence(ages, A_max, SCR, SRR)
  if (method == "waldcc") {
    ic_sp <- IC_SP_Waldcc(sp, n, conf.level)
  } else {
    ic_sp <- IC_SP(sp, n, conf.level, method)
  }

  ic_scr <- IC_SCR(ic_sp, SRR, ages, A_max, limits)
  amplitude <- ic_scr[2] - ic_scr[1]

  iter <- 0

  while (amplitude > ampl && iter < max_iter) {
    n <- n + 1
    if (method == "waldcc") {
      ic_sp <- IC_SP_Waldcc(sp, n, conf.level)
    } else {
      ic_sp <- IC_SP(sp, n, conf.level, method)
    }
    ic_scr <- IC_SCR(ic_sp, SRR, ages, A_max, limits)
    amplitude <- ic_scr[2] - ic_scr[1]
    iter <- iter + 1
  }

  if (iter == max_iter) {
    warning("Maximum number of iteration reached")
  }
  return(list(n = n, ci_sp = ic_sp, ci_scr = ic_scr))
}
