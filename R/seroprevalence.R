#' Calculation of Seropositivity Probability
#'
#' This function calculates the probability of seropositivity based on the age and the seroconversion and seroreversion rates,
#'  using a reversible catalytic model.
#'
#' @param SCR Seroconversion Rate
#' @param SRR Seroreversion Rate.
#' @param t Age for which we want to calculate the probability of seropositivity.
#' @return The probability of seropositivity for age `t`.
#' @examples
#' prob_seropositive(0.03, 0.01, 45)
#' @references
#' For more information on the reversible catalytic model, see \url{https://rdcu.be/d0d69}
#' @export

prob_seropositive <- function(SCR,SRR,t){
  gamma <- SCR+SRR
  ps <- (SCR/gamma)*(1-exp(-gamma*t))
  return (ps)
}
#' Seroprevalence Calculation
#'
#' Calculates the seroprevalence considering an age distribution and a reversible
#' catalytic model.
#'
#' @param ages Vector with the proportions of different ages in the population (age structure).
#' @param A_max Maximum age considered in the population.
#' @param SCR Seroconversion rate.
#' @param SRR Seroreversion rate.
#' @return The total seroprevalence weighted by the age distribution.
#' @examples
#' A_max <- 80
#' age_distribution <- rep(1 / A_max, A_max)
#' seroprevalence(age_distribution, A_max, 0.03, 0.01)

#' @export

seroprevalence <- function(ages, A_max, SCR,SRR){
  SP <- 0
  for (i in 1:A_max){
    SP = SP + ages[i] * prob_seropositive(SCR,SRR,i)
  }
  return (SP)
}
