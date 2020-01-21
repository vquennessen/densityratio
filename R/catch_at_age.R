#' Catch at age
#'
#' \code{catch} returns the catch at age vector. It is based on the continuous
#' (Baranov) formulation of of selectivity, such that fishing mortality (FM) and
#' natural mortality (M) are proportional to fish abundance (N) and act
#' simultaneously  and uniformly throughout the year, i.e. dN/dt = -(M+F)*N.
#'
#' @param a temporary numeric value, the current area.
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param nm temporary numeric value, the current natural mortality estimate.
#' @param FM numeric array that corresponds to the fishing mortality at each
#'    age in each area, at each timestep, under all control rules, with all
#'    estimates of natural mortality.
#' @param Nat_mortality numeric vector, the estimates of natural mortality.
#' @param N numeric array, the number of individuals at each age, in each
#'    area, at each timestep, under each control rule, and for each estimate of
#'    natural mortality.
#' @param A numeric value, the number of total areas in the model. Default
#'    value is 5.
#' @param Fb numeric value, the historical fishing effort for the fished species.
#' @param E numeric array, the relative fishing effort displayed in each area,
#'    at each time step, under each control rule, and for each natural mortality
#'    estimate.
#' @param Catch numeric array, the number of individuals caught at each age, in
#'    each area, at each timestep, under each control rule, for each estimate of
#'    natural mortality.
#'
#' @return a numeric array with an updated vector of catch at ages given the
#'    current area, timestep, control rule, and estimate of natural mortality.
#' @export
#'
#' @examples
#' FM <- array(rep(0.2, 34*5*70*6*3), c(34, 5, 70, 6, 3))
#' N <- array(rep(10, 34*5*70*6*3), c(34, 5, 70, 6, 3))
#' E <- array(rep(1, 5*70*6*3), c(5, 70, 6, 3))
#' Catch <- array(rep(2, 34*5*70*6*3), c(34, 5, 70, 6, 3))
#' catch(a = 1, t = 1, cr = 1, nm = 1, FM, Nat_mortality = c(0.09, 0.14, 0.19),
#'    N, A = 5, Fb = 0.2, E, Catch)

catch <- function(a, t, cr, nm, FM, Nat_mortality, N, A, Fb, E, Catch) {

  # calculate the coefficient
  coeff <- FM[ , a, t, cr, nm]/(Nat_mortality[nm] + FM[ , a, t, cr, nm])

  # calculate catch at age
  Catch[ , a, t, cr, nm] <- coeff * N[ , a, t, cr, nm] * exp(-1*Nat_mortality[nm]
                                                             - FM[ , a, t, cr, nm])

  return(Catch[, a, t, cr, nm])

}
