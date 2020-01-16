#' Fishing mortality
#'
#' \code{fishing_mortality} returns a  numeric vector that represents the
#'    fishing mortality for all ages, on the interval (0, 1) in a specific area,
#'    at a specific time step, under a certain control rule, with a certain
#'    estimate of natural mortality
#'
#' @param a temporary numeric value, the current area
#' @param t temporary numeric value, the current time step
#' @param cr temporary numeric value, the current control rule
#' @param nm temporary numeric value, the current natural mortality estimate
#' @param FM numeric matrix, the values of fishing mortality for all previous
#'    areas, timesteps, control rules, and natural mortality estimates
#' @param A numeric value, the number of total areas in the model. The default
#'    value is 5.
#' @param Fb numeric value, the historical fishing effort for the fished species
#' @param E numeric matrix, the relative fishing effort displayed in each area,
#'    at each time step, under each control rule, and for each natural mortality
#'    estimate
#' @param S numeric vector, the selectivities at age from age at recruitment to
#'    maximum age, on the interval (0, 1)
#'
#' @return a numeric vector that corresponds to the fishing mortality at each
#'    age in a certain area, at a certain timestep, under a certain control
#'    rule, with a certain estimate of natural mortality
#' @export
#'
#' @examples
#' FM <- array(rep(0, 34*5*70*6*3), c(34, 5, 70, 6, 3))
#' E <- array(rep(0, 5*70*6*3), c(5, 70, 6, 3))
#' L <- length_age(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#' L2 = 47.95, K = 0.2022, All_ages = FALSE)
#' S <- selectivity(Rec_age = 2, Max_age = 35, L,
#' Fleets = c('sport', 'hook', 'trawl'), A50_up = c(2, 5, 10),
#' A50_down = c(6, 16, 35), Alpha = c(0.33, 0.6, 0.64),
#' F_fin = c(0.25, 0.06, 1), Beta = c(1.2, 0.6, 0), Cf = c(0.71, 0.28, 0.01))
#' fish_mort(a = 1, t = 1, cr = 1, nm = 1, FM, A = 5, Fb = 0.2, E, S)
fish_mort <- function(a, t, cr, nm, FM, A, Fb, E, S) {

  # Catchability (Vulnerability to fishing gear)
  # Based on Babcock & MacCall (2011): Eq. (6)
  vulnerability <- (A*Fb)/(sum(E[, 1, 1, 1]))

  # Selectivity as a matrix
  # dimensions = age * 1
  selectivity <- array(S, c(length(S), 1))

  # Effort as a matrix
  # Dimensions = area * time * CR
  effort <- E[a, t, cr, nm]

  # Fishing mortality
  # Based on Babcock & MacCall (2011): Eq. (5)
  # Dimensions = age * area * time * CR
  FM[, a, t, cr, nm] <- vulnerability * selectivity * effort

  return(FM[, a, t, cr, nm])

}
