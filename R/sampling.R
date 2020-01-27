#' Sampling
#'
#' \code{sampling} samples from the current population to update the Count
#'    array.
#'
#' @param a temporary numeric value, the current area.
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param nm temporary numeric value, the current natural mortality estimate.
#' @param Delta numeric value, the proportion of positive transects divided by
#'    depletion, also known as the constant of proportionality
#' @param Gamma numeric value, the average value of a positive transect divided
#'    by depletion
#' @param Abundance_all numeric array, the total number of individuals in each
#'    area, at each timestep, under all control rules, with all estimates of
#'    natural mortality.
#' @param Abundance_mature numeric array, the number of mature individuals in
#'    each area, at each timestep, under all control rules, with all estimates
#'    of natural mortality.
#' @param Transects numerical value, the number of sampling transects conducted
#'    in each area to estimate density ratio. Default value is 24.
#' @param X numeric value, the average value of individuals seen during positive
#'    transects.
#' @param Count numeric array, the number of individuals estimated to be in each
#'    area, at each timestep, under each control rule, for each estimate of
#'    natural mortality, for both all individuals and just mature individuals.
#' @param NuS numeric array, the sampling normal variable pulled from a normal
#'    distribution with mean 0 and sd Sigma_S.
#' @param A numeric value, the number of total areas in the model. Default value
#'    is 5.
#'
#' @return a two-dimentional numeric array with updated observed counts for the
#'    fished species in this particular area, at this particular timestep,
#'    under this particular control rule, with this particular estimate of
#'    natural mortality.
#' @export
#'
#' @importFrom stats rbinom
#'
#' @examples
#' Abundance_all <- array(rep(340, 5*70*6*3), c(5, 70, 6, 3))
#' Abundance_mature <- array(rep(280, 5*70*6*3), c(5, 70, 6, 3))
#' Count <- array(rep(50, 5*70*24*2*6*3), c(5, 70, 24, 2, 6, 3))
#' NuS <- array(stats::rnorm(5*70*6*3, 0, 0.89), c(5, 70, 6, 3))
#' sampling(a = 1, t = 51, cr = 1, nm = 1, Delta = 1.6, Gamma = 31.6,
#'    Abundance_all, Abundance_mature, Transects = 24, X = 15.42, Count, NuS,
#'    A = 5)
sampling <- function(a, t, cr, nm, Delta, Gamma, Abundance_all,
                     Abundance_mature, Transects = 24, X, Count, NuS, A = 5) {

  # Total population size across all areas
  total_all <- sum(Abundance_all[, t, cr, nm])
  total_mature <- sum(Abundance_mature[, t, cr, nm])

  # Calculate odds ratio of seeing a fish
  # Based on Babcock & MacCall (2011): Eq. (12)
  odds_all <-  (Delta * Abundance_all[, t, cr, nm]) / (total_all / A)
  odds_mature <-  (Delta * Abundance_mature[, t, cr, nm]) / (total_mature / A)

  # Calculate probability based on odds ratio
  p_all <- 1 / (1 + exp(odds_all))
  p_mature <- 1 / (1 + exp(odds_mature))

  # Determine if species is seen at least once
  # Dimensions = 1 * transects
  presence_all <- rbinom(Transects, 1, p_all)
  presence_mature <- rbinom(Transects, 1, p_mature)

  # Calculate species count given transects with positive visuals
  Count[a, t, , 1, cr, nm] <- presence_all*(Gamma*Abundance_all[a, t, cr, nm]*exp(NuS[a, t, cr, nm]))
  Count[a, t, , 2, cr, nm] <- presence_mature*(Gamma*Abundance_mature[a, t, cr, nm]*exp(NuS[a, t, cr, nm]))

  return(Count[a, t, , , cr, nm])
}
