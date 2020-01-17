#' Recruitment
#'
#' \code{recruitment} calculates the number of new recruits entering the
#'    population
#'
#' @param a temporary numeric value, the current area
#' @param t temporary numeric value, the current time step
#' @param cr temporary numeric value, the current control rule
#' @param nm temporary numeric value, the current natural mortality estimate
#' @param SSB numeric matrix, gives the spawning stock biomass of the whole
#'    stock for each area, at each timestep, under each control rule, and for
#'    each estimate of natural mortality
#' @param A numeric value, the number of total areas in the model. The default
#'    value is 5.
#' @param R0 numeric value, set arbitrarily, the unfished recruitment. The
#'    default value is 1e+5.
#' @param H numeric value, the steepness of the stock-recruitment curve
#' @param B0 numeric value, set arbitrarily, the unfished biomass
#' @param Eps numeric matrix, recruitment error terms
#' @param Sigma_R numeric value, recruitment standard deviation
#' @param Rec_age numeric value, gives age at which fish has entered the fishery
#'    in years
#' @param Recruitment_mode character value, values can be 'closed' (if the
#'    recruits in each area originate from adults in that area) or 'pool' (if
#'    the recruits in each area come from a pool of larvae produced by all
#'    reproducing individuals in all areas) - the default value is 'pool'
#'
#' @return a numeric value representing the number of new recruits coming into
#'    the population in area a, at timestep t, under control rule cr, with an
#'    estimate of natural mortality of nm
#' @export
#'
#' @examples
#' SSB <- array(rep(10, 5*70*6*3), c(5, 70, 6, 3))
#' NuR <- array(rnorm(5*70*6*3, 0, 0.5), c(5, 70, 6, 3))
#' Eps <- epsilon(A = 5, TimeT = 70, CR = 6, NM = 3, NuR, Rho_R = 0)
#' recruitment(a = 1, t = 3, cr = 1, nm = 2, SSB, A = 5, R0 = 1e+5, H = 0.7,
#'    B0 = 1e+5/1.1, Eps, Sigma_R = 0.5, Rec_age = 2, Recruitment_mode = 'pool')
recruitment = function(a, t, cr, nm, SSB, A = 5, R0 = 1e+5, H, B0, Eps, Sigma_R,
                       Rec_age, Recruitment_mode) {

  # Recruitment
  # Based on Babcock & MacCall (2011): Eq. (3)
  # Dimensions = 1 * 1

  # closed recruitment indicates that recruits are only offspring of adults from
  # the same area
  if (Recruitment_mode == 'closed') { ssb <- SSB[a, t - Rec_age, cr, nm]}

  # pool recruitment indicates that recruits come from a larval pool produced by
  # adults from all areas
  if (Recruitment_mode == 'pool') { ssb <- sum(SSB[, t - Rec_age, cr, nm] / A)}

  adjR0 <- R0 / A
  adjB0 <- B0 / A

  R1 <- (0.8 * adjR0 * H * ssb) / (0.2 * adjB0 * (1 - H) + (H - 0.2) * ssb)
  recruits <- R1 * (exp(Eps[a, t, cr, nm] - Sigma_R^2 / 2))

  return(recruits)

}
