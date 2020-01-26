#' Movement
#'
#' \code{movement} redistributes individuals amongst areas.
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param nm temporary numeric value, the current natural mortality estimate.
#' @param N numeric array, the number of individuals at each age, in each
#'    area, at each timestep, under each control rule, and for each estimate of
#'    natural mortality.
#' @param A numeric value, the number of total areas in the model. Default
#'    value is 5.
#' @param AMP numeric value, adult movement proportion, the fraction of
#'    individuals that move from one area to an adjacent area from one timestep
#'    to the next. Default value is 0.1.
#'
#' @return numeric array of updated N (numbers at age, area, timestep, control
#'    rule, and estimate of natural mortality)
#'
#' @examples
#' N <- array(rep(10, 34*5*70*6*3), c(34, 5, 70, 6, 3))
#' movement(t = 1, cr = 1, nm = 1, N, A = 5, AMP = 0.10)
movement <- function(t, cr, nm, N, A, AMP = 0.1) {

  # First area to second area
  N[, 1, t, cr, nm] <- (1 - AMP)*N[, 1, t, cr, nm] + AMP*N[, 2, t, cr, nm]

  # Intermediate areas to adjacent areas
  for (a in 2:(A-1)) {
    N[, a, t, cr, nm] <- (1 - 2*AMP)*N[, a, t, cr, nm] +
      AMP*(N[, a - 1, t, cr, nm] + N[, a + 1, t, cr, nm])
  }

  # Last area to next to last area
  N[, A, t, cr, nm] <- (1 - AMP)*N[, A, t, cr, nm] + AMP*N[, A - 1, t, cr, nm]

  return(N)

}
