#' Effort allocation
#'
#' \code{effort_allocation} redistributes fishing effort for the next timestep.
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param nm temporary numeric value, the current natural mortality estimate.
#' @param Allocation character value, how effort is to be allocated. Values can
#'    be:
#'    'IFD' - the ideal free distribution. Effort is allocated proportional to
#'       the yield caught in each area in the previous timestep.
#'    'equal' - Effort is allocated equally between all areas.
#'    Default value is 'IFD'.
#' @param E numeric array, the relative fishing effort displayed in each area,
#'    at each time step, under each control rule, and for each natural mortality
#'    estimate.
#' @param Yield numeric array, the yield caught in each area, at each timestep,
#'    under each control rule, and for each estimate of natural mortality, in kg.
#' @param Time1 numeric value, the number of years to run the model before a
#'    marine reserve is implemented. Default value is 50.
#' @param Inside numeric vector, the area(s) inside the marine reserve. Default
#'    value is c(3).
#' @param Outside numeric vector, the area(s) outside the marine reserve.
#'    Default value is c(1, 2, 4, 5).
#'
#' @return a numeric array updated with the fishing effort allocation for the
#'    next timestep.
#' @export
#'
#' @examples
#' E <- array(rep(1, 5*70*6*3), c(5, 70, 6, 3))
#' Yield <- array(rep(2458, 5*70*6*3), c(5, 70, 6, 3))
#' effort_allocation(t = 2, cr = 1, nm = 1, Allocation = 'IFD', E, Yield,
#'    Time1 = 50, Inside = c(3), Outside = c(1, 2, 4, 5))
effort_allocation <- function(t, cr, nm, Allocation = 'IFD', E, Yield,
                              Time1 = 50, Inside = c(3),
                              Outside = c(1, 2, 4, 5)) {

  # number of areas not in a reserve
  outs <- length(Outside)
  ins <- length(Inside)
  all <- outs + ins

  # If effort is allocated using the ideal free distribution, effort for one
  # year depends on the distribution of yield from the previous year
  if (Allocation == 'IFD') {

      prop_yield_out <- Yield[Outside, t - 1, cr, nm] /
        sum(Yield[Outside, t - 1, cr, nm])

      E[Outside, t, cr, nm] <- sum(E[, t - 1, cr, nm])*prop_yield_out
      E[Inside, t, cr, nm] <- 0

  # Otherwise, distribute effort equally between the four areas outside the
  # marine reserve, regardless of yield
  } else if (Allocation == 'equal') {

    if (t < Time1) {

      E[, t, cr] <- rep(sum(E[, t, cr, nm])/all, all)

    } else if (t >= Time1) {

    E[Outside, t, cr, nm] <- rep(sum(E[, t - 1, cr, nm])/outs, outs)
    E[Inside, t, cr, nm] <- 0

    }

  }

  return(E)

}
