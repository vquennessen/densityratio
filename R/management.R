#' Management
#'
#' \code{management} calculates the next timestep's fishing effort allowed,
#'    given the control rule and estimated density ratio.
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param E numeric array, the relative fishing effort displayed in each area,
#'    at each time step, under each control rule, and for each natural mortality
#'    estimate.
#' @param DR numeric value, the calculated density ratio at a particular
#'    timestep and under a particular control rule.
#' @param target_DR numeric value, the target density ratio, below which effort
#'    is allowed to increase, and above which effort is decreased for the next
#'    timestep.
#' @param floor_DR numeric value, the density ratio below which effort is set
#'    back to 10\% of the original value. Default value is 0.20.
#' @param effort_inc_allowed numeric value, the percent increase allowed in
#'    effort if the density ratio is below the target density ratio, on the
#'    interval (0, 1). Default value is 0.10.
#' @param Time1 numeric value, the number of years to run the model before a
#'    marine reserve is implemented. Default value is 50.
#'
#' @return a numeric vector of fishing effort for the next timestep, under the
#'    specific control rule, with a specific estimate of natural mortality.
#'
#' @examples
#' E <- array(rep(1, 5*70*6*3), c(5, 70, 6, 3))
#' management(t = 51, cr = 1, E, DR = 0.8, target_DR = 0.6, floor_DR = 0.2,
#'    effort_inc_allowed = 0.10, Time1 = 50)

management <- function(t, cr, E, DR, target_DR, floor_DR = 0.2,
                       effort_inc_allowed = 0.10, Time1 = 50) {

  # If the control rule is based on effort, and the density ratio is higher
  # than the target density ratio, allow effort in each area to increase by the
  # allowed effort increase value (typically 10%)
  if (DR > target_DR) {

    E[, t + 1, cr, ] <- E[, t, cr, ]*(1 - effort_inc_allowed)

  # If the density ratio is lower than the target density ratio but
  # greater than the floor density ratio, allow effort in each area to decrease
  # by the allowed effort increase value (typically 10%)
  } else if (DR <= target_DR & DR > floor_DR) {

    E[, t + 1, cr, ] <- E[, t, cr, ]*(1 + effort_inc_allowed)

  # Finally, if the density ratio is below the floor density ratio, effort is
  # decreased back down to 10% of the original value
  } else if (DR <= floor_DR) {

    E[, t + 1, cr, ] <- E[, Time1, cr, ]*0.10

  }

  return(E[, t + 1, cr, ])

}
