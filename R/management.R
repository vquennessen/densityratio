#' Management
#'
#' \code{management} calculates the next timestep's fishing effort allowed,
#'    given the control rule and estimated density ratio.
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param fdr temporary numeric value, the current final target density ratio.
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
#' @export
#'
#' @examples
#' A = 5; TimeT = 70; CR = 6;  FDR = 4
#' E <- array(rep(0.2, A*TimeT*CR*FDR), c(A, TimeT, CR, FDR))
#' management(t = 51, cr = 1, fdr = 1, E, DR = 0.8, target_DR = 0.6,
#'    floor_DR = 0.2, effort_inc_allowed = 0.10, Time1 = 50)

management <- function(t, cr, fdr, E, DR, target_DR, floor_DR = 0.2,
                       effort_inc_allowed = 0.10, Time1) {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (fdr %% 1 != 0) {stop('fdr must be an integer value.')}
  if (!is.numeric(E)) {stop('E must be a numeric array.')}
  if (!is.numeric(DR)) {stop('DR must be a numeric value.')}
  if (!is.numeric(target_DR)) {stop('target_DR must be a numeric value.')}
  if (!is.numeric(floor_DR)) {stop('floor_DR must be a numeric value.')}
  if (!is.numeric(effort_inc_allowed)) {
    stop('effort_inc_allowed must be a numeric value.')}
  if (Time1 %% 1 != 0) {stop('Time1 must be an integer value.')}

  # acceptable values
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (fdr <= 0) {stop('fdr must be greater than 0.')}
  if (sum(E < 0) > 0) {stop('All values in E must be greater than or equal to 0.')}
  if (DR <= 0) {stop('DR must be greater than 0.')}
  if (target_DR <= 0) {stop('target_DR must be greater than 0.')}
  if (floor_DR < 0) {stop('floor_DR must be greater than or equal to 0.')}
  if (effort_inc_allowed < 0) {
    stop('effort_inc_allowed must be greater than or equal to 0.')}
  if (Time1 <= 0) {stop('Time1 must be greater than 0.')}

  # relational values
  if (t > dim(E)[2]) {stop('The given "t" value is too high for E.')}
  if (cr > dim(E)[3]) {stop('The given "cr" value is too high for E.')}
  if (fdr > dim(E)[4]) {stop('The given "fdr" value is too high for E.')}

  ##############################################################################

  # If the density ratio is higher than the target density ratio, allow effort
  # in each area to increase by the amount allowed (typically 10%)
  if (DR > target_DR) {

    E[, t + 1, cr, fdr] <- E[, t, cr, fdr]*(1 + effort_inc_allowed)

  # If the density ratio is lower than the target density ratio but
  # greater than the floor density ratio, allow effort in each area to decrease
  # by the allowed effort increase value (set to 10% here)
  } else if (DR <= target_DR & DR > floor_DR) {

    E[, t + 1, cr, fdr] <- 0.90 * E[, t, cr, fdr]

  # Finally, if the density ratio is below the floor density ratio, effort is
  # decreased back down to 10% of the original value
  } else if (DR <= floor_DR) {

    E[, t + 1, cr, fdr] <- E[, Time1, cr, fdr]*0.10

  }

  return(E[, t + 1, cr, fdr])

}
