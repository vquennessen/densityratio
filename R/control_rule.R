#' Control rule
#'
#' \code{control_rule} determines how to calculate the density ratio and manage
#'    the fishery (i.e. set fishing effort for the next timestep).
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param A numeric value, the number of total areas in the model. Default
#'    value is 5.
#' @param E numeric array, the relative fishing effort displayed in each area,
#'    at each time step, under each control rule, and for each natural mortality
#'    estimate.
#' @param Count numeric array, the number of individuals estimated to be in each
#'    area, at each timestep, under each control rule, for each estimate of
#'    natural mortality, for both all individuals and just mature individuals.
#' @param Time1 numeric value, the number of years to run the model before a
#'    marine reserve is implemented. Default value is 50.
#' @param TimeT numeric value, the number of years to run the model total.
#'    Default value is 70.
#' @param Transects numerical value, the number of sampling transects conducted
#'    in each area to estimate density ratio. Default value is 24.
#' @param Nat_mortality numeric vector, the estimates of natural mortality.
#' @param Final_DR numeric value, the final target density ratio.
#' @param Inside numeric vector, the area(s) inside the marine reserve. Default
#'    is c(3).
#' @param Outside numeric vector, the area(s) outside the marine reserve.
#'    Default is c(1, 2, 4, 5).
#' @param Areas_sampled character value, the areas to be sampled to calculate
#'    density ratio. Values can be:
#'    'all' - sample all areas.
#'    'far' - sample only areas farthest from the marine reserve.
#'    Default value is 'all'.
#' @param Ind_sampled character value, the individuals to be sampled to
#'    calculate density ratio. Values can be:
#'    'all' - sample all individuals.
#'    'mature' - sample only mature individuals.
#'    Default value is 'all'.
#' @param Years_sampled numeric value, the number of years of sampling upon
#'    which to base the estimate of density ratio. Default value is 1.
#'
#' @return a numeric vector of fishing effort for the next timestep, under the
#'    specific control rule, with a specific estimate of natural mortality.
#' @export
#'
#' @examples
#' E <- array(rep(1, 5*70*6*3), c(5, 70, 6, 3))
#' Count <- array(rep(5, 5*70*24*2*6*3), c(5, 70, 24, 2, 6, 3))
#' control_rule(t = 51, cr = 1, A = 5, E, Count, Time1 = 50, TimeT = 70,
#'    Transects = 24, Nat_mortality = c(0.09, 0.14, 0.19), Final_DR = 0.6,
#'    Inside = c(3), Outside = c(1, 2, 4, 5), Areas_sampled = 'all',
#'    Ind_sampled = 'all', Years_sampled = 1)
control_rule <- function(t, cr, A = 5, E, Count, Time1 = 50, TimeT = 70,
                         Transects = 24, Nat_mortality, Final_DR, Inside = c(3),
                         Outside = c(1, 2, 4, 5), Areas_sampled = 'all',
                         Ind_sampled = 'all', Years_sampled = 1) {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (!is.numeric(E)) {stop('E must be a numeric array.')}
  if (!is.numeric(Count)) {stop('Count must be a numeric array.')}
  if (Time1 %% 1 != 0) {stop('Time1 must be an integer value.')}
  if (TimeT %% 1 != 0) {stop('TimeT must be an integer value.')}
  if (Transects %% 1 != 0) {stop('Transects must be an integer value.')}
  if (!is.numeric(Nat_mortality)) {stop('Nat_mortality must be a numeric vector.')}
  if (!is.numeric(Final_DR)) {stop('Final_DR must be a numeric value.')}
  if (sum(Inside %% 1 != 0) != 0) {stop('Inside must be a vector of integers.')}
  if (sum(Outside %% 1 != 0) != 0) {stop('Outside must be a vector of integers.')}
  if (!is.character(Areas_sampled)) {stop('Areas_sampled must be a character value.')}
  if (!is.character(Ind_sampled)) {stop('Ind_sampled must be a character value.')}
  if (Years_sampled %% 1 != 0) {stop('Years_sampled must be an integer value.')}

  # acceptable values
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (A <= 0) {stop('A must be greater than 0.')}
  if (sum(E < 0) > 0) {stop('All values in E must be greater than or equal to 0.')}
  if (sum(Count < 0) > 0) {
    stop('All values in Count must be greater than or equal to 0.')}
  if (Time1 <= 0) {stop('Time1 must be greater than 0.')}
  if (TimeT <= 0) {stop('TimeT must be greater than 0.')}
  if (Transects <= 0) {stop('Transects must be greater than 0.')}
  if (sum(Nat_mortality <= 0) > 0 || sum(Nat_mortality > 1) > 0) {
    stop('All values in Nat_mortality must be between 0 and 1.')}
  if (Final_DR <= 0) {stop('Final_DR must be greater than 0.')}
  if (sum(Inside < 0) > 0) {
    stop('All values in Inside must be greater than or equal to 0.')}
  if (sum(Outside < 0) > 0) {
    stop('All values in Outside must be greater than or equal to 0.')}
  if (Areas_sampled != 'far' && Areas_sampled != 'all') {
    stop('Areas_sampled must be either "far" or "all".')}
  if (Ind_sampled != 'mature' && Ind_sampled != 'all') {
    stop('Ind_sampled must be either "mature" or "all".')}
  if (Years_sampled <= 0) {stop('Years_sampled must be greater than 0.')}

  # relational values
  if (sum(Inside > A) > 0) {
    stop('All values in Inside must be less than or equal to A.')}
  if (sum(Outside > A) > 0) {
    stop('All values in Outside must be less than or equal to A.')}
  if (sum(intersect(Inside, Outside)) > 0) {
    stop('Areas cannot both be inside and outside the marine reserve.')}
  if (Time1 >= TimeT) {stop('TimeT must be greater than Time1.')}
  if(dim(E)[1] != dim(Count)[1] || dim(E)[1] != A) {
    stop('E or Count has an incorrect number of areas.')}
  if(dim(E)[2] != dim(Count)[2] || dim(E)[2] != TimeT) {
    stop('E or Count has an incorrect number of time steps.')}
  if(dim(Count)[3] != Transects) {stop('Count has the wrong number of transects.')}
  if(dim(E)[3] != dim(Count)[5]) {
    stop('E or Count has an incorrect number of control rules.')}
  if(dim(E)[4] != dim(Count)[6]) {
    stop('E or Count has an incorrect number of natural mortality estimates.')}
  if (t > dim(E)[2]) {stop('The given "t" value is too high for E.')}
  if (cr > dim(E)[3]) {stop('The given "cr" value is too high for E.')}
  if (length(Nat_mortality) > dim(E)[4]) {
    stop('Incorrect number of natural mortality estimates.')}

  ##############################################################################

  if (cr < 4) { nm = cr } else { nm = cr - 3 }

  # static control rules, with constant target density ratios
  if (cr < 4) {

    DR <- density_ratio(t, cr, nm, A, Count, Years_sampled, Areas_sampled,
                        Ind_sampled, Transects, Inside, Outside)

    # calculate effort at the next timestep
    E[, t + 1, cr, ] <- management(t, cr, E, DR, target_DR = Final_DR,
                                   floor_DR = 0.2, effort_inc_allowed = 0.1,
                                   Time1)
  }

  # transient control rules with shifting target density ratios
  else {

    target_DR <- transient_DR(Time1, TimeT, Final_DR, Nat_mortality, nm)

    # calculate density ratio
    DR <- density_ratio(t, cr, nm, A, Count, Years_sampled, Areas_sampled,
                        Ind_sampled, Transects, Inside, Outside)

    # calculate effort at the next timestep
    E[, t + 1, cr, ] <- management(t, cr, E, DR, target_DR[t - Time1 + 1],
                                   floor_DR = 0.2, effort_inc_allowed = 0.1,
                                   Time1)
  }

  return(E[, t + 1, cr, ])

}
