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

  ###### Error handling ########################################################

  # classes of variables
  if (a %% 1 != 0) {stop('a must be an integer value.')}
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (nm %% 1 != 0) {stop('nm must be an integer value.')}
  if (!is.numeric(Delta)) {stop('Delta must be a numeric value.')}
  if (!is.numeric(Gamma)) {stop('Gamma must be a numeric value.')}
  if (!is.numeric(Abundance_all)) {
    stop('Abundance_all must be a numeric array.')}
  if (!is.numeric(Abundance_mature)) {
    stop('Abundance_mature must be a numeric array.')}
  if (Transects %% 1 != 0) {stop('Transects must be an integer value.')}
  if (!is.numeric(X)) {stop('X must be a numeric value.')}
  if (!is.numeric(Count)) {stop('Count must be a numeric array.')}
  if (!is.numeric(NuS)) {stop('NuS must be a numeric array.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}

  # acceptable values
  if (a <= 0) {stop('a must be greater than 0.')}
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (nm <= 0 || nm > 3) {
    stop('nm must be greater than 0 and less than or equal to 3.')}
  if (Delta <= 0) {stop('Delta must be greater than 0.')}
  if (Gamma <= 0) {stop('Gamma must be greater than 0.')}
  if (sum(Abundance_all < 0) > 0) {
    stop('All values in Abundance_all must be greater than or equal to 0.')}
  if (sum(Abundance_mature < 0) > 0) {
    stop('All values in Abundance_mature must be greater than or equal to 0.')}
  if (Transects <= 0) {stop('Transects must be greater than 0.')}
  if (X < 0) {stop('X must be greater than or equal to 0.')}
  if (sum(Count < 0) > 0) {
    stop('All values in Count must be greater than or equal to 0.')}
  if (A <= 0) {stop('A must be greater than 0.')}

  # relational values
  if (dim(Abundance_all)[1] != dim(Abundance_mature)[1] || dim(Abundance_all)[1] != A) {
    stop('Abundance_all or Abundance_mature has an incorrect number of areas.')}
  if (dim(Abundance_all)[2] != dim(Abundance_mature)[2]) {
    stop('Abundance_all or Abundance_mature has an incorrect number of time steps.')}
  if (dim(Abundance_all)[3] != dim(Abundance_mature)[3]) {
    stop('Abundance_all or Abundance_mature has an incorrect number of control rules.')}
  if (dim(Abundance_all)[4] != dim(Abundance_mature)[4]) {
    stop('Abundance_all or Abundance_mature has an incorrect number of natural
         mortality estimates.')}
  if (dim(Count)[1] != dim(NuS)[1] || dim(Abundance_all)[1] != A) {
    stop('Count or NuS has an incorrect number of areas.')}
  if (dim(Count)[2] != dim(NuS)[2]) {
    stop('Count or NuS has an incorrect number of time steps.')}
  if (dim(Count)[3] != Transects) {stop('Count has the wrong number of transects.')}
  if (dim(Count)[5] != dim(NuS)[3]) {
    stop('Count or NuS has an incorrect number of control rules.')}
  if (dim(Count)[6] != dim(NuS)[4]) {
    stop('Count or NuS has an incorrect number of natural mortality estimates.')}
  if (a > A) {stop('The given "a" value is too high.')}
  if (t > dim(Abundance_all)[2] || t > dim(Count)[2]) {
    stop('The given "t" value is too high for Abundance_all or Count.')}
  if (cr > dim(Abundance_all)[3]|| t > dim(Count)[3]) {
    stop('The given "cr" value is too high for Abundance_all or Count.')}
  if (nm > dim(Abundance_all)[4] || t > dim(Count)[4]) {
    stop('The given "nm" value is too high for Abundance_all or Count.')}

  ##############################################################################

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
