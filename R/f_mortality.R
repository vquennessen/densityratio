#' Fishing mortality
#'
#' \code{f_mortality} returns a  numeric vector that represents the
#'    fishing mortality for all ages, on the interval (0, 1) in a specific area,
#'    at a specific time step, under a certain control rule, with a certain
#'    estimate of natural mortality
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param nm temporary numeric value, the current natural mortality estimate.
#' @param fdr temporary numeric value, the current final target density ratio.
#' @param FM numeric array, the values of fishing mortality for all ages,
#'    areas, timesteps, control rules, and natural mortality estimates.
#' @param A numeric value, the number of total areas in the model. Default
#'    value is 5.
#' @param Fb numeric value, the historical fishing effort for the fished species.
#' @param E numeric matrix, the relative fishing effort displayed in each area,
#'    at each time step, under each control rule, and for each natural mortality
#'    estimate.
#' @param S numeric vector, the selectivities at age from age at recruitment to
#'    maximum age, on the interval (0, 1).
#'
#' @return a numeric vector that corresponds to the fishing mortality at each
#'    age in a certain area, at a certain timestep, under a certain control
#'    rule, with a certain estimate of natural mortality.
#' @export
#'
#' @examples
#' n = 34; A = 5; TimeT = 70; CR = 6; NM = 3; FDR = 4
#' FM <- array(rep(0, n*A*TimeT*CR*NM*FDR), c(n, A, TimeT, CR, NM, FDR))
#' E <- array(rep(0.2, A*TimeT*CR*NM*FDR), c(A, TimeT, CR, NM, FDR))
#' L <- length_age(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#'    L2 = 47.95, K = 0.2022, All_ages = FALSE)
#' S <- selectivity(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#'    L2 = 47.95, K = 0.2022, Fleets = c('sport', 'hook', 'trawl'),
#'    A50_up = c(2, 5, 10), A50_down = c(6, 16, 35), Alpha = c(0.33, 0.6, 0.64),
#'    F_fin = c(0.25, 0.06, 1), Beta = c(1.2, 0.6, 0), Cf = c(0.71, 0.28, 0.01))
#' f_mortality(t = 1, cr = 1, nm = 1, fdr = 1, FM, A = 5, Fb = 0.2, E, S)
f_mortality <- function(t, cr, nm, fdr, FM, A, Fb, E, S) {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (nm %% 1 != 0) {stop('nm must be an integer value.')}
  if (fdr %% 1 != 0) {stop('fdr must be an integer value.')}
  if (!is.numeric(FM)) {stop('FM must be a numeric array.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (!is.numeric(Fb)) {stop('Fb must be a numeric value.')}
  if (!is.numeric(E)) {stop('E must be a numeric array.')}
  if (!is.numeric(S)) {stop('S must be a numeric vector.')}

  # acceptable values
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (nm <= 0 || nm > 3) {
    stop('nm must be greater than 0 and less than or equal to 3.')}
  if (fdr <= 0) {stop('fdr must be greater than 0.')}
  if (sum(FM < 0) > 0) {
    stop('All values in FM must be greater than or equal to 0.')}
  if (A <= 0) {stop('A must be greater than 0.')}
  if (Fb < 0) {stop('Fb must be greater than or equal to 0.')}
  if (sum(E < 0) > 0) {stop('All values in E must be greater than or equal to 0.')}
  if (sum(S < 0) > 0 || sum(S > 1) > 0) {
    stop('All values in S must be between 0 and 1.')}

  # relational values
  if(dim(FM)[1] != length(S)) {
    stop('FM or S has an incorrect number of age classes.')}
  if(dim(FM)[2] != dim(E)[1] || dim(FM)[2] != A) {
    stop('FM or E has an incorrect number of areas.')}
  if(dim(FM)[3] != dim(E)[2]) {
    stop('FM or E has an incorrect number of time steps.')}
  if(dim(FM)[4] != dim(E)[3]) {
    stop('FM or E has an incorrect number of control rules.')}
  if(dim(FM)[5] != dim(E)[4]) {
    stop('FM or E has an incorrect number of natural mortality estimates.')}
  if(dim(FM)[6] != dim(E)[5]) {
    stop('FM or E has an incorrect number of final density ratios.')}
  if (t > dim(FM)[3]) {stop('The given "t" value is too high for FM.')}
  if (cr > dim(FM)[4]) {stop('The given "cr" value is too high for FM.')}
  if (nm > dim(FM)[5]) {stop('The given "nm" value is too high for FM.')}
  if (fdr > dim(FM)[6]) {stop('The given "fdr" value is too high for FM.')}

  ##############################################################################

  # Catchability (Vulnerability to fishing gear)
  # Based on Babcock & MacCall (2011): Eq. (6)
  vulnerability <- (A*Fb)/(sum(E[, 1, 1, 1, 1]))

  # Selectivity as a matrix
  # dimensions = age * 1
  selectivity <- array(S, c(length(S), 1))

  # Effort as a matrix
  # Dimensions = area * time * CR * NM * FDR
  effort <- E[, t, cr, nm, fdr]

  # Fishing mortality
  # Based on Babcock & MacCall (2011): Eq. (5)
  # Dimensions = age * area * time * CR
  FM[, , t, cr, nm, fdr] <- vulnerability * selectivity %*% effort

  return(FM[, , t, cr, nm, fdr])

}
