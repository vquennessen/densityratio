#' Catch at age
#'
#' \code{catch} returns the catch at age vector. It is based on the continuous
#' (Baranov) formulation of of selectivity, such that fishing mortality (FM) and
#' natural mortality (M) are proportional to fish abundance (N) and act
#' simultaneously  and uniformly throughout the year, i.e. dN/dt = -(M+F)*N.
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param ENM numeric value, the nm value that represents the 'true' population.
#' @param fdr temporary numeric value, the current final target density ratio.
#' @param FM numeric array that corresponds to the fishing mortality at each
#'    age in each area, at each timestep, under all control rules, with all
#'    estimates of natural mortality.
#' @param M numeric value, the natural mortality on the interval (0, 1).
#' @param N numeric array, the number of individuals at each age, in each
#'    area, at each timestep, under each control rule, and for each estimate of
#'    natural mortality.
#' @param A numeric value, the number of total areas in the model. Default
#'    value is 5.
#' @param Fb numeric value, the historical fishing effort for the fished species.
#' @param E numeric array, the relative fishing effort displayed in each area,
#'    at each time step, under each control rule, and for each natural mortality
#'    estimate.
#' @param Catch numeric array, the number of individuals caught at each age, in
#'    each area, at each timestep, under each control rule, for each estimate of
#'    natural mortality.
#'
#' @return a numeric array with an updated vector of catch at ages given the
#'    current area, timestep, control rule, and estimate of natural mortality.
#' @export
#'
#' @examples
#' n = 34; A = 5; TimeT = 70; CR = 6; FDR = 4; NM = 3
#' FM <- array(rep(0.2, n*A*TimeT*CR*FDR), c(n, A, TimeT, CR, FDR))
#' N <- array(rep(10, n*A*TimeT*CR*FDR*NM), c(n, A, TimeT, CR, FDR, NM))
#' E <- array(rep(1, A*TimeT*CR*FDR), c(A, TimeT, CR, FDR))
#' Catch <- array(rep(2, n*A*TimeT*CR*FDR), c(n, A, TimeT, CR, FDR))
#' catch(t = 1, cr = 1, ENM = 2, fdr = 1, FM, M = 0.14, N, A = 5, Fb = 0.2,
#'    E, Catch)

catch <- function(t, cr, ENM, fdr, FM, M, N, A, Fb, E, Catch) {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (ENM %% 1 != 0) {stop('ENM must be an integer value.')}
  if (fdr %% 1 != 0) {stop('fdr must be an integer value.')}
  if (!is.numeric(FM)) {stop('FM must be a numeric array.')}
  if (!is.numeric(M)) {stop('M must be a numeric value.')}
  if (!is.numeric(N)) {stop('N must be a numeric array.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (!is.numeric(Fb)) {stop('Fb must be a numeric value.')}
  if (!is.numeric(E)) {stop('E must be a numeric array.')}
  if (!is.numeric(Catch)) {stop('Catch must be a numeric array.')}

  # acceptable values
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (ENM <= 0 || ENM > 3) {
    stop('ENM must be greater than 0 and less than or equal to 3.')}
  if (fdr <= 0) {stop('fdr must be greater than 0.')}
  if (sum(FM < 0) > 0) {
    stop('All values in FM must be greater than or equal to 0.')}
  if (M <= 0 || M > 1) {stop('M must be between 0 and 1.')}
  if (sum(N < 0) > 0) {stop('All values in N must be greater than or equal to 0.')}
  if (A <= 0) {stop('A must be greater than 0.')}
  if (Fb < 0) {stop('Fb must be greater than or equal to 0.')}
  if (sum(E < 0) > 0) {stop('All values in E must be greater than or equal to 0.')}
  if (sum(Catch < 0) > 0) {
    stop('All values in Catch must be greater than or equal to 0.')}

  # relational values
  if(dim(N)[1] != dim(Catch)[1] || dim(N)[1] != dim(FM)[1]) {
    stop('N, FM, or Catch has an incorrect number of age classes.')}
  if(dim(N)[2] != dim(E)[1] || dim(N)[2] != dim(Catch)[2] || dim(N)[2] != dim(FM)[2]) {
    stop('N, E, or Catch has an incorrect number of areas.')}
  if(dim(N)[3] != dim(E)[2] || dim(N)[3] != dim(Catch)[3] || dim(N)[3] != dim(FM)[3]) {
    stop('N, E, FM, or Catch has an incorrect number of time steps.')}
  if(dim(N)[4] != dim(E)[3] || dim(N)[4] != dim(Catch)[4] || dim(N)[4] != dim(FM)[4]) {
    stop('N, E, FM, or Catch has an incorrect number of control rules.')}
  if(dim(N)[5] != dim(E)[4] || dim(N)[5] != dim(Catch)[4] || dim(N)[5] != dim(FM)[5]) {
    stop('N, E, FM, or Catch has an incorrect number of final density ratios.')}
  if (t > dim(N)[3]) {stop('The given "t" value is too high for N.')}
  if (cr > dim(N)[4]) {stop('The given "cr" value is too high for N.')}
  if (ENM > dim(N)[5]) {stop('The given ENM value is too high for N.')}

  ##############################################################################

  # calculate the coefficient
  coeff <- FM[ , , t, cr, fdr]/(M + FM[ , , t, cr, fdr])

  # calculate catch at age
  Catch[ , , t, cr, fdr] <- coeff * N[ , , t, cr, fdr, ENM] * exp(-1*M - FM[ , , t, cr, fdr])

  return(Catch[, , t, cr, fdr])

}
