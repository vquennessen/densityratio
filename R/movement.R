#' Movement
#'
#' \code{movement} redistributes individuals amongst areas.
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param fdr temporary numeric value, the current final target density ratio.
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
#' @export
#'
#' @examples
#' n = 34; A = 5; TimeT = 70; CR = 6; FDR = 4
#' N <- array(rep(10, n*A*TimeT*CR*FDR), c(n, A, TimeT, CR, FDR))
#' movement(t = 1, cr = 1, fdr = 1, N, A = 5, AMP = 0.1)
movement <- function(t, cr, fdr, N, A, AMP = 0.1) {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (fdr %% 1 != 0) {stop('fdr must be an integer value.')}
  if (!is.numeric(N)) {stop('N must be a numeric array.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (!is.numeric(AMP)) {stop('AMP must be a numeric value.')}

  # acceptable values
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (fdr <= 0) {stop('fdr must be greater than 0.')}
  if (sum(N < 0) > 0) {
    stop('All values in N must be greater than or equal to 0.')}
  if (A <= 0) {stop('A must be greater than 0.')}
  if (AMP < 0) {stop('AMP must be greater than or equal to 0.')}

  # relational values
  if (A != dim(N)[2]) {stop('N has the wrong number of areas.')}
  if (t > dim(N)[3]) {stop('The given "t" value is too high for N.')}
  if (cr > dim(N)[4]) {stop('The given "cr" value is too high for N.')}
  if (fdr > dim(N)[5]) {stop('The given "fdr" value is too high for N.')}

  ##############################################################################

  # First area to second area
  N[, 1, t, cr, fdr] <- (1 - AMP)*N[, 1, t, cr, fdr] + AMP*N[, 2, t, cr, fdr]

  # Intermediate areas to adjacent areas
  for (a in 2:(A-1)) {
    N[, a, t, cr, fdr] <- (1 - 2*AMP)*N[, a, t, cr, fdr] +
      AMP*(N[, a - 1, t, cr, fdr] + AMP*N[, a + 1, t, cr, fdr])
  }

  # Last area to next to last area
  N[, A, t, cr, fdr] <- (1 - AMP)*N[, A, t, cr, fdr] + AMP*N[, A - 1, t, cr, fdr]

  return(N)

}
