#' Epsilon
#'
#' \code{epsilon} returns a vector of error terms for recruitment, based on
#' equation 4 in Babcock & MacCall (2011).
#'
#' @param A numeric value, the number of total areas in the model. Default value
#'    is 5.
#' @param TimeT numeric value, the number of years to run the model total.
#'    Default value is 70.
#' @param CR numeric value, the number of control rules to be compared. Default
#'    value is 6.
#' @param FDR numeric value, the total number of final density ratios.
#' @param NuR numeric vector, the recruitment random normal variable, pulled
#'    from a normal distribution of mean 0 and standard deviation equal to
#'    Sigma_R.
#' @param Rho_R numeric value, the recruitment autocorrelation on the interval
#'    (-1, 1). Default value is 0.
#'
#' @return a numeric vector of recruitment error terms, of dimensions
#'    A \* timeT \* CR.
#' @export
#'
#' @examples
#' A = 5; TimeT = 70; CR = 6; FDR = 4
#' NuR <- array(stats::rnorm(A*TimeT*CR*FDR, 0, 0.5), c(A, TimeT, CR, FDR))
#' epsilon(A, TimeT, CR, FDR, NuR, Rho_R = 0)
epsilon <- function (A = 5, TimeT = 70, CR = 6, FDR, NuR, Rho_R = 0) {

  ###### Error handling ########################################################

  # classes of variables
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (TimeT %% 1 != 0) {stop('TimeT must be an integer value.')}
  if (CR %% 1 != 0) {stop('CR must be an integer value.')}
  if (FDR %% 1 != 0) {stop('FDR must be an integer value.')}
  if (!is.numeric(NuR)) {stop('NuR must be a numeric array.')}
  if (!is.numeric(Rho_R)) {stop('Rho_R must be a numeric array.')}

  # acceptable values
  if (A <= 0) {stop('A must be greater than 0.')}
  if (TimeT <= 0) {stop('TimeT must be greater than 0.')}
  if (CR < 1) {stop('CR must be greater than or equal to 1.')}
  if (FDR <= 0) {stop('FDR must be greater than 0.')}
  if (Rho_R < -1 || Rho_R > 1) {stop('Rho_R must be between -1 and 1.')}

  # relational values
  if(dim(NuR)[1] != A) {stop('NuR has an incorrect number of areas.')}
  if(dim(NuR)[2] != TimeT) {stop('NuR has an incorrect number of time steps.')}
  if(dim(NuR)[3] != CR) {stop('NuR has an incorrect number of control rules.')}
  if(dim(NuR)[4] != FDR) {
    stop('NuR has an incorrect number of final density ratios.')}

  ##############################################################################

  # initialize epsilon vector
  # Dimensions = area * timeT * CR * M values (3)
  Eps <- array(rep(0, A*TimeT*CR*FDR), c(A, TimeT, CR, FDR))

  # initial value
  Eps[, 1, , ] <- NuR[, 1, , ]*sqrt(1 + Rho_R^2)

  # fill in rest of epsilon vector
  for (a in 1:A) {
    for (t in 2:TimeT) {
      for (cr in 1:CR) {
        for (fdr in 1:FDR) {
          Eps[a, t, cr, fdr] <- Rho_R*Eps[a, t - 1, cr, fdr] +
            NuR[a, t, cr, fdr]*sqrt(1 + Rho_R^2)
        }
      }
    }
  }


  return(Eps)

}
