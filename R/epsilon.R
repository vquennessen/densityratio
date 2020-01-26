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
#' @param NM numeric value, the total number of estimated values of natural
#'    mortality. Default value is 3.
#' @param NuR numeric vector, the recruitment random normal variable, pulled
#'    from a normal distribution of mean 0 and standard deviation equal to
#'    Sigma_R.
#' @param Rho_R numeric value, the recruitment autocorrelation on the interval
#'    (-1, 1). Default value is 0.
#'
#' @return a numeric vector of recruitment error terms, of dimensions
#'    A \* timeT \* CR \* NM.
#'
#' @examples
#' NuR <- array(stats::rnorm(5*70*6*3, 0, 0.5), c(5, 70, 6, 3))
#' epsilon(A = 5, TimeT = 70, CR = 6, NM = 3, NuR, Rho_R = 0)
epsilon <- function (A = 5, TimeT = 70, CR = 6, NM = 3, NuR, Rho_R = 0) {

  # initialize epsilon vector
  # Dimensions = area * timeT * CR * M values (3)
  Eps <- array(rep(0, A*TimeT*CR*NM), c(A, TimeT, CR, NM))

  # eps[, 1, ]
  Eps[, 1, , ] <- NuR[, 1, , ]*sqrt(1 + Rho_R^2)

  # fill in rest of epsilon vector
  for (a in 1:A) {
    for (t in 2:TimeT) {
      for (cr in 1:CR) {
        for (nm in 1:NM) {
        Eps[a, t, cr, nm] <- Rho_R*Eps[a, t-1, cr, nm] + NuR[a, t, cr, nm]*sqrt(1 + Rho_R^2)
        }
      }
    }
  }

  return(Eps)

}
