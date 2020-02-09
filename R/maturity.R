#' Maturity at age
#'
#' @description \code{maturity} calculates the proportion of individuals at each
#'    age class, from age at recruitment to maximum age, that are expected to be
#'    mature.
#'
#' @param Rec_age numeric value, the age at recruitment, in years.
#' @param Max_age numeric value, the maximum age or total lifespan, in years.
#' @param K_mat numeric value, the slope of the maturity curve.
#' @param L numeric vector, the length at age vector, in cm.
#' @param L50 numeric value, the length at 50\% maturity, in cm.
#'
#' @return a numeric vector of fraction mature at ages, from age at recruitment
#'    to maximum age, on the interval (0, 1).
#' @export
#'
#' @examples
#' L <- length_age(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#'    L2 = 47.95, K = 0.2022, All_ages = FALSE)
#' maturity(Rec_age = 2, Max_age = 35, K_mat = -0.4103, L, L50 = 39.53)
maturity = function(Rec_age, Max_age, K_mat, L, L50) {

  ###### Error handling ########################################################

  # classes of variables
  if (!is.numeric(Rec_age)) {stop('Rec_age must be a numeric value.')}
  if (!is.numeric(Max_age)) {stop('Max_age must be a numeric value.')}
  if (!is.numeric(K_mat)) {stop('K_mat must be a numeric value.')}
  if (!is.numeric(L)) {stop('L must be a numeric vector.')}
  if (!is.numeric(L50)) {stop('L50 must be a numeric value.')}

  # acceptable values
  if (Rec_age <= 0) {stop('Rec_age must be greater than 0.')}
  if (K_mat >= 0) {stop('K_mat must be less than 0.')}
  if (sum(L <= 0) > 0) {stop('All values in L must be greater than 0.')}
  if (L50 <= 0) {stop('L50 must be greater than 0.')}

  ##############################################################################

  # relational values
  if (Rec_age >= Max_age) {stop('Rec_age must be less than Max_age.')}

  # number of age classes
  ages <- Rec_age:Max_age
  num <- length(ages)

  # Initialize fraction mature at age vector
  mature <- array(rep(0, num), c(1, num))

  # Calculate fraction mature at age
  mature <- (1)/(1 + exp(K_mat*(L - L50)))

  return(mature)

}
