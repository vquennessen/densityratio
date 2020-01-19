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

  # number of age classes
  ages <- Rec_age:Max_age
  num <- length(ages)

  # Initialize fraction mature at age vector
  mature <- array(rep(0, num), c(1, num))

  # Calculate fraction mature at age
  mature <- (1)/(1 + exp(K_mat*(L - L50)))

  return(mature)

}
