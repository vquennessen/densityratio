#' Fraction mature at age
#'
#' @param Rec_age numeric value, gives age at which fish has entered the fishery
#'    in years
#' @param Max_age numeric value, gives maximum age of fish or total lifespan in
#'    years
#' @param K_mat numeric value, slope of maturity curve
#' @param L numeric vector, gives length at age vector, in cm
#' @param L50 numeric value, gives length at 50% maturity, in cm
#'
#' @return a numeric vector of fractions mature at age, or the proportion of
#'    individuals that are expected to be mature at any age, on the interval
#'    (0, 1).
#' @export
#'
#' @examples
#' L <- length_at_age(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#' L2 = 47.95, K = 0.2022, All_ages = FALSE)
#' frac_mat(Rec_age = 3, Max_age = 35, K_mat = -0.4103, L, L50 = 39.53)

frac_mat = function(Rec_age, Max_age, K_mat, L, L50) {

  # Initialize fraction mature at age vector
  mature <- array(rep(0, Num), c(1, Num))

  # Calculate fraction mature at age
  mature <- (1)/(1 + exp(K_mat*(L - L50)))

  return(mature)

}
