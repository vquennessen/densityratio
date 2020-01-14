#' Weight at age
#'
#' \code{weight} returns a vector of weights at age, in kg.
#'
#' @param L numeric vector, gives length at age vector, in cm
#' @param WA numeric value, coefficient in weight at length equation, where
#'    weight = WA * length ^ WB
#' @param WB numeric value, exponent in weight at length equation, where
#'    weight = WA * length ^ WB
#'
#' @return a numeric vector of weights at ages, from age at recruitment to
#'    maximum age
#' @export
#'
#' @examples
#' L <- length_age(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#' L2 = 47.95, K = 0.2022, All_ages = FALSE)
#' weight(L, WA = 1.68e-5, WB = 3)
weight = function(L, WA, WB) {

  # Weight at age
  # Based on Babcock & MacCall (2011): Eq. (11)
  weights <- WA*L^WB

  return(weights)

}
