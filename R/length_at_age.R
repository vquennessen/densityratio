#' Calculates length at age
#'
#' @param rec_age numeric value, gives age at which fish has entered the fishery
#' @param max_age numeric value
#' @param L1f numeric value
#' @param L2f numeric value
#' @param Kf numeric value
#' @param a1f numeric value
#' @param a2f numeric value
#' @param all_ages logical vector, if T, gives length at ages not included in
#' population dynamics
#'
#' @return A vector of lengths at ages
#' @export
#'
#' @examples
#' length_at_age(3, 40, 17, 49, -0.17, 1, 25, FALSE)

length_at_age = function(rec_age, max_age, L1f, L2f, Kf, a1f, a2f, all_ages = F) {

  if (all_ages == T) {
    ages <- 1:max_age } else {
      ages <- rec_age:max_age
    }

  L_inf <- L1f + (L2f - L1f)/(1 - exp(-1*Kf*(a2f - a1f)))
  L <- L_inf + (L1f - L_inf)*exp(-1*Kf*(ages - a1f))

  return(L)

}
