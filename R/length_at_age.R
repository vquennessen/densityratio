#' Length at age
#'
#' \code{length_at_age} returns a vector of lengths at ages, in centimeters.
#'
#' @param Rec_age numeric, gives age at which fish has entered the fishery in years
#' @param Max_age numeric, gives maximum age of fish or total lifespan in years
#' @param A1 numeric, age 1 of female fish in years
#' @param L1 numeric, gives the length of females at age 1 in cm
#' @param A2 numeric, age 2 of female fish in years
#' @param L2 numeric, gives the length of females at age 2 in cm
#' @param K numeric, gives the von Bertalanffy growth parameter for females
#' @param All_ages logical, determines if the returned vector should contain lengths for ages before recruitment
#' @return A vector of lengths at ages, from rec_age to max_age
#' @export
#'
#' @examples
#' length_at_age(Rec_age = 3, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#' L2 = 47.95, K = 0.2022, All_ages = FALSE)

length_at_age = function(Rec_age, Max_age, A1, L1, A2, L2, K, All_ages = F) {

  if (All_ages == T) {
    Ages <- 1:Max_age } else {
      Ages <- Rec_age:Max_age
    }

  L_inf <- L1 + (L2 - L1)/(1 - exp(-1*K*(A2 - A1)))
  lengths <- L_inf + (L1 - L_inf)*exp(-1*K*(Ages - A1))

  return(lengths)

}
