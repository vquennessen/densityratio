#' Length at age
#'
#' \code{length_age} returns a vector of lengths at ages, in centimeters.
#'
#' @param Rec_age numeric value, the age at recruitment, in years.
#' @param Max_age numeric value, the maximum age or total lifespan, in years.
#' @param A1 numeric value, age 1 of female fish in years.
#' @param L1 numeric value, the length of females at age 1, in cm.
#' @param A2 numeric value, age 2 of female fish in years.
#' @param L2 numeric value, the length of females at age 2, in cm.
#' @param K numeric value, the von Bertalanffy growth parameter for females.
#' @param All_ages logical value, should the returned vector contain lengths for
#'    ages before recruitment.
#'
#' @return A numeric vector of lengths at ages, from age at recruitment to
#'    maximum age.
#' @export
#'
#' @examples
#' length_age(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#'    L2 = 47.95, K = 0.2022, All_ages = FALSE)

length_age = function(Rec_age, Max_age, A1, L1, A2, L2, K, All_ages = F) {

  ###### Error handling ########################################################

  # classes of variables
  if (Rec_age %% 1 != 0) {stop('Rec_age must be an integer value.')}
  if (Max_age %% 1 != 0) {stop('Max_age must be an integer value.')}
  if (A1 %% 1 != 0) {stop('A1 must be an integer value.')}
  if (!is.numeric(L1)) {stop('L1 must be a numeric value.')}
  if (A2 %% 1 != 0) {stop('A2 must be an integer value.')}
  if (!is.numeric(L2)) {stop('L2 must be a numeric value.')}
  if (!is.numeric(K)) {stop('K must be a numeric value.')}
  if (!is.logical(All_ages)) {stop('All_ages must be a logical value.')}

  # acceptable values
  if (Rec_age <= 0) {stop('Rec_age must be greater than 0.')}
  if (A1 <= 0) {stop('A1 must be greater than 0.')}
  if (L1 <= 0) {stop('L1 must be greater than 0.')}
  if (K <= 0) {stop('K must be greater than 0.')}

  # relational values
  if (Rec_age >= Max_age) {stop('Rec_age must be less than Max_age.')}
  if (A1 >= A2) {stop('A1 must be less than A2.')}
  if (L1 >= L2) {stop('L1 must be less than L2.')}

  ##############################################################################

  # should the lengths vector include ages before recruitment?
  if (All_ages == T) {
    ages <- 1:Max_age } else {
      ages <- Rec_age:Max_age
    }

  # calculate maximum length given von Bertalanffy equation
  L_inf <- L1 + (L2 - L1)/(1 - exp(-1*K*(A2 - A1)))

  # calculate length at age (cm)
  lengths <- L_inf + (L1 - L_inf)*exp(-1*K*(ages - A1))

  return(lengths)

}
