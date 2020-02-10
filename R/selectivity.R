#' Selectivity at age
#'
#' \code{selectivity} returns a vector of selectivities at age, on the
#'    interval (0, 1).
#'
#' @param Fleets character vector, the names of all fishing fleets that
#'    contribute to selectivity at age.
#' @param Max_age numeric value, the maximum age or total lifespan, in years.
#' @param Rec_age numeric value, the age at recruitment, in years.
#' @param Alpha numeric vector, the alpha values for each fishing fleet in
#'    Fleets, the relative slope of the upcurve on the interval (0, 1).
#' @param A50_up numeric vector, the age values at which selectivity is equal to
#'    0.5 on the upcurve for each fishing fleet in Fleets on the interval (0, 1).
#' @param A50_down numeric vector, the age values at which selectivity is equal
#'    to 0.5 on the downcurvefor  each fishing fleet in Fleets on the interval
#'    (0, 1).
#' @param F_fin numeric vector, the final selectivity values for each fishing
#'    fleet in Fleets on the interval (0, 1).
#' @param Beta numeric vector, the beta values for each fishing fleet in Fleets,
#'    the relative slope of the downcurve if dome-shaped on the interval (0, 1).
#' @param Cf numeric vector, the fraction of the whole fishery represented by
#' @param A1 numeric value, age 1 of female fish in years.
#' @param L1 numeric value, the length of females at age 1, in cm.
#' @param A2 numeric value, age 2 of female fish in years.
#' @param L2 numeric value, the length of females at age 2, in cm.
#' @param K numeric value, the von Bertalanffy growth parameter for females.
#'    each fleet, on the interval (0, 1).
#' @return A numeric vector of selectivities at age, from recruitment to maximum
#'    age, on the interval (0, 1).
#' @export
#'
#' @examples
#' selectivity(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#'    L2 = 47.95, K = 0.2022, Fleets = c('sport', 'hook', 'trawl'),
#'    A50_up = c(2, 5, 10), A50_down = c(6, 16, 35), Alpha = c(0.33, 0.6, 0.64),
#'    F_fin = c(0.25, 0.06, 1), Beta = c(1.2, 0.6, 0), Cf = c(0.71, 0.28, 0.01))

selectivity <- function(Rec_age, Max_age, A1, L1, A2, L2, K, Fleets, A50_up,
                        A50_down, Alpha, F_fin, Beta, Cf) {

  ###### Error handling ########################################################

  # classes of variables
  if (Rec_age %% 1 != 0) {stop('Rec_age must be an integer value.')}
  if (Max_age %% 1 != 0) {stop('Max_age must be an integer value.')}
  if (A1 %% 1 != 0) {stop('A1 must be an integer value.')}
  if (!is.numeric(L1)) {stop('L1 must be a numeric value.')}
  if (A2 %% 1 != 0) {stop('A2 must be an integer value.')}
  if (!is.numeric(L2)) {stop('L2 must be a numeric value.')}
  if (!is.numeric(K)) {stop('K must be a numeric value.')}
  if (!is.character(Fleets)) {stop('Fleets must be a character vector.')}
  if (A50_up %% 1 != 0) {stop('A50_up must be an integer value.')}
  if (A50_down %% 1 != 0) {stop('A50_down must be an integer value.')}
  if (!is.numeric(Alpha)) {stop('Alpha must be a numeric vector.')}
  if (!is.numeric(F_fin)) {stop('F_fin must be a numeric vector.')}
  if (!is.numeric(Beta)) {stop('Beta must be a numeric vector.')}
  if (!is.numeric(Cf)) {stop('Cf must be a numeric vector.')}

  # acceptable values
  if (Rec_age <= 0) {stop('Rec_age must be greater than 0.')}
  if (A1 <= 0) {stop('A1 must be greater than 0.')}
  if (L1 <= 0) {stop('L1 must be greater than 0.')}
  if (K <= 0) {stop('K must be greater than 0.')}
  if (sum(A50_up <= 0) > 0) {stop('All values in A50_up must be greater than 0.')}
  if (sum(A50_down < 0) > 0) {
    stop('All values in A50_down must be greater than 0.')}
  if (sum(Alpha < 0) > 0) {
    stop('All values in Alpha must be greater than 0.')}
  if (sum(F_fin < 0) > 0) {
    stop('All values in F_fin must be greater than or equal to 0.')}
  if (sum(Beta < 0) > 0) {
    stop('All values in Beta must be greater than or equal to 0.')}
  if (sum(Cf <= 0) > 0) {stop('All values in Cf must be greater than 0.')}

  # relational values
  if (Rec_age >= Max_age) {stop('Rec_age must be less than Max_age.')}
  if (A1 >= A2) {stop('A1 must be less than A2.')}
  if (L1 >= L2) {stop('L1 must be less than L2.')}
  if (length(A50_up) != length(Fleets)) {
    stop('A50_up must have the same number of elements as Fleets.')}
  if (length(A50_down) != length(Fleets)) {
    stop('A50_down must have the same number of elements as Fleets.')}
  if (length(Alpha) != length(Fleets)) {
    stop('Alpha must have the same number of elements as Fleets.')}
  if (length(F_fin) != length(Fleets)) {
    stop('F_fin must have the same number of elements as Fleets.')}
  if (length(Beta) != length(Fleets)) {
    stop('Beta must have the same number of elements as Fleets')}
  if (length(Cf) != length(Fleets)) {
    stop('Cf must have the same number of elements as Fleets.')}
  ##############################################################################

  # length at age for all ages
  L <- length_age(Rec_age, Max_age, A1, L1, A2, L2, K, All_ages = T)

  # Calculated values
  ages <- Rec_age:Max_age
  num <- length(ages)

  # length of fleet vector
  f <- length(Fleets)
  l <- length(L)

  # translate A50_up and A50_down to lengths instead of ages
  L50up <- L[A50_up]
  L50down <- L[A50_down]

  # initialize upcurves and downcurves
  upcurve <- array(rep(NA, f*num), c(f, num))
  downcurve <- array(rep(NA, f*num), c(f, num))

  # initialize selectivity at age array
  # dimensions = age * fleet
  S <- array(rep(0, f*num), c(f, num))

  for (i in 1:f) {
    upcurve[i, ages - Rec_age + 1] <- 1 / (1 + exp(-1*Alpha[i]*(L[Rec_age:Max_age] - L50up[i])))

    downcurve[i, ages - Rec_age + 1] <- 1 -
      (1 - F_fin[i]) / (1 + exp(-1*Beta[i]*(L[Rec_age:Max_age] - L50down[i])))

    # if (Beta[i] == 0) { downcurve[i, ages + 1] <- rep(1, num + 1) }

    for (a in 1:num) {
      S[i, a] <- min(upcurve[i, a], downcurve[i, a])
    }

    S[i, ] <- Cf[i]*S[i, ]

  }

  #### plot selectivities to double check they're right #####

  # colors <- rainbow(n = f, s = 1, v = 1, start = 0, end = max(1, f - 1)/f,
  # alpha = 1)
  #
  # plot(ages, S[1, ], type = 'l', lwd = 2, col = colors[1],
  #      ylim = c(0, 1),
  #      main = "Vic's Attempt",
  #      xlab = 'Age (year)',
  #      ylab = 'Selectivity',
  #      xlim = c(0, 40))
  #
  # for (i in 2:f) {
  # lines(ages, S[i, ], type = 'l', lwd = 2, col = colors[i])
  # }
  #
  # legend(x = 'topright', Fleets, lwd = 2, cex = 0.8, col = colors)

  selectivity <- colSums(S)

  return (selectivity)

}
