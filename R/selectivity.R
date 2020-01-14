#' Selectivity at age
#'
#' \code{selectivity} returns a vector of selectivities at age, on the
#'    interval (0, 1).
#'
#' @param Fleets character vector, the names of all fishing fleets that
#'    contribute to selectivity at age on the interval (0, 1)
#' @param Max_age numeric value, gives maximum age of fish or total lifespan in
#'    years
#' @param Rec_age numeric value, gives age at which fish has entered the fishery
#'    in years
#' @param Alpha numeric vector, the alpha values for each fishing fleet in
#'    Fleets on the interval (0, 1), gives the relative slope of the upcurve
#' @param A50_up numeric vector, the age values at which selectivity is equal to
#'    0.5 for the upcurve each fishing fleet in Fleets on the interval (0, 1)
#' @param A50_down numeric vector, the age values at which selectivity is equal
#'    to 0.5 for the downcurve each fishing fleet in Fleets on the interval
#'    (0, 1)
#' @param F_fin numeric vector, the final selectivity values for fleets on the
#'    interval (0, 1)
#' @param Beta numeric vector, the beta values for each fishing fleet in Fleets
#'    on the interval (0, 1), gives the relative slope of the downcurve if
#'    dome-shaped
#' @param Cf numeric vector, the fraction of the whole fishery represented by
#'    each fleet
#' @param L numeric vector, gives length at age vector, in cm
#'
#' @return A vector of selectivities at age, from recruitment to maximum age, on
#'    the interval (0, 1).
#' @export
#'
#' @examples
#' L <- length_age(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#' L2 = 47.95, K = 0.2022, All_ages = FALSE)
#' selectivity(Rec_age = 2, Max_age = 35, L,
#' Fleets = c('sport', 'hook', 'trawl'), A50_up = c(2, 5, 10),
#' A50_down = c(6, 16, 35), Alpha = c(0.33, 0.6, 0.64),
#' F_fin = c(0.25, 0.06, 1), Beta = c(1.2, 0.6, 0), Cf = c(0.71, 0.28, 0.01))

selectivity <- function(Rec_age, Max_age, L, Fleets, A50_up, A50_down, Alpha,
                        F_fin, Beta, Cf) {

  # Calculated values
  ages <- Rec_age:Max_age
  num <- length(ages)

  # length of fleet vector
  f <- length(Fleets)
  l <- length(L)

  # translate A50_up and A50_down to lengths instead of ages
  L50up <- L[A50_up - Rec_age + 1]
  L50down <- L[A50_down - Rec_age + 1]

  # initialize upcurves and downcurves
  upcurve <- array(rep(NA, f*num), c(f, num))
  downcurve <- array(rep(NA, f*num), c(f, num))

  # initialize selectivity at age array
  # dimensions = age * fleet
  S <- array(rep(0, f*num), c(f, num))

  for (i in 1:f) {
    upcurve[i, ages - Rec_age + 1] <- 1 / (1 + exp(-1*Alpha[i]*(L - L50up[i])))

    downcurve[i, ages - Rec_age + 1] <- 1 -
      (1 - F_fin[i]) / (1 + exp(-1*Beta[i]*(L - L50down[i])))

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
