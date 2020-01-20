#' Density ratio
#'
#' \code{density_ratio} estimates the density ratio for a stock in line with a
#'    specific control rule.
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param nm temporary numeric value, the current natural mortality estimate.
#' @param A numeric value, the number of total areas in the model. Default
#'    value is 5.
#' @param Count numeric array, the number of individuals estimated to be in each
#'    area, at each timestep, under each control rule, for each estimate of
#'    natural mortality, for both all individuals and just mature individuals.
#' @param Years_sampled numeric value, the number of years of sampling upon
#'    which to base the estimate of density ratio. Default value is 1.
#' @param Areas_sampled character value, the areas to be sampled to calculate
#'    density ratio. Values can be:
#'    'all' - sample all areas.
#'    'far' - sample only areas farthest from the marine reserve.
#'    Default value is 'all'.
#' @param Ind_sampled character value, the individuals to be sampled to
#'    calculate density ratio. Values can be:
#'    'all' - sample all individuals.
#'    'mature' - sample only mature individuals.
#'    Default value is 'all'.
#' @param Transects numeric value, the number of sampling transects conducted
#'    in each area to estimate density ratio. Default value is 24.
#' @param Inside numeric vector, the area(s) inside the marine reserve. Default
#'    is c(3).
#' @param Outside numeric vector, the area(s) outside the marine reserve.
#'    Default is c(1, 2, 4, 5).
#'
#' @return numeric value, the calculated density ratio for a specific timestep,
#'    under a specific control rule, with a specific estimate of natural
#'    mortality.
#' @export
#'
#' @examples
#' Count <- array(rep(5, 5*70*24*2*6*3), c(5, 70, 24, 2, 6, 3))
#' density_ratio(t = 2, cr = 1, nm = 1, A = 5, Count, Years_sampled = 1,
#'    Areas_sampled = 'all', Ind_sampled = 'all', Transects = 24,
#'    Inside = c(3), Outside = c(1, 2, 4, 5))
density_ratio <- function (t, cr, nm, A, Count, Years_sampled = 1,
                           Areas_sampled = 'all', Ind_sampled = 'all',
                           Transects = 24, Inside, Outside) {

  # sample all fish or just mature fish
  if (Ind_sampled == 'all') {
    ind <- 1
  } else if (Ind_sampled == 'mature') {
    ind <- 2
  }

  # calculate count inside marine reserve
  count_in <- Count[Inside, t - 1, , ind, cr, nm]

  # time sampled = 1 or 3 years
  if (Years_sampled == 1) {
    years <- t - 1
  } else {
    years <- (t - Years_sampled):(t - 1)
  }

  # calculate counts outside marine reserve
  if (Areas_sampled == 'far') {
    count_out <- Count[c(1, A), years, , ind, cr, nm]
  } else if (Areas_sampled == 'all') {
    count_out <- Count[Outside, years, , ind, cr, nm]
  }

  # calculate density ratio
  DR <- (sum(count_out)/(Transects*nrow(count_out)))/(sum(count_in)/Transects)

  return(DR)

}
