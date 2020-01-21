#' True density ratio
#'
#' \code{true_DR} returns the true density ratio for the current timestep,
#'    updating the Density_ratio array.
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param Abundance_all numeric array, the total number of individuals in each
#'    area, at each timestep, under all control rules, with all estimates of
#'    natural mortality.
#' @param Inside numeric vector, the area(s) inside the marine reserve. Default
#'    is c(3).
#' @param Outside numeric vector, the area(s) outside the marine reserve.
#'    Default is c(1, 2, 4, 5).
#' @param Density_ratio numeric array, the true density ratios at each timestep,
#'    under each control rule, and for each estimate of natural mortality.
#' @param Time1 numeric value, the number of years to run the model before a
#'    marine reserve is implemented. Default value is 50.
#'
#' @return a numeric array that updates the true density ratios for the most
#'    current timestep and control rule.
#' @export
#'
#' @examples
#' Abundance_all <- array(rep(3400, 5*70*6*3), c(5, 70, 6, 3))
#' Density_ratio <- array(rep(0, (20 + 1)*6), c(20 + 1, 6))
#' true_DR(t = 51, cr = 1, Abundance_all, Inside = c(3), Outside = c(1, 2, 4, 5),
#'    Density_ratio, Time1 = 50)
true_DR <- function(t, cr, Abundance_all, Inside = c(3), Outside = c(1, 2, 4, 5),
                    Density_ratio, Time1 = 50) {

  # Density of fish outside marine reserve(s)
  Outside_density <- sum(Abundance_all[Outside, t, cr, 2]) / length(Outside)

  # Density of fish inside marine reserve(s)
  Inside_density <- sum(Abundance_all[Inside, t, cr, 2]) / length(Inside)

  # True density ratio
  Density_ratio[t - Time1 + 1, cr] <- Outside_density / Inside_density

  return(Density_ratio)
}
