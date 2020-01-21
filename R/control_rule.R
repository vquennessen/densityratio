#' Control rule
#'
#' \code{control_rule} determines how to calculate the density ratio and manage
#'    the fishery (i.e. set fishing effort for the next timestep).
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param A numeric value, the number of total areas in the model. Default
#'    value is 5.
#' @param E numeric array, the relative fishing effort displayed in each area,
#'    at each time step, under each control rule, and for each natural mortality
#'    estimate.
#' @param Count numeric array, the number of individuals estimated to be in each
#'    area, at each timestep, under each control rule, for each estimate of
#'    natural mortality, for both all individuals and just mature individuals.
#' @param Time1 numeric value, the number of years to run the model before a
#'    marine reserve is implemented. Default value is 50.
#' @param TimeT numeric value, the number of years to run the model total.
#'    Default value is 70.
#' @param Transects numerical value, the number of sampling transects conducted
#'    in each area to estimate density ratio. Default value is 24.
#' @param Nat_mortality numeric vector, the estimates of natural mortality.
#' @param Final_DR numeric value, the final target density ratio.
#' @param Inside numeric vector, the area(s) inside the marine reserve. Default
#'    is c(3).
#' @param Outside numeric vector, the area(s) outside the marine reserve.
#'    Default is c(1, 2, 4, 5).
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
#' @param Years_sampled numeric value, the number of years of sampling upon
#'    which to base the estimate of density ratio. Default value is 1.
#'
#' @return a numeric vector of fishing effort for the next timestep, under the
#'    specific control rule, with a specific estimate of natural mortality.
#' @export
#'
#' @examples
#' E <- array(rep(1, 5*70*6*3), c(5, 70, 6, 3))
#' Count <- array(rep(5, 5*70*24*2*6*3), c(5, 70, 24, 2, 6, 3))
#' control_rule(t = 51, cr = 1, A = 5, E, Count, Time1 = 50, TimeT = 70,
#'    Transects = 24, Nat_mortality = c(0.09, 0.14, 0.19), Final_DR = 0.6,
#'    Inside = c(3), Outside = c(1, 2, 4, 5), Areas_sampled = 'all',
#'    Ind_sampled = 'all', Years_sampled = 1)
control_rule <- function(t, cr, A = 5, E, Count, Time1 = 50, TimeT = 70,
                         Transects = 24, Nat_mortality, Final_DR, Inside = c(3),
                         Outside = c(1, 2, 4, 5), Areas_sampled = 'all',
                         Ind_sampled = 'all', Years_sampled = 1) {


  if (cr < 4) { nm = cr } else { nm = cr - 3 }

  # static control rules, with constant target density ratios
  if (cr < 4) {

    DR <- density_ratio(t, cr, nm, A, Count, Years_sampled, Areas_sampled,
                        Ind_sampled, Transects, Inside, Outside)

    # calculate effort at the next timestep
    E[, t + 1, cr, ] <- management(t, cr, E, DR, target_DR = Final_DR,
                                   floor_DR = 0.2, effort_inc_allowed = 0.1,
                                   Time1)
  }

  # transient control rules with shifting target density ratios
  else {

    target_DR <- transient_DR(Time1, TimeT, Final_DR, Nat_mortality, nm)

    # calculate density ratio
    DR <- density_ratio(t, cr, nm, A, Count, Years_sampled, Areas_sampled,
                        Ind_sampled, Transects, Inside, Outside)

    # calculate effort at the next timestep
    E[, t + 1, cr, ] <- management(t, cr, E, DR, target_DR[t - Time1 + 1],
                                   floor_DR = 0.2, effort_inc_allowed = 0.1,
                                   Time1)
  }

  return(E[, t + 1, cr, ])

}
