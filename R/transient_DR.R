#' Transient DR
#'
#' \code{transient_DR} calculates the target density ratios, over all timesteps
#'    with management, assuming a transient control rule is being used.
#'
#' @param Time1 numeric value, the number of years to run the model before a
#'    marine reserve is implemented. Default value is 50.
#' @param TimeT numeric value, the number of years to run the model total.
#'    Default value is 70.
#' @param Final_DR numeric value, the final target density ratio.
#' @param Nat_mortality numeric vector, the estimates of natural mortality.
#' @param nm temporary numeric value, the current natural mortality estimate.
#'
#' @return a numeric vector, the target density ratios over all timesteps with
#'    transient control rules
#' @export
#'
#' @examples
#' transient_DR(Time1 = 50, TimeT = 70, Final_DR = 0.6,
#'    Nat_mortality = c(0.09, 0.14, 0.19), nm = 1)
transient_DR <- function(Time1 = 50, TimeT = 70, Final_DR, Nat_mortality, nm) {

  # calculate target_DR based on transient timescales

  # set timesteps
  years <- 0:(TimeT - Time1)

  # calculate moving DR vector
  target_DR <- 1 - (1 - Final_DR)*(1 - exp(-1 * Nat_mortality[nm] * years))

  return(target_DR)

}
