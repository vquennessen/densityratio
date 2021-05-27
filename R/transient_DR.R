#' Transient DR
#'
#' \code{transient_DR} calculates the target density ratios, over all timesteps
#'    with management, assuming a transient control rule is being used.
#'
#' @param Time1 numeric value, the number of years to run the model before a
#'    marine reserve is implemented. Default value is 50.
#' @param TimeT numeric value, the number of years to run the model total.
#'    Default value is 70.
#' @param Final_DRs numeric vector, the final target density ratios.
#' @param Nat_mortality numeric vector, the estimates of natural mortality.
#' @param cr temporary numeric value, the current control rule.
#' @param fdr temporary numeric value, the current final target density ratio.
#'
#' @return a numeric vector, the target density ratios over all timesteps with
#'    transient control rules
#' @export
#'
#' @examples
#' transient_DR(Time1 = 50, TimeT = 70, Final_DRs = c(0.6, 0.8),
#'    Nat_mortality = c(0.17, 0.12, 0.22), cr = 1, fdr = 1)
transient_DR <- function(Time1 = 50, TimeT = 70, Final_DRs, Nat_mortality, cr,
                         fdr) {

  ###### Error handling ########################################################

  # classes of variables
  if (Time1 %% 1 != 0) {stop('Time1 must be an integer value.')}
  if (TimeT %% 1 != 0) {stop('TimeT must be an integer value.')}
  if (!is.numeric(Final_DRs)) {stop('Final_DRs must be a numeric vector.')}
  if (!is.numeric(Nat_mortality)) {stop('Nat_mortality must be a numeric vector.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (fdr %% 1 != 0) {stop('fdr must be an integer value.')}

  # acceptable values
  if (Time1 <= 0) {stop('Time1 must be greater than 0.')}
  if (TimeT <= 0) {stop('TimeT must be greater than 0.')}
  if (sum(Final_DRs <= 0) > 0) {
    stop('All values in Final_DRs must be greater than 0.')}
  if (sum(Nat_mortality <= 0) > 0 || sum(Nat_mortality > 1) > 0) {
    stop('All values in Nat_mortality must be between 0 and 1.')}
  if (cr <= 0 || cr > 3) {
    stop('cr must be greater than 0 and less than or equal to 3.')}
  if (fdr <= 0) {stop('fdr must be greater than 0.')}

  # relational values
  if (Time1 >= TimeT) {stop('TimeT must be greater than Time1.')}

  ##############################################################################

  # calculate target_DR based on transient timescales

  # set whether control rule uses correct (1, 2), low (3, 4), or high (5, 6)
  # estimate of natural mortality
  j <- cr / 2

  # set timesteps
  years <- 0:(TimeT - Time1)

  # calculate moving DR vector
  target_DR <- 1 - (1 - Final_DRs[fdr])*(1 - exp(-1 * Nat_mortality[j] * years))

  return(target_DR)

}
