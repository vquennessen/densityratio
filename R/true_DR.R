#' True density ratio
#'
#' \code{true_DR} returns the true density ratio for the current timestep,
#'    updating the Density_ratio array.
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param fdr temporary numeric value, the current final target density ratio.
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
#' @param ENM numeric value, the nm value that represents the 'true' population.
#'
#' @return a numeric array that updates the true density ratios for the most
#'    current timestep and control rule.
#' @export
#'
#' @examples
#' A = 5; TimeT = 70; CR = 6; NM = 3; FDR = 4; Time2= 20
#' Abundance_all <- array(rep(3400, A*TimeT*CR*NM*FDR), c(A, TimeT, CR, NM, FDR))
#' Density_ratio <- array(rep(0, (Time2 + 1)*CR*FDR), c(Time2 + 1, CR, FDR))
#' true_DR(t = 51, cr = 1, fdr = 1, Abundance_all, Inside = 3,
#'    Outside = c(1, 2, 4, 5), Density_ratio, Time1 = 50, ENM = 2)
true_DR <- function(t, cr, fdr, Abundance_all, Inside = 3,
                    Outside = c(1, 2, 4, 5), Density_ratio, Time1 = 50, ENM) {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (fdr %% 1 != 0) {stop('fdr must be an integer value.')}
  if (!is.numeric(Abundance_all)) {
    stop('Abundance_all must be a numeric array.')}
  if (sum(Inside %% 1 != 0) != 0) {stop('Inside must be a vector of integers.')}
  if (sum(Outside %% 1 != 0) != 0) {stop('Outside must be a vector of integers.')}
  if (!is.numeric(Density_ratio)) {
    stop('Density_ratio must be a numeric array.')}
  if (Time1 %% 1 != 0) {stop('Time1 must be an integer value.')}

  # acceptable values
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (fdr <= 0) {stop('fdr must be greater than 0.')}
  if (sum(Abundance_all < 0) > 0) {
    stop('All values in Abundance_all must be greater than or equal to 0.')}
  if (sum(Inside < 0) > 0) {
    stop('All values in Inside must be greater than or equal to 0.')}
  if (sum(Outside < 0) > 0) {
    stop('All values in Outside must be greater than or equal to 0.')}
  if (sum(Density_ratio < 0) > 0) {
    stop('All values in Density_ratio must be greater than or equal to 0.')}
  if (Time1 <= 0) {stop('Time1 must be greater than 0.')}

  # relational values
  if (sum(intersect(Inside, Outside)) > 0) {
    stop('Areas cannot both be inside and outside the marine reserve.')}
  if (dim(Abundance_all)[1] != length(Inside) + length(Outside)) {
    stop('Abundance_all has the wrong number of areas.')}
 if (t > dim(Abundance_all)[2]) {
   stop('Abundance_all has the wrong number of time steps.')}
  if (cr > dim(Abundance_all)[3]) {
    stop('Abundance_all has the wrong number of control rules.')}
  if (fdr> dim(Abundance_all)[5]) {
    stop('Abundance_all has the wrong number of final density ratios.')}

  ##############################################################################

  # Density of fish outside marine reserve(s)
  Outside_density <- sum(Abundance_all[Outside, t, cr, ENM, fdr]) / length(Outside)

  # Density of fish inside marine reserve(s)
  Inside_density <- sum(Abundance_all[Inside, t, cr, ENM, fdr]) / length(Inside)

  # True density ratio
  Density_ratio[t - Time1 + 1, cr, fdr] <- Outside_density / Inside_density

  return(Density_ratio[t - Time1 + 1, cr, fdr])
}
