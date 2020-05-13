#' True density ratio
#'
#' \code{true_DR} returns the true density ratio for the current timestep,
#'    updating the Density_ratio array.
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param fdr temporary numeric value, the current final target density ratio.
#' @param Abundance numeric array, the total number of individuals in each area,
#'    at each timestep, under all control rules, with all estimates of natural
#'    mortality.
#' @param Inside numeric vector, the area(s) inside the marine reserve. Default
#'    value is 3.
#' @param Outside numeric vector, the area(s) outside the marine reserve.
#'    Default value is c(1, 2, 4, 5).
#' @param Density_ratio numeric array, the true density ratios at each timestep,
#'    under each control rule, and for each estimate of natural mortality.
#' @param Ind_sampled character value, the individuals to be sampled to
#'    calculate density ratio. Values can be:
#'    'all' - sample all individuals.
#'    'mature' - sample only mature individuals.
#'    Default value is 'all'.
#'
#' @return a numeric array that updates the true density ratios for the most
#'    current timestep and control rule.
#' @export
#'
#' @examples
#' A = 5; TimeT = 70; CR = 2; FDR = 4; Time2 = 20
#' Abundance <- array(rep(3400, A*TimeT*CR*FDR*1), c(A, TimeT, CR, FDR, 1))
#' Density_ratio <- array(rep(0, TimeT*CR*FDR), c(TimeT, CR, FDR))
#' true_DR(t = 51, cr = 1, fdr = 1, Abundance, Inside = 3,
#'    Outside = c(1, 2, 4, 5), Density_ratio, Ind_sampled = 'all')
true_DR <- function(t, cr, fdr, Abundance, Inside = 3, Outside = c(1, 2, 4, 5),
                    Density_ratio, Ind_sampled) {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (fdr %% 1 != 0) {stop('fdr must be an integer value.')}
  if (!is.numeric(Abundance)) {stop('Abundance must be a numeric array.')}
  if (sum(Inside %% 1 != 0) != 0) {stop('Inside must be a vector of integers.')}
  if (sum(Outside %% 1 != 0) != 0) {stop('Outside must be a vector of integers.')}
  if (!is.numeric(Density_ratio)) {
    stop('Density_ratio must be a numeric array.')}
  if (!is.character(Ind_sampled) && !is.null(Ind_sampled)) {
    stop('Ind_sampled must be a character value or NULL.')}

  # acceptable values
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (fdr <= 0) {stop('fdr must be greater than 0.')}
  if (sum(Abundance < 0) > 0) {
    stop('All values in Abundance must be greater than or equal to 0.')}
  if (sum(Inside < 0) > 0) {
    stop('All values in Inside must be greater than or equal to 0.')}
  if (sum(Outside < 0) > 0) {
    stop('All values in Outside must be greater than or equal to 0.')}
  if (sum(Density_ratio < 0) > 0) {
    stop('All values in Density_ratio must be greater than or equal to 0.')}

  # relational values
  if (sum(intersect(Inside, Outside)) > 0) {
    stop('Areas cannot both be inside and outside the marine reserve.')}
  if (dim(Abundance)[1] != length(Inside) + length(Outside)) {
    stop('Abundance has the wrong number of areas.')}
  if (t > dim(Abundance)[2]) {
    stop('Abundance has the wrong number of time steps.')}
  if (cr > dim(Abundance)[3]) {
    stop('Abundance has the wrong number of control rules.')}
  if (fdr > dim(Abundance)[4]) {
    stop('Abundance has the wrong number of final density ratios.')}

  ##############################################################################

  i <- ifelse(is.null(Ind_sampled), 2, ifelse(Ind_sampled == 'all', 1, 2))

  # Density of fish outside marine reserve(s)
  Outside_density <- sum(Abundance[Outside, t, cr, fdr, i]) / length(Outside)

  # Density of fish inside marine reserve(s)
  Inside_density <- sum(Abundance[Inside, t, cr, fdr, i]) / length(Inside)

  # True density ratio
  Density_ratio[t, cr, fdr] <- Outside_density / Inside_density

  return(Density_ratio[t, cr, fdr])

}
