#' Density ratio
#'
#' \code{density_ratio} estimates the density ratio for a stock in line with a
#'    specific control rule.
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param fdr temporary numeric value, the current final target density ratio.
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
#'    'far' - sample only the first and last areas (assuming the reserve is
#'       somewhere in the middle.
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
#' A = 5; TimeT = 70; CR = 6; FDR = 4; Transects = 24
#' Count <- array(rep(5, A*TimeT*Transects*2*CR*FDR),
#'    c(A, TimeT, Transects, 2, CR, FDR))
#' density_ratio(t = 51, cr = 1, fdr = 1, A = 5, Count, Years_sampled = 1,
#'    Areas_sampled = 'all', Ind_sampled = 'all', Transects = 24, Inside = 3,
#'    Outside = c(1, 2, 4, 5))
density_ratio <- function (t, cr, fdr, A, Count, Years_sampled = 1,
                           Areas_sampled = 'all', Ind_sampled = 'all',
                           Transects = 24, Inside, Outside) {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (fdr %% 1 != 0) {stop('fdr must be an integer value.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (!is.numeric(Count)) {stop('Count must be a numeric array.')}
  if (Years_sampled %% 1 != 0 && !is.null(Years_sampled)) {
    stop('Years_sampled must be an integer value or NULL.')}
  if (!is.character(Areas_sampled) && !is.null(Areas_sampled)) {
    stop('Areas_sampled must be a character value or NULL.')}
  if (!is.character(Ind_sampled) && !is.null(Ind_sampled)) {
    stop('Ind_sampled must be a character value or NULL.')}
  if (Transects %% 1 != 0) {stop('Transects must be an integer value.')}
  if (sum(Inside %% 1 != 0) != 0) {stop('Inside must be a vector of integers.')}
  if (sum(Outside %% 1 != 0) != 0) {stop('Outside must be a vector of integers.')}

  # acceptable values
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (fdr <= 0) {stop('fdr must be greater than 0.')}
  if (A <= 0) {stop('A must be greater than 0.')}
  if (sum(Count < 0) > 0) {
    stop('All values in Count must be greater than or equal to 0.')}
  if (Years_sampled <= 0 && !is.null(Years_sampled)) {
    stop('Years_sampled must be greater than 0 or NULL.')}
  if (is.numeric(Years_sampled) && Years_sampled <= 0) {
    stop('Years_sampled must be greater than 0 or NULL.')}
  if (is.character(Areas_sampled) && Areas_sampled != 'far' &&
      Areas_sampled != 'all' ) {
    stop('Areas_sampled must be either "far" or "all" or NULL.')}
  if (is.character(Ind_sampled) && Ind_sampled != 'mature' &&
      Ind_sampled != 'all') {
    stop('Ind_sampled must be either "mature" or "all" or NULL.')}
  if (Transects <= 0) {stop('Transects must be greater than 0.')}
  if (sum(Inside < 0) > 0) {
    stop('All values in Inside must be greater than or equal to 0.')}
  if (sum(Outside < 0) > 0) {
    stop('All values in Outside must be greater than or equal to 0.')}

  # relational values
  if (sum(Inside > A) > 0) {
    stop('All values in Inside must be less than or equal to A.')}
  if (sum(Outside > A) > 0) {
    stop('All values in Outside must be less than or equal to A.')}
  if (sum(intersect(Inside, Outside)) > 0) {
    stop('Areas cannot both be inside and outside the marine reserve.')}
  if(dim(Count)[1] != A) {stop('Count has an incorrect number of areas.')}
  if(dim(Count)[3] != Transects) {stop('Count has the wrong number of transects.')}
  if (t > dim(Count)[2]) {stop('The given "t" value is too high for Count.')}
  if (cr > dim(Count)[3]) {stop('The given "cr" value is too high for Count.')}
  if (fdr > dim(Count)[6]) {stop('The given "fdr" value is too high for Count.')}

  ##############################################################################

  # sample all fish or just mature fish
  if (Ind_sampled == 'all') { ind <- 1 } else { ind <- 2 }

  # time sampled = 1 or more years
  if (Years_sampled == 1) { years <- t - 1
  } else { years <- (t - Years_sampled):(t - 1) }

  # calculate count inside marine reserve
  count_in <- Count[Inside, years, , ind, cr, fdr]

  # calculate counts outside marine reserve
  if (Areas_sampled == 'far') {
    count_out <- Count[c(1, A), years, , ind, cr, fdr]
  } else if (Areas_sampled == 'all') {
    count_out <- Count[Outside, years, , ind, cr, fdr]
  }

  # calculate density ratio
  DR <- (sum(count_out)/(length(count_out)))/(sum(count_in)/length(count_in))

  return(DR)

}
