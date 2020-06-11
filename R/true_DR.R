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
#' @param BM logical value, are the control rules from Babcock and MacCall 2011?
#'    Default value is FALSE.
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
#' @param A numeric value, the number of total areas in the model. Default value
#'    is 5.
#'
#' @return a numeric array that updates the true density ratios for the most
#'    current timestep and control rule.
#' @export
#'
#' @examples
#' A = 5; TimeT = 70; CR = 6; NM = 2; FDR = 4; Time2 = 20
#' Abundance <- array(rep(3400, A*TimeT*CR*NM*FDR*1),
#'    c(A, TimeT, CR, NM, FDR, 1))
#' Density_ratio <- array(rep(0, TimeT*CR*FDR), c(TimeT, CR, FDR))
#' true_DR(t = 51, cr = 1, fdr = 1, Abundance, Inside = 3,
#'    Outside = c(1, 2, 4, 5), Density_ratio, BM = FALSE, Years_sampled = 1,
#'    Areas_sampled = 'all', Ind_sampled = 'all', A = 5)
true_DR <- function(t, cr, fdr, Abundance, Inside = 3, Outside = c(1, 2, 4, 5),
                    Density_ratio, BM = FALSE, Years_sampled = 1,
                    Areas_sampled = 'all', Ind_sampled = 'all', A = 5) {

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
  if (!is.logical(BM)) {stop('BM must be a logical value.')}
  if (Years_sampled %% 1 != 0 && !is.null(Years_sampled)) {
    stop('Years_sampled must be an integer value or NULL.')}
  if (!is.character(Areas_sampled) && !is.null(Areas_sampled)) {
    stop('Areas_sampled must be a character value or NULL.')}
  if (!is.character(Ind_sampled) && !is.null(Ind_sampled)) {
    stop('Ind_sampled must be a character value or NULL.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}

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
  if (A <= 0) {stop('A must be greater than 0.')}

  # relational values
  if (sum(intersect(Inside, Outside)) > 0) {
    stop('Areas cannot both be inside and outside the marine reserve.')}
  if (dim(Abundance)[1] != length(Inside) + length(Outside)) {
    stop('Abundance has the wrong number of areas.')}
  if (t > dim(Abundance)[2]) {
    stop('Abundance has the wrong number of time steps.')}
  if (cr > dim(Abundance)[3]) {
    stop('Abundance has the wrong number of control rules.')}
  if (fdr > dim(Abundance)[5]) {
    stop('Abundance has the wrong number of final density ratios.')}

  ##############################################################################

  if (BM == TRUE) {

    Years_sampled <- ifelse(cr == 2, 3, 1)
    Areas_sampled <- ifelse(cr == 4, 'far', 'all')
    Ind_sampled <- ifelse(cr == 5, 'mature', 'all')

    # sample all fish or just mature fish
    if (Ind_sampled == 'all') { ind <- 1 } else { ind <- 2 }

    # time sampled = 1 or more years
    if (Years_sampled == 1) { years <- t - 1
    } else { years <- (t - Years_sampled):(t - 1) }

    # calculate density inside marine reserve
    Inside_count <- sum(Abundance[Inside, years, , ind, cr, 1, fdr])

    # calculate counts outside marine reserve
    if (Areas_sampled == 'far') {
      Outside_count <- sum(Abundance[c(1, A), years, , ind, cr, 1, fdr])
    } else if (Areas_sampled == 'all') {
      Outside_count <- sum(Abundance[Outside, years, , ind, cr, 1, fdr])
    }

    # Density of fish outside marine reserve(s)
    Outside_density <- Outside_count / length(Outside)

    # Density of fish inside marine reserve(s)
    Inside_density <- Inside_count / length(Inside)

    # True density ratio
    Density_ratio[t, cr, fdr] <- Outside_density / Inside_density

  } else {

    # Density of fish outside marine reserve(s)
    Outside_density <- sum(Abundance[Outside, t, cr, 1, fdr, 1]) / length(Outside)

    # Density of fish inside marine reserve(s)
    Inside_density <- sum(Abundance[Inside, t, cr, 1, fdr, 1]) / length(Inside)

    # True density ratio
    Density_ratio[t, cr, fdr] <- Outside_density / Inside_density

  }



  return(Density_ratio[t, cr, fdr])
}
