#' Effort allocation
#'
#' \code{effort_allocation} redistributes fishing effort for the next timestep.
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param NM  numeric value, the total number of natural mortality estimates.
#' @param fdr temporary numeric value, the current final target density ratio.
#' @param Allocation character value, how effort is to be allocated. Values can
#'    be:
#'    'IFD' - the ideal free distribution. Effort is allocated proportional to
#'       the yield caught in each area in the previous timestep.
#'    'equal' - Effort is allocated equally between all areas.
#'    Default value is 'IFD'.
#' @param E numeric array, the relative fishing effort displayed in each area,
#'    at each time step, under each control rule, and for each natural mortality
#'    estimate.
#' @param Yield numeric array, the yield caught in each area, at each timestep,
#'    under each control rule, and for each estimate of natural mortality, in kg.
#' @param Time1 numeric value, the number of years to run the model before a
#'    marine reserve is implemented. Default value is 50.
#' @param Inside numeric vector, the area(s) inside the marine reserve. Default
#'    value is c(3).
#' @param Outside numeric vector, the area(s) outside the marine reserve.
#'    Default value is c(1, 2, 4, 5).
#'
#' @return a numeric array updated with the fishing effort allocation for the
#'    next timestep.
#' @export
#'
#' @examples
#' A = 5; TimeT = 70; CR = 6; NM = 3; FDR = 4
#' E <- array(rep(0.2, A*TimeT*CR*NM*FDR), c(A, TimeT, CR, NM, FDR))
#' Yield <- array(rep(2458, A*TimeT*CR*NM*FDR), c(A, TimeT, CR, NM, FDR))
#' effort_allocation(t = 51, cr = 1, NM, fdr = 1, Allocation = 'IFD', E,
#'    Yield, Time1 = 50, Inside = 3, Outside = c(1, 2, 4, 5))
effort_allocation <- function(t, cr, NM, fdr, Allocation = 'IFD', E, Yield,
                              Time1 = 50, Inside = c(3),
                              Outside = c(1, 2, 4, 5)) {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (NM %% 1 != 0) {stop('NM must be an integer value.')}
  if (fdr %% 1 != 0) {stop('fdr must be an integer value.')}
  if (!is.character(Allocation)) {stop('Allocation must be a character value.')}
  if (!is.numeric(E)) {stop('E must be a numeric array.')}
  if (!is.numeric(Yield)) {stop('Yield must be a numeric array.')}
  if (Time1 %% 1 != 0) {stop('Time1 must be an integer value.')}
  if (sum(Inside %% 1 != 0) != 0) {stop('Inside must be a vector of integers.')}
  if (sum(Outside %% 1 != 0) != 0) {stop('Outside must be a vector of integers.')}

  # acceptable values
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (NM <= 0 || NM > 3) {
    stop('NM must be greater than 0 and less than or equal to 3.')}
  if (fdr <= 0) {stop('fdr must be greater than 0.')}
  if (Allocation != 'IFD' && Allocation != 'equal') {
    stop('Allocation must be either "IFD" or "equal".')}
  if (sum(E < 0) > 0) {stop('All values in E must be greater than or equal to 0.')}
  if (sum(Yield < 0) > 0) {
    stop('All values in Yield must be greater than or equal to 0.')}
  if (Time1 <= 0) {stop('Time1 must be greater than 0.')}
  if (sum(Inside < 0) > 0) {
    stop('All values in Inside must be greater than or equal to 0.')}
  if (sum(Outside < 0) > 0) {
    stop('All values in Outside must be greater than or equal to 0.')}

  # relational values
  if(dim(E)[1] != dim(Yield)[1]) {
    stop('E or Yield has an incorrect number of areas.')}
  if(dim(E)[2] != dim(Yield)[2]) {
    stop('E or Yield has an incorrect number of time steps.')}
  if(dim(E)[3] != dim(Yield)[3]) {
    stop('E or Yield has an incorrect number of control rules.')}
  if(dim(E)[4] != dim(Yield)[4]) {
    stop('E or Yield has an incorrect number of natural mortality estimates.')}
  if(dim(E)[5] != dim(Yield)[5]) {
    stop('E or Yield has an incorrect number of final target density ratios.')}
  if (t > dim(E)[2]) {stop('The given "t" value is too high for E.')}
  if (cr > dim(E)[3]) {stop('The given "cr" value is too high for E.')}
  if (NM > dim(E)[4]) {stop('The given "nm" value is too high for E.')}
  if (fdr > dim(E)[5]) {stop('The given "fdr" value is too high for E.')}
  if (sum(intersect(Inside, Outside)) > 0) {
    stop('Areas cannot both be inside and outside the marine reserve.')}

  ##############################################################################

  if (NM > 1) {

  # number of areas not in a reserve
  outs <- length(Outside)
  ins <- length(Inside)
  all <- outs + ins

  if (t == Time1) {

    E[Outside, t, cr, , fdr] <- array(rep(rep(colSums(E[, t - 1, cr, , fdr]) / outs,
                                              outs), NM), c(outs, NM))
    E[Inside, t, cr, , fdr] <- 0

    # If effort is allocated using the ideal free distribution, effort for one
    # year depends on the distribution of yield from the previous year
  } else if (t != Time1 & Allocation == 'IFD') {

    yield <- Yield[, t - 1, cr, , fdr]
    prop_yield <- sweep(yield, MARGIN = 2, STATS = colSums(yield), FUN = '/')
    colSums_E <- array(rep(colSums(E[, t, cr, , fdr]), each = all), c(all, NM))
    E[, t, cr, , fdr] <- colSums_E * prop_yield

    # Otherwise, distribute effort equally between the four areas outside the
    # marine reserve, regardless of yield
  } else if (t < Time1 & Allocation == 'equal') {

    E[, t, cr, , fdr] <- array(rep(colSums(E[, t, cr, , fdr])/all,
                                   each = all), c(all, NM))

  } else if (t > Time1 & Allocation == 'equal') {

    E[Outside, t, cr, , fdr] <- array(rep(colSums(E[, t, cr, , fdr])/outs,
                                          each = outs), c(outs, NM))
    E[Inside, t, cr, , fdr] <- 0

  }

  } else if (NM == 1) {

    # number of areas not in a reserve
    outs <- length(Outside)
    ins <- length(Inside)
    all <- outs + ins

    if (t == Time1) {

      E[Outside, t, cr, 1, fdr] <- rep(colSums(E[, t - 1, cr, 1, fdr]) / outs, outs)
      E[Inside, t, cr, 1, fdr] <- 0

      # If effort is allocated using the ideal free distribution, effort for one
      # year depends on the distribution of yield from the previous year
    } else if (t != Time1 & Allocation == 'IFD') {

      prop_yield <- Yield[, t - 1, cr, 1, fdr] / sum(Yield[, t - 1, cr, 1, fdr])
      E[, t, cr, 1, fdr] <- sum(E[, t, cr, 1, fdr]) * prop_yield

      # Otherwise, distribute effort equally between the four areas outside the
      # marine reserve, regardless of yield
    } else if (t < Time1 & Allocation == 'equal') {

      E[, t, cr, 1, fdr] <- rep(sum(E[, t, cr, 1, fdr])/all, each = all)

    } else if (t > Time1 & Allocation == 'equal') {

      E[Outside, t, cr, 1, fdr] <- rep(colSums(E[, t, cr, , fdr])/outs, each = outs)
      E[Inside, t, cr, 1, fdr] <- 0

    }

}

  return(E[, t, cr, , fdr])

}
