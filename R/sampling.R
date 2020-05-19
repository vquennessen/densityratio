#' Sampling
#'
#' \code{sampling} samples from the current population to update the Count
#'    array.
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param nm temporary numeric value, the current natural mortality estimate.
#' @param fdr temporary numeric value, the current final target density ratio.
#' @param Delta numeric value, the proportion of positive transects divided by
#'    depletion, also known as the constant of proportionality
#' @param Gamma numeric value, the average value of a positive transect divided
#'    by depletion
#' @param Abundance numeric array, the total number of all and/or mature
#'    individuals in each area, at each timestep, under all control rules, with
#'    all estimates of natural mortality.
#' @param Transects numerical value, the number of sampling transects conducted
#'    in each area to estimate density ratio. Default value is 24.
#' @param X numeric value, the average value of individuals seen during positive
#'    transects.
#' @param Count numeric array, the number of individuals estimated to be in each
#'    area, at each timestep, under each control rule, for each estimate of
#'    natural mortality, for both all individuals and just mature individuals.
#' @param NuS numeric array, the sampling normal variable pulled from a normal
#'    distribution with mean 0 and sd Sigma_S.
#' @param A numeric value, the number of total areas in the model. Default value
#'    is 5.
#' @param Ind_sampled character value, the individuals to be sampled to
#'    calculate density ratio. Values can be:
#'    'all' - sample all individuals.
#'    'mature' - sample only mature individuals.
#'    Default value is 'all'.
#'
#' @return a two-dimentional numeric array with updated observed counts for the
#'    fished species in this particular area, at this particular timestep,
#'    under this particular control rule, with this particular estimate of
#'    natural mortality.
#' @export
#'
#' @importFrom stats rbinom
#'
#' @examples
#' n = 34; A = 5; TimeT = 70; CR = 6; NM = 3; FDR = 4; Transects = 24
#' Abundance <- array(rep(340, A*TimeT*CR*NM*FDR*1),
#'    c(A, TimeT, CR, NM, FDR, 1))
#' Count <- array(rep(50, A*TimeT*Transects*2*CR*NM*FDR),
#'    c(A, TimeT, Transects, 2, CR, NM, FDR))
#' NuS <- array(rnorm(A*TimeT*CR*NM*FDR*1, 0, 0.89), c(A, TimeT, CR, NM, FDR, 1))
#' sampling(t = 51, cr = 1, nm = 1, fdr = 1, Delta = 1.6, Gamma = 31.6,
#'    Abundance, Transects, X = 15.42, Count, NuS, A, Ind_sampled = 'all')
sampling <- function(t, cr, nm, fdr, Delta, Gamma, Abundance, Transects = 24,
                     X, Count, NuS, A = 5, Ind_sampled = 'all') {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (nm %% 1 != 0) {stop('nm must be an integer value.')}
  if (fdr %% 1 != 0) {stop('fdr must be an integer value.')}
  if (!is.numeric(Delta)) {stop('Delta must be a numeric value.')}
  if (!is.numeric(Gamma)) {stop('Gamma must be a numeric value.')}
  if (!is.numeric(Abundance)) {stop('Abundance must be a numeric array.')}
  if (Transects %% 1 != 0) {stop('Transects must be an integer value.')}
  if (!is.numeric(X)) {stop('X must be a numeric value.')}
  if (!is.numeric(Count)) {stop('Count must be a numeric array.')}
  if (!is.numeric(NuS)) {stop('NuS must be a numeric array.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (!is.character(Ind_sampled) && !is.null(Ind_sampled)) {
    stop('Ind_sampled must be a character value or NULL.')}

  # acceptable values
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (nm <= 0 || nm > 3) {
    stop('nm must be greater than 0 and less than or equal to 3.')}
  if (fdr <= 0) {stop('fdr must be greater than 0.')}
  if (Delta <= 0) {stop('Delta must be greater than 0.')}
  if (Gamma <= 0) {stop('Gamma must be greater than 0.')}
  if (sum(Abundance < 0) > 0) {
    stop('All values in Abundance must be greater than or equal to 0.')}
  if (Transects <= 0) {stop('Transects must be greater than 0.')}
  if (X < 0) {stop('X must be greater than or equal to 0.')}
  if (sum(Count < 0) > 0) {
    stop('All values in Count must be greater than or equal to 0.')}
  if (A <= 0) {stop('A must be greater than 0.')}
  if (is.character(Ind_sampled) && Ind_sampled != 'mature' &&
      Ind_sampled != 'all') {
    stop('Ind_sampled must be either "mature" or "all" or NULL.')}

  # relational values
  if (dim(Abundance)[1] != A || dim(Count)[1] != A) {
    stop('Abundance or Count has an incorrect number of areas.')}
  if (t > dim(Abundance)[2] || t > dim(Count)[2]) {
    stop('The given "t" value is too high for Abundance or Count.')}
  if (cr > dim(Abundance)[3]|| cr > dim(Count)[5]) {
    stop('The given "cr" value is too high for Abundance or Count.')}
  if (nm > dim(Abundance)[4] || nm > dim(Count)[6]) {
    stop('The given "nm" value is too high for Abundance or Count.')}
  if (fdr > dim(Abundance)[5] || fdr > dim(Count)[7]) {
    stop('The given "fdr" value is too high for Abundance or Count.')}

  ##############################################################################

  # Calculate probability of detection based on odds ratio
  # Based on Babcock & MacCall (2011): Eq. (12)
  A_all <- Abundance[, t - 1, cr, nm, fdr, 1]
  total_all <- sum(A_all)
  odds_all <-  (Delta * A_all) / (total_all / A)
  p_all <- 1 / (1 + exp(odds_all))

  # Determine if species is detected at least once, and replace a random 0 with
  #   a 1 if all zeros to prevent errors in calculating density ratio
  # Dimensions = 1 * transects
  presence_all <- array(rbinom(Transects, 1, p_all), c(Transects, 1))
  if (sum(presence_all) == 0) {r <- sample(1:Transects, 1); presence_all[r] = 1}

  # Calculate species count given transects with positive visuals
  nus <- NuS[, t - 1, cr, nm, fdr, 1]
  All <- Gamma*Abundance[, t - 1, cr, nm, fdr, 1]*exp(nus)
  Count[, t, , 1, cr, nm, fdr] <- presence_all %*% All

  if (Ind_sampled == 'mature' || is.null(Ind_sampled)) {

    # Calculate probability of detection based on odds ratio
    A_mature <- Abundance[, t - 1, cr, nm, fdr, 2]
    total_mature <- sum(A_mature)
    odds_mature <-  (Delta * A_mature) / (total_mature / A)
    p_mature <- 1 / (1 + exp(odds_mature))

    # Determine if species is detected at least once, and replace a random 0
    #   with a 1 if all zeros to prevent errors in calculating density ratio
    presence_mature <- array(rbinom(Transects, 1, p_mature), c(Transects, 1))
    if (sum(presence_mature) == 0) {
      r <- sample(1:Transects, 1); presence_mature[r] = 1}

    # Calculate species count given transects with positive visuals
    nus <- NuS[, t - 1, cr, nm, fdr, 2]
    Mature <- Gamma*A_mature*exp(nus)
    Count[, t, , 2, cr, nm, fdr] <- presence_mature %*% Mature

  }

  return(Count[, t, , , cr, nm, fdr])

}
