#' Recruitment
#'
#' \code{recruitment} calculates the number of new recruits entering the
#'    population
#'
#' @param a temporary numeric value, the current area.
#' @param t temporary numeric value, the current time step .
#' @param cr temporary numeric value, the current control rule .
#' @param nm temporary numeric value, the current natural mortality estimate.
#' @param SSB numeric array, the spawning stock biomass of the whole stock for
#'    each area, at each timestep, under each control rule, and for each
#'    estimate of natural mortality, in kg.
#' @param A numeric value, the number of total areas in the model. Default
#'    value is 5.
#' @param R0 numeric value, set arbitrarily, the unfished recruitment. Default
#'    value is 1e+5.
#' @param H numeric value, the steepness of the stock-recruitment curve.
#' @param B0 numeric value, the unfished biomass, in kg.
#' @param Eps numeric matrix, the recruitment error terms.
#' @param Sigma_R numeric value, the recruitment standard deviation.
#' @param Rec_age numeric value, the age at recruitment, in years.
#' @param Recruitment_mode character value, values can be:
#'    'closed' - the recruits in each area originate from adults in that area.
#'    'pool' - the recruits in each area come from a pool of larvae produced by
#'       adults in all areas.
#'    Default value is 'pool'.
#'
#' @return a numeric value representing the number of new recruits coming into
#'    the population in area a, at timestep t, under control rule cr, with an
#'    estimate of natural mortality of nm.
#' @export
#'
#' @examples
#' SSB <- array(rep(10, 5*70*6*3), c(5, 70, 6, 3))
#' NuR <- array(rnorm(5*70*6*3, 0, 0.5), c(5, 70, 6, 3))
#' Eps <- epsilon(A = 5, TimeT = 70, CR = 6, NM = 3, NuR, Rho_R = 0)
#' recruitment(a = 1, t = 3, cr = 1, nm = 2, SSB, A = 5, R0 = 1e+5, H = 0.65,
#'    B0 = 1e+5/1.1, Eps, Sigma_R = 0.5, Rec_age = 2, Recruitment_mode = 'pool')
recruitment = function(a, t, cr, nm, SSB, A = 5, R0 = 1e+5, H, B0, Eps, Sigma_R,
                       Rec_age, Recruitment_mode) {

  ###### Error handling ########################################################

  # classes of variables
  if (a %% 1 != 0) {stop('a must be an integer value.')}
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (nm %% 1 != 0) {stop('nm must be an integer value.')}
  if (!is.numeric(SSB)) {stop('SSB must be a numeric array.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (R0 %% 1 != 0) {stop('R0 must be an integer value.')}
  if (!is.numeric(H)) {stop('H must be a numeric value.')}
  if (!is.numeric(B0)) {stop('B0 must be a numeric value.')}
  if (!is.numeric(Eps)) {stop('Eps must be a numeric array.')}
  if (!is.numeric(Sigma_R)) {stop('Sigma_R must be a numeric array.')}
  if (Rec_age %% 1 != 0) {stop('Rec_age must be an integer value.')}
  if (!is.character(Recruitment_mode)) {
    stop('Recruitment mode must be a character value.')}

  # acceptable values
  if (a <= 0) {stop('a must be greater than 0.')}
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (nm <= 0 || nm > 3) {
    stop('nm must be greater than 0 and less than or equal to 3.')}
  if (sum(SSB < 0) > 0) {
    stop('All values in SSB must be greater than or equal to 0.')}
  if (A <= 0) {stop('A must be greater than 0.')}
  if (R0 <= 0) {stop('R0 must be greater than 0.')}
  if (H <= 0 || H > 1) {stop('H must be between 0 and 1.')}
  if (B0 <= 0) {stop('B0 must be greater than 0.')}
  if (Sigma_R <= 0) {stop('Sigma_R must be greater than 0.')}
  if (Rec_age <= 0) {stop('Rec_age must be greater than 0.')}
  if (Recruitment_mode != 'pool' && Recruitment_mode != 'closed') {
    stop('Recruitment_mode must be either "pool" or "closed".')}

  # relational values
  if(dim(SSB)[1] != dim(Eps)[1]) {
    stop('SSB or Eps has an incorrect number of areas.')}
  if(dim(SSB)[2] != dim(Eps)[2]) {
    stop('SSB or Eps has an incorrect number of time steps.')}
  if(dim(SSB)[3] != dim(Eps)[3]) {
    stop('SSB or Eps has an incorrect number of control rules.')}
  if(dim(SSB)[4] != dim(Eps)[4]) {
    stop('SSB or Eps has an incorrect number of natural mortality estimates.')}
  if (a > dim(SSB)[1]) {stop('The given "a" value is too high for SSB.')}
  if (t > dim(SSB)[2]) {stop('The given "t" value is too high for SSB.')}
  if (cr > dim(SSB)[3]) {stop('The given "cr" value is too high for SSB.')}
  if (nm > dim(SSB)[4]) {stop('The given "nm" value is too high for SSB.')}
  if (a > A) {stop('a must be less than or equal to A.')}

  ##############################################################################

  # Recruitment
  # Based on Babcock & MacCall (2011): Eq. (3)
  # Dimensions = 1 * 1

  # closed recruitment indicates that recruits are only offspring of adults from
  # the same area
  if (Recruitment_mode == 'closed') { ssb <- SSB[a, t - Rec_age, cr, nm]}

  # pool recruitment indicates that recruits come from a larval pool produced by
  # adults from all areas
  if (Recruitment_mode == 'pool') { ssb <- sum(SSB[, t - Rec_age, cr, nm] / A)}

  adjR0 <- R0 / A
  adjB0 <- B0 / A

  R1 <- (0.8 * adjR0 * H * ssb) / (0.2 * adjB0 * (1 - H) + (H - 0.2) * ssb)
  recruits <- R1 * (exp(Eps[a, t, cr, nm] - Sigma_R^2 / 2))

  return(recruits)

}
