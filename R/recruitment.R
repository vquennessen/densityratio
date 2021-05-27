#' Recruitment
#'
#' \code{recruitment} calculates the number of new recruits entering the
#'    population
#'
#' @param t temporary numeric value, the current time step .
#' @param cr temporary numeric value, the current control rule .
#' @param fdr temporary numeric value, the current final target density ratio.
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
#'    'regional_DD' - larvae experience regional density dependence before
#'       settling evenly across all areas
#'    'local_DD' - larvae experience local density dependence before settling
#'       evely across all areas
#'    Default value is 'pool'.
#' @param LDP numeric value, the larval drift proportion, the proportion of
#'    larvae that drift from one area to an adjacent area before settling.
#'    Default value is 0.1.
#'
#' @return a numeric value representing the number of new recruits coming into
#'    the population in area a, at timestep t, under control rule cr.
#' @export
#'
#' @examples
#' A = 5; TimeT = 70; CR = 6; FDR = 4
#' SSB <- array(rep(10, A*TimeT*CR*FDR), c(A, TimeT, CR, FDR))
#' NuR <- array(rnorm(A*TimeT*CR*FDR, 0, 0.5), c(A, TimeT, CR, FDR))
#' Eps <- epsilon(A, TimeT, CR, FDR, NuR, Rho_R = 0)
#' recruitment(t = 3, cr = 1, fdr = 1, SSB, A, R0 = 1e+5, H = 0.65,
#'    B0 = 1e+5/1.1, Eps, Sigma_R = 0.5, Rec_age = 2, Recruitment_mode = 'pool',
#'    LDP = 0.1)
recruitment = function(t, cr, fdr, SSB, A = 5, R0 = 1e+5, H, B0, Eps,
                       Sigma_R, Rec_age, Recruitment_mode, LDP = 0.1) {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (fdr %% 1 != 0) {stop('fdr must be an integer value.')}
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
  if (!is.numeric(LDP)) {stop('LDP must be a numeric value.')}

  # acceptable values
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (fdr <= 0) {stop('fdr must be greater than 0.')}
  if (sum(SSB < 0) > 0) {
    stop('All values in SSB must be greater than or equal to 0.')}
  if (A <= 0) {stop('A must be greater than 0.')}
  if (R0 <= 0) {stop('R0 must be greater than 0.')}
  if (H <= 0 || H > 1) {stop('H must be between 0 and 1.')}
  if (B0 <= 0) {stop('B0 must be greater than 0.')}
  if (Sigma_R <= 0) {stop('Sigma_R must be greater than 0.')}
  if (Rec_age <= 0) {stop('Rec_age must be greater than 0.')}
  if (Recruitment_mode != 'pool' && Recruitment_mode != 'closed' &&
      Recruitment_mode != 'regional_DD' && Recruitment_mode != 'local_DD') {
    stop('Recruitment_mode must be either "pool", "closed", "regional_DD", or
         "local_DD".')}
  if (LDP < 0) {stop('LDP must be greater than or equal to 0.')}

  # relational values
  if(dim(SSB)[1] != dim(Eps)[1]) {
    stop('SSB or Eps has an incorrect number of areas.')}
  if(dim(SSB)[2] != dim(Eps)[2]) {
    stop('SSB or Eps has an incorrect number of time steps.')}
  if(dim(SSB)[3] != dim(Eps)[3]) {
    stop('SSB or Eps has an incorrect number of control rules.')}
  if(dim(SSB)[4] != dim(Eps)[4]) {
    stop('SSB or Eps has an incorrect number of final density ratios.')}
  if (t > dim(SSB)[2]) {stop('The given "t" value is too high for SSB.')}
  if (cr > dim(SSB)[3]) {stop('The given "cr" value is too high for SSB.')}
  if (fdr > dim(SSB)[4]) {stop('The given "fdr" value is too high for SSB.')}
  ##############################################################################

  # adjust R0 and B0 per area
  adjR0 <- R0 / A
  adjB0 <- B0 / A

    recruits <- rep(0, A)
    ssb <- SSB[, t - Rec_age, cr, fdr]
    sum_ssb <- sum(ssb)

    # closed recruitment indicates that recruits are only offspring of adults
    # from the same area
      if (Recruitment_mode == 'closed') {

        R1 <- (0.8 * adjR0 * H * ssb) / (0.2 * adjB0 * (1 - H) + (H - 0.2) * ssb)

        # 'pool' recruitment indicates that recruits come from a larval pool
        # produced by adults from all areas
      } else if (Recruitment_mode == 'pool') {

        nume <- 0.8 * adjR0 * H * sum_ssb / A
        denom <- 0.2 * adjB0 * (1 - H) + (H - 0.2) * ssb

        R1 <- nume / denom

        # regional / stock larval density dependence and recruits distributed
        # evenly across areas; OR larvae distributed evenly across areas then
        # local density dependence in each area
      } else if (Recruitment_mode == 'regional_DD') {

        nume <- 0.8 * R0 * H * sum_ssb
        denom <- A * (0.2 * B0 * (1 - H) + (H - 0.2) * sum_ssb)

        R1 <- rep(nume / denom, each = A)

        # larval density dependence within areas and recruitment in equal
        # amounts in each area
      } else if (Recruitment_mode == 'local_DD') {

        nume <- 0.8 * adjR0 * H * ssb
        denom <- 0.2 * adjB0 * (1 - H) + (H - 0.2) * ssb

        R1 <- 1/A * sum(nume / denom)

      }

      recruits <- R1 * (exp(Eps[, t, cr, fdr] - Sigma_R^2 / 2))

      # larval movement if there are multiple areas and LDP != 0
      if (A > 1 && LDP > 0) {

        # First area to second area
        recruits[1] <- (1 - LDP)*recruits[1] + LDP*recruits[2]

        # Intermediate areas to adjacent areas
        for (a in 2:(A-1)) {
          recruits[a] <- (1 - 2*LDP)*recruits[a] +
            LDP*(recruits[a - 1] + recruits[a + 1])
        }

        # Last area to next to last area
        recruits[A] <- (1 - LDP)*recruits[A] + LDP*recruits[A - 1]

      }

  return(recruits)

}
