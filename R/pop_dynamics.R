#' Population Dynamics
#'
#' \code{pop_dynamics} returns updated values for fishing mortality, numbers at
#'    age, abundance of all and mature individuals, biomass, and spawning stock
#'    biomass in a specific area, at a specific time step, under a specific
#'    control rule, and with a specific estimate of natural mortality
#'
#' @param a temporary numeric value, the current area.
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param nm temporary numeric value, the current natural mortality estimate.
#' @param Rec_age numeric value, the age at recruitment, in years.
#' @param Max_age numeric value, the maximum age or total lifespan, in years.
#' @param SSB numeric array,  the spawning stock biomass of the whole stock for
#'    each area, at each timestep, under each control rule, and for each
#'    estimate of natural mortality.
#' @param N numeric array, the number of individuals at each age, in each
#'    area, at each timestep, under each control rule, and for each estimate of
#'    natural mortality.
#' @param W numeric vector, the estimated weight at age from age at recruitment
#'    to maximum age, in kg.
#' @param Mat numeric vector, the estimated fraction of individuals mature at
#'    each age, from age at recruitment to maximum age, on the interval (0, 1).
#' @param A numeric value, the number of total areas in the model. Default
#'    value is 5.
#' @param R0 numeric value, set arbitrarily, the unfished recruitment. Default
#'    value is 1e+5.
#' @param H numeric value, the steepness of the stock-recruitment curve.
#' @param B0 numeric value, the unfished biomass.
#' @param Eps numeric array, the recruitment error terms.
#' @param Sigma_R numeric value, the recruitment standard deviation.
#' @param Fb numeric value, the historical fishing effort for the fished species.
#' @param E numeric array, the relative fishing effort displayed in each area,
#'    at each time step, under each control rule, and for each natural mortality
#'    estimate.
#' @param S numeric vector, the selectivities at age from age at recruitment to
#'    maximum age, on the interval (0, 1).
#' @param NM numeric value, the total number of estimated values of natural
#'    mortality. Default value is 3.
#' @param FM numeric array that corresponds to the fishing mortality at each
#'    age in each area, at each timestep, under all control rules, with all
#'    estimates of natural mortality.
#' @param A50_mat numeric value, the first age at which 50\% or more individuals
#'    are estimated to be mature. on the interval (Rec_age, Max_age).
#' @param Abundance_all numeric array, the total number of individuals in each
#'    area, at each timestep, under all control rules, with all estimates of
#'    natural mortality.
#' @param Abundance_mature numeric array, the number of mature individuals in
#'    each area, at each timestep, under all control rules, with all estimates
#'    of natural mortality.
#' @param Biomass numeric array, the total biomass in each area, at each time
#'    step, under all control rules, with all estimates of natural mortality,
#'    in kg.
#' @param Fishing logical value, is fishing occurring? Default value is TRUE.
#' @param Nat_mortality numeric vector, the estimates of natural mortality.
#' @param Recruitment_mode character value, values can be:
#'    'closed' - the recruits in each area originate from adults in that area.
#'    'pool' - the recruits in each area come from a pool of larvae produced by
#'       adults in all areas.
#'    Default value is 'pool'.
#'
#' @return numeric arrays, with updated values of fishing mortality (FM),
#'    numbers at age (N), Abundance_all, Abundance_mature, Biomass, and spawning
#'    stock biomass (SSB).
#' @export
#'
#' @examples
#' SSB <- array(rep(548, 5*70*6*3), c(5, 70, 6, 3))
#' Biomass <- array(rep(568, 5*70*6*3), c(5, 70, 6, 3))
#' Abundance_all <- array(rep(3400, 5*70*6*3), c(5, 70, 6, 3))
#' Abundance_mature <- array(rep(2800, 5*70*6*3), c(5, 70, 6, 3))
#' N <- array(rep(10, 34*5*70*6*3), c(34, 5, 70, 6, 3))
#' FM <- array(rep(0.2, 34*5*70*6*3), c(34, 5, 70, 6, 3))
#' L <- length_age(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#'    L2 = 47.95, K = 0.2022, All_ages = FALSE)
#' W <- weight(L, WA = 1.68e-5, WB = 3)
#' Mat <- maturity(Rec_age = 2, Max_age = 35, K_mat = -0.4103, L, L50 = 39.53)
#' NuR <- array(rnorm(5*70*6*3, 0, 0.5), c(5, 70, 6, 3))
#' Eps <- epsilon(A = 5, TimeT = 70, CR = 6, NM = 3, NuR, Rho_R = 0)
#' E <- array(rep(1, 5*70*6*3), c(5, 70, 6, 3))
#' S <- selectivity(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#'    L2 = 47.95, K = 0.2022, Fleets = c('sport', 'hook', 'trawl'),
#'    A50_up = c(2, 5, 10), A50_down = c(6, 16, 35), Alpha = c(0.33, 0.6, 0.64),
#'    F_fin = c(0.25, 0.06, 1), Beta = c(1.2, 0.6, 0), Cf = c(0.71, 0.28, 0.01))
#' pop_dynamics(a = 1, t = 3, cr = 1, nm = 2, Rec_age = 2, Max_age = 35, SSB, N,
#'    W, Mat, A = 5, R0 = 1e+5, H = 0.65, B0 = 1e+5/1.1, Eps, Sigma_R = 0.5,
#'    Fb = 0.2, E, S, NM = 3, FM, A50_mat = 8, Abundance_all, Abundance_mature,
#'    Biomass, Fishing = TRUE, Nat_mortality = c(0.09, 0.14, 0.19),
#'    Recruitment_mode = 'pool')
pop_dynamics <- function(a, t, cr, nm, Rec_age, Max_age, SSB, N, W, Mat, A = 5,
                         R0 = 1e+5, H, B0, Eps, Sigma_R, Fb, E, S, NM = 3, FM,
                         A50_mat, Abundance_all, Abundance_mature, Biomass,
                         Fishing = T, Nat_mortality, Recruitment_mode = 'pool') {

  ###### Error handling ########################################################

  # classes of variables
  if (a %% 1 != 0) {stop('a must be an integer value.')}
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (nm %% 1 != 0) {stop('nm must be an integer value.')}
  if (Rec_age %% 1 != 0) {stop('Rec_age must be an integer value.')}
  if (Max_age %% 1 != 0) {stop('Max_age must be an integer value.')}
  if (!is.numeric(SSB)) {stop('SSB must be a numeric array.')}
  if (!is.numeric(N)) {stop('N must be a numeric array.')}
  if (!is.numeric(W)) {stop('W must be a numeric vector.')}
  if (!is.numeric(Mat)) {stop('Mat must be a numeric vector.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (R0 <= 0) {stop('R0 must be greater than 0.')}
  if (H <= 0 || H > 1) {stop('H must be between 0 and 1.')}
  if (!is.numeric(B0)) {stop('B0 must be a numeric value.')}
  if (!is.numeric(Eps)) {stop('Eps must be a numeric array.')}
  if (!is.numeric(Sigma_R)) {stop('Sigma_R must be a numeric array.')}
  if (!is.numeric(Fb)) {stop('Fb must be a numeric value.')}
  if (!is.numeric(E)) {stop('E must be a numeric array.')}
  if (!is.numeric(S)) {stop('S must be a numeric vector.')}
  if (NM %% 1 != 0) {stop('NM must be an integer value.')}
  if (!is.numeric(FM)) {stop('FM must be a numeric array.')}
  if (A50_mat %% 1 != 0) {stop('A50_mat must be an integer value.')}
  if (!is.numeric(Abundance_all)) {
    stop('Abundance_all must be a numeric array.')}
  if (!is.numeric(Abundance_mature)) {
    stop('Abundance_mature must be a numeric array.')}
  if (!is.numeric(Biomass)) {stop('Biomass must be a numeric array.')}
  if (!is.logical(Fishing)) {stop('Fishing must be a logical value.')}
  if (!is.numeric(Nat_mortality)) {stop('Nat_mortality must be a numeric vector.')}
  if (!is.character(Recruitment_mode)) {
    stop('Recruitment mode must be a character value.')}

  # acceptable values
  if (a <= 0) {stop('a must be greater than 0.')}
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (nm <= 0 || nm > 3) {
    stop('nm must be greater than 0 and less than or equal to 3.')}
  if (Rec_age <= 0) {stop('Rec_age must be greater than 0.')}
  if (sum(SSB < 0) > 0) {stop('All values in SSB must be greater than or equal to 0.')}
  if (sum(N < 0) > 0) {stop('All values in N must be greater than or equal to 0.')}
  if (sum(W <= 0) > 0) {stop('All values in W must be greater than 0.')}
  if (sum(Mat <= 0) > 0 || sum(Mat > 1) > 0) {
    stop('All values in Mat must be between 0 and 1.')}
  if (A <= 0) {stop('A must be greater than 0.')}
  if (R0 <= 0) {stop('R0 must be greater than 0.')}
  if (H <= 0 || H > 1) {stop('H must be between 0 and 1.')}
  if (B0 <= 0) {stop('B0 must be greater than 0.')}
  if (Sigma_R <= 0) {stop('Sigma_R must be greater than 0.')}
  if (Fb < 0) {stop('Fb must be greater than or equal to 0.')}
  if (sum(E < 0) > 0) {stop('All values in E must be greater than or equal to 0.')}
  if (sum(S < 0) > 0) {stop('All values in S must be greater than or equal to 0.')}
  if (NM != 1 && NM != 3) {stop('NM must be equal to 1 or 3.')}
  if (sum(FM < 0) > 0 || sum(FM > 1) > 0) {
    stop('All values in FM must be between 0 and 1.')}
  if (A50_mat <= 0) {stop('A50_mat must be greater than 0.')}
  if (sum(Abundance_all < 0) > 0) {
    stop('All values in Abundance_all must be greater than or equal to 0.')}
  if (sum(Abundance_mature < 0) > 0) {
    stop('All values in Abundance_mature must be greater than or equal to 0.')}
  if (sum(Biomass < 0) > 0) {
    stop('All values in Biomass must be greater than or equal to 0.')}
  if (sum(Nat_mortality <= 0) > 0 || sum(Nat_mortality > 1) > 0) {
    stop('All values in Nat_mortality must be between 0 and 1.')}
  if (Recruitment_mode != 'pool' && Recruitment_mode != 'closed') {
    stop('Recruitment_mode must be either "pool" or "closed".')}

  # relational values
  if (Rec_age >= Max_age) {stop('Rec_age must be less than Max_age.')}
  if(dim(N)[1] != dim(FM)[1]) {
    stop('N or FM has an incorrect number of age classes.')}
  if(dim(N)[2] != dim(SSB)[1] || dim(N)[2] != dim(FM)[2] || dim(N)[2] != dim(E)[1]) {
    stop('N, SSB, FM, or E has an incorrect number of areas.')}
  if(dim(N)[3] != dim(SSB)[2] || dim(N)[3] != dim(FM)[3]|| dim(N)[3] != dim(E)[2]) {
    stop('N, SSB, FM, or E has an incorrect number of time steps.')}
  if(dim(N)[4] != dim(SSB)[3] || dim(N)[4] != dim(FM)[4]|| dim(N)[4] != dim(E)[3]) {
    stop('N, SSB, FM, or E has an incorrect number of control rules.')}
  if(dim(N)[5] != dim(SSB)[4] || dim(N)[5] != dim(FM)[5]|| dim(N)[5] != dim(E)[4]) {
    stop('N, SSB, FM, or E has an incorrect number of values in Nat_mortality.')}
  if (A != dim(N)[2]) {stop('N has the wrong number of areas.')}
  if (t > dim(N)[3]) {stop('The given "t" value is too high for N.')}
  if (cr > dim(N)[4]) {stop('The given "cr" value is too high for N.')}
  if (NM != dim(N)[5]) {stop('N has the wrong number of natural mortality estimates.')}
  if (a > A) {stop('a must be less than or equal to A.')}
  if (dim(Abundance_all)[1] != dim(Abundance_mature)[1] || dim(Abundance_all)[1] != A) {
    stop('Abundance_all or Abundance_mature has an incorrect number of areas.')}
  if (dim(Abundance_all)[2] != dim(Abundance_mature)[2]) {
    stop('Abundance_all or Abundance_mature has an incorrect number of time steps.')}
  if (dim(Abundance_all)[3] != dim(Abundance_mature)[3]) {
    stop('Abundance_all or Abundance_mature has an incorrect number of control rules.')}
  if (dim(Abundance_all)[4] != dim(Abundance_mature)[4]) {
    stop('Abundance_all or Abundance_mature has an incorrect number of natural
         mortality estimates.')}
  if (length(W) != length(S) || length(W) != length(Mat)) {
    stop('W, S, or Mat has the wrong number of age classes.')}

  ##############################################################################

  # range of ages
  ages <- Rec_age:Max_age
  num <- length(ages)

  # Calculate fishing mortality
  if (Fishing == T) {
    FM[, a, t, cr, nm] <- f_mortality(a, t, cr, nm, FM, A, Fb, E, S)
  } else { FM[, a, t, cr, nm] <- 0 }

  ##### Step population foward in time

  # Calculate recruitment and add recruits to population
  N[1, a, t, cr, nm] <- recruitment(a, t, cr, nm, SSB, A, R0, H, B0, Eps,
                                    Sigma_R, Rec_age, Recruitment_mode)

  # Ages rec_age + 1 to max_age - 1
  for (i in 2:(num - 1)) {
    N[i, a, t, cr, nm] <- N[i - 1, a, t - 1, cr, nm] * exp(-1 * (FM[i - 1, a, t - 1, cr, nm] + Nat_mortality[nm]))
  }

  # Final age bin
  N[num, a, t, cr, nm] <- N[num - 1, a, t - 1, cr, nm] * exp(-1 * (FM[num - 1, a, t - 1, cr, nm] + Nat_mortality[nm])) +
    N[num, a, t - 1, cr, nm] * exp(-1 * (FM[num, a, t - 1, cr, nm] + Nat_mortality[nm]))

  # Calculate abundance of all fish
  Abundance_all[a, t, cr, nm] <- sum(N[, a, t, cr, nm])

  # Calculate abundance of mature fish
  Abundance_mature[a, t, cr, nm] <- sum(N[A50_mat:(Max_age - Rec_age + 1), a, t, cr, nm])

  # Calculate biomass of all fish
  Biomass[a, t, cr, nm] <- sum(N[, a, t, cr, nm] * W)

  # Calculate spawning stock biomass
  SSB[a, t, cr, nm] <- sum(N[, a, t, cr, nm]*W*Mat)

  output <- list(FM[, a, t, cr, nm], N[, a, t, cr, nm],
                 Abundance_all[a, t, cr, nm], Abundance_mature[a, t, cr, nm],
                 Biomass[a, t, cr, nm], SSB[a, t, cr, nm])

  return(output)

}
