#' Population Dynamics
#'
#' \code{pop_dynamics} returns updated values for fishing mortality, numbers at
#'    age, abundance of all and mature individuals, biomass, and spawning stock
#'    biomass in a specific area, at a specific time step, under a specific
#'    control rule, and with a specific estimate of natural mortality
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param nm temporary numeric value, the current natural mortality estimate.
#' @param NM numeric value, the total number of natural mortality estimates.
#' @param fdr temporary numeric value, the current final target density ratio.
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
#' @param Fb numeric value, the historical fishing effort for the fished species.
#' @param E numeric array, the relative fishing effort displayed in each area,
#'    at each time step, under each control rule, and for each natural mortality
#'    estimate.
#' @param S numeric vector, the selectivities at age from age at recruitment to
#'    maximum age, on the interval (0, 1).
#' @param FM numeric array that corresponds to the fishing mortality at each
#'    age in each area, at each timestep, under all control rules, with all
#'    estimates of natural mortality.
#' @param A50_mat numeric value, the first age at which 50\% or more individuals
#'    are estimated to be mature. on the interval (Rec_age, Max_age).
#' @param Abundance numeric array, the total number of all and/or mature
#'    individuals in each area, at each timestep, under all control rules, with
#'    all estimates of natural mortality.
#' @param Biomass numeric array, the total biomass in each area, at each time
#'    step, under all control rules, with all estimates of natural mortality,
#'    in kg.
#' @param Fishing logical value, is fishing occurring? Default value is TRUE.
#' @param Nat_mortality numeric vector, the estimates of natural mortality.
#' @param R numeric vector, the estimated recruits coming into each area at the
#'    current time step, under the current control rule, and with the current
#'    estimate of natural mortality.
#' @param Ind_sampled character value, the individuals to be sampled to
#'    calculate density ratio. Values can be:
#'    'all' - sample all individuals.
#'    'mature' - sample only mature individuals.
#'    Default value is 'all'.
#'
#' @return numeric arrays, with updated values of fishing mortality (FM),
#'    numbers at age (N), Abundance, Biomass, and spawning stock biomass (SSB).
#' @export
#'
#' @examples
#' n = 34; A = 5; TimeT = 70; CR = 6; NM = 3; FDR = 4
#' SSB <- array(rep(548, A*TimeT*CR*NM*FDR), c(A, TimeT, CR, NM, FDR))
#' Biomass <- array(rep(568, A*TimeT*CR*NM*FDR), c(A, TimeT, CR, NM, FDR))
#' Abundance <- array(rep(3400, A*TimeT*CR*NM*FDR*1), c(A, TimeT, CR, NM, FDR, 1))
#' N <- array(rep(10, n*A*TimeT*CR*NM*FDR), c(n, A, TimeT, CR, NM, FDR))
#' FM <- array(rep(0.2, n*A*TimeT*CR*NM*FDR), c(n, A, TimeT, CR, NM, FDR))
#' L <- length_age(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#'    L2 = 47.95, K = 0.2022, All_ages = FALSE)
#' W <- weight(L, WA = 1.68e-5, WB = 3)
#' Mat <- maturity(Rec_age = 2, Max_age = 35, K_mat = -0.4103, L, L50 = 39.53)
#' E <- array(rep(1, A*TimeT*CR*NM*FDR), c(A, TimeT, CR, NM, FDR))
#' S <- selectivity(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#'    L2 = 47.95, K = 0.2022, Fleets = c('sport', 'hook', 'trawl'),
#'    A50_up = c(2, 5, 10), A50_down = c(6, 16, 35), Alpha = c(0.33, 0.6, 0.64),
#'    F_fin = c(0.25, 0.06, 1), Beta = c(1.2, 0.6, 0), Cf = c(0.71, 0.28, 0.01))
#' NuR <- array(rnorm(A*TimeT*CR*NM*FDR, 0, 0.5), c(A, TimeT, CR, NM, FDR))
#' Eps <- epsilon(A = 5, TimeT = 70, CR = 6, NM = 3, FDR = 4, NuR, Rho_R = 0)
#' R <- recruitment(t = 3, cr = 1, NM, fdr = 1, SSB, A = 5, R0 = 1e+5,
#'    H = 0.65, B0 = 1e+5/1.1, Eps, Sigma_R = 0.5, Rec_age = 2,
#'    Recruitment_mode = 'pool', LDP = 0.1)
#' pop_dynamics(t = 3, cr = 1, nm = 1, fdr = 1, Rec_age = 2, Max_age = 35, SSB,
#'    N, W, Mat, A = 5, Fb = 0.2, E, S, FM, NM, A50_mat = 8, Biomass,
#'    Abundance, Fishing = TRUE, Nat_mortality = c(0.09, 0.14, 0.19), R,
#'    Ind_sampled = 'all')
pop_dynamics <- function(t, cr, nm, fdr, Rec_age, Max_age, SSB, N, W, Mat,
                         A = 5, Fb, E, S, FM, NM, A50_mat, Biomass,
                         Abundance, Fishing = T, Nat_mortality, R,
                         Ind_sampled = 'all') {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (nm %% 1 != 0) {stop('nm must be an integer value.')}
  if (NM %% 1 != 0) {stop('NM must be an integer value.')}
  if (fdr %% 1 != 0) {stop('fdr must be an integer value.')}
  if (Rec_age %% 1 != 0) {stop('Rec_age must be an integer value.')}
  if (Max_age %% 1 != 0) {stop('Max_age must be an integer value.')}
  if (!is.numeric(SSB)) {stop('SSB must be a numeric array.')}
  if (!is.numeric(N)) {stop('N must be a numeric array.')}
  if (!is.numeric(W)) {stop('W must be a numeric vector.')}
  if (!is.numeric(Mat)) {stop('Mat must be a numeric vector.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (!is.numeric(Fb)) {stop('Fb must be a numeric value.')}
  if (!is.numeric(E)) {stop('E must be a numeric array.')}
  if (!is.numeric(S)) {stop('S must be a numeric vector.')}
  if (!is.numeric(FM)) {stop('FM must be a numeric array.')}
  if (A50_mat %% 1 != 0) {stop('A50_mat must be an integer value.')}
  if (!is.numeric(Abundance)) {stop('Abundance must be a numeric array.')}
  if (!is.numeric(Biomass)) {stop('Biomass must be a numeric array.')}
  if (!is.logical(Fishing)) {stop('Fishing must be a logical value.')}
  if (!is.numeric(Nat_mortality)) {stop('Nat_mortality must be a numeric vector.')}

  # acceptable values
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (nm <= 0 || nm > 3) {
    stop('nm must be greater than 0 and less than or equal to 3.')}
  if (NM <= 0 || NM > 3) {
    stop('NM must be greater than 0 and less than or equal to 3.')}
  if (fdr <= 0) {stop('fdr must be greater than 0.')}
  if (Rec_age <= 0) {stop('Rec_age must be greater than 0.')}
  if (sum(SSB < 0) > 0) {stop('All values in SSB must be greater than or equal to 0.')}
  if (sum(N < 0) > 0) {stop('All values in N must be greater than or equal to 0.')}
  if (sum(W <= 0) > 0) {stop('All values in W must be greater than 0.')}
  if (sum(Mat <= 0) > 0 || sum(Mat > 1) > 0) {
    stop('All values in Mat must be between 0 and 1.')}
  if (A <= 0) {stop('A must be greater than 0.')}
  if (Fb < 0) {stop('Fb must be greater than or equal to 0.')}
  if (sum(E < 0) > 0) {stop('All values in E must be greater than or equal to 0.')}
  if (sum(S < 0) > 0) {stop('All values in S must be greater than or equal to 0.')}
  if (sum(FM < 0) > 0) {
    stop('All values in FM must be greater than or equal to 0.')}
  if (A50_mat <= 0) {stop('A50_mat must be greater than 0.')}
  if (sum(Abundance < 0) > 0) {
    stop('All values in Abundance must be greater than or equal to 0.')}
  if (sum(Biomass < 0) > 0) {
    stop('All values in Biomass must be greater than or equal to 0.')}
  if (sum(Nat_mortality <= 0) > 0 || sum(Nat_mortality > 1) > 0) {
    stop('All values in Nat_mortality must be between 0 and 1.')}
  if (sum(R < 0) > 0) {stop('All values in R must be greater than or equal to 0.')}

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
  if(dim(N)[6] != dim(SSB)[5] || dim(N)[6] != dim(FM)[6]|| dim(N)[6] != dim(E)[5]) {
    stop('N, SSB, FM, or E has an incorrect number of final density ratios.')}
  if (dim(N)[2] != A) {stop('N has the wrong number of areas.')}
  if (t > dim(N)[3]) {stop('The given "t" value is too high for N.')}
  if (cr > dim(N)[4]) {stop('The given "cr" value is too high for N.')}
  if (dim(N)[5] != NM) {stop('N has the wrong number of natural mortality estimates.')}
  if (fdr > dim(N)[6]) {stop('The given "fdr" value is too high for N.')}
  if (dim(Abundance)[1] != A) {
    stop('Abundance has an incorrect number of areas.')}
  if (t > dim(Abundance)[2]) {
    stop('Given "t" value is too high for Abundance.')}
  if (cr > dim(Abundance)[3]) {
    stop('Given "cr" value is too high for Abundance.')}
  if (dim(Abundance)[4] != NM) {
    stop('Abundance has an incorrect number of natural mortality estimates.')}
  if (fdr > dim(Abundance)[5]) {
    stop('Given "fdr" value is too high for Abundance.')}
  if (length(W) != length(S) || length(W) != length(Mat)) {
    stop('W, S, or Mat has the wrong number of age classes.')}

  ##############################################################################

  # for (nm in 1:NM) {

    # range of ages
    ages <- Rec_age:Max_age
    num <- length(ages)

    # Calculate fishing mortality
    if (Fishing == T) {
      FM[, , t, cr, nm, fdr] <- f_mortality(t, cr, nm, fdr, FM, A, Fb, E, S)
    } else { FM[, , t, cr, nm, fdr] <- 0 }

    ##### Step population foward in time

    # Calculate recruitment and add recruits to population
    if (NM > 1) {N[1, , t, cr, nm, fdr] <- R[1:A, nm]
    } else {N[1, , t, cr, nm, fdr] <- R[1:A]}

    # N[1, , t, cr, nm, fdr] <- R[1:A, nm]

    # Ages rec_age + 1 to max_age - 1
    for (i in 2:(num - 1)) {
      N[i, , t, cr, nm, fdr] <- N[i - 1, , t - 1, cr, nm, fdr] * exp(-1 * (FM[i - 1, , t - 1, cr, nm, fdr] + Nat_mortality[nm]))
    }

    # Final age bin
    N[num, , t, cr, nm, fdr] <- N[num - 1, , t - 1, cr, nm, fdr] * exp(-1 * (FM[num - 1, , t - 1, cr, nm, fdr] + Nat_mortality[nm])) +
      N[num, , t - 1, cr, nm, fdr] * exp(-1 * (FM[num, , t - 1, cr, nm, fdr] + Nat_mortality[nm]))

    if (A > 1) {

      # Calculate biomass of all fish
      Biomass[, t, cr, nm, fdr] <- colSums(N[, , t, cr, nm, fdr] * W)

      # Calculate spawning stock biomass
      SSB[, t, cr, nm, fdr] <- colSums(N[, , t, cr, nm, fdr]*W*Mat)

      # Abundance of all / mature individuals
      Abundance[, t, cr, nm, fdr, 1] <- colSums(N[, , t, cr, nm, fdr])
      if (Ind_sampled == 'mature' || is.null(Ind_sampled)) {
        Abundance[, t, cr, nm, fdr, 2] <- colSums(N[, , t, cr, nm, fdr])
      }

    } else if (A == 1) {

      # Calculate biomass of all fish
      Biomass[1, t, cr, nm, fdr] <- sum(N[, 1, t, cr, nm, fdr] * W)

      # Calculate spawning stock biomass
      SSB[1, t, cr, nm, fdr] <- sum(N[, 1, t, cr, nm, fdr]*W*Mat)

      # Abundance of all / mature individuals
      Abundance[1, t, cr, nm, fdr, 1] <- sum(N[, 1, t, cr, nm, fdr])
      if (Ind_sampled == 'mature' || is.null(Ind_sampled)) {
        Abundance[1, t, cr, nm, fdr, 2] <- sum(N[ages, 1, t, cr, nm, fdr])
      }

    }

  # }

  output <- list(FM[, , t, cr, nm, fdr],
                 N[, , t, cr, nm, fdr],
                 Biomass[, t, cr, nm, fdr],
                 SSB[, t, cr, nm, fdr],
                 Abundance[, t, cr, nm, fdr, ])

  return(output)

}
