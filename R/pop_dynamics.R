#' Population Dynamics
#'
#' \code{pop_dynamics} returns updated values for fishing mortality, numbers at
#'    age, abundance of all and mature individuals, biomass, and spawning stock
#'    biomass in a specific area, at a specific time step, under a specific
#'    control rule, and with a specific estimate of natural mortality
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param NM  numeric value, the total number of natural mortality estimates.
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
#' n = 34; A = 5; TimeT = 70; CR = 6; NM = 1; FDR = 4
#' SSB <- array(rep(548, A*TimeT*CR*FDR*NM), c(A, TimeT, CR, FDR, NM))
#' Biomass <- array(rep(568, A*TimeT*CR*FDR), c(A, TimeT, CR, FDR))
#' Abundance <- array(rep(3400, A*TimeT*CR*FDR*1*NM), c(A, TimeT, CR, FDR, 1, NM))
#' N <- array(rep(10, n*A*TimeT*CR*FDR*NM), c(n, A, TimeT, CR, FDR, NM))
#' FM <- array(rep(0.2, n*A*TimeT*CR*FDR), c(n, A, TimeT, CR, FDR))
#' L <- length_age(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#'    L2 = 47.95, K = 0.2022, All_ages = FALSE)
#' W <- weight(L, WA = 1.68e-5, WB = 3)
#' Mat <- maturity(Rec_age = 2, Max_age = 35, K_mat = -0.4103, L, L50 = 39.53)
#' E <- array(rep(1, A*TimeT*CR*FDR), c(A, TimeT, CR, FDR))
#' S <- selectivity(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#'    L2 = 47.95, K = 0.2022, Fleets = c('sport', 'hook', 'trawl'),
#'    A50_up = c(2, 5, 10), A50_down = c(6, 16, 35), Alpha = c(0.33, 0.6, 0.64),
#'    F_fin = c(0.25, 0.06, 1), Beta = c(1.2, 0.6, 0), Cf = c(0.71, 0.28, 0.01))
#' NuR <- array(rnorm(A*TimeT*CR*FDR*NM, 0, 0.5), c(A, TimeT, CR, FDR, NM))
#' Eps <- epsilon(A, TimeT, CR, FDR, NM, NuR, Rho_R = 0)
#' R <- recruitment(t = 3, cr = 1, NM, fdr = 1, SSB, A, R0 = 1e+5, H = 0.65,
#'    B0 = 1e+5/1.1, Eps, Sigma_R = 0.5, Rec_age = 2,
#'    Recruitment_mode = 'pool', LDP = 0.1)
#' pop_dynamics(t = 3, cr = 1, NM, fdr = 1, Rec_age = 2, Max_age = 35, SSB,
#'    N, W, Mat, A = 5, Fb = 0.2, E, S, FM, A50_mat = 8, Biomass,
#'    Abundance, Fishing = TRUE, Nat_mortality = 0.14, R, Ind_sampled = 'all')
pop_dynamics <- function(t, cr, NM, fdr, Rec_age, Max_age, SSB, N, W, Mat,
                         A = 5, Fb, E, S, FM, A50_mat, Biomass,
                         Abundance, Fishing = T, Nat_mortality, R,
                         Ind_sampled = 'all') {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
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
    stop('N, SSB, FM, or E has an incorrect number of final density ratios.')}
  if (dim(N)[2] != A | dim(Abundance)[1] != A) {
    stop('N or Abundance has the wrong number of areas.')}
  if (t > dim(N)[3] | t > dim(Abundance)[2]) {
    stop('The given "t" value is too high for N or Abundance.')}
  if (cr > dim(N)[4] | cr > dim(Abundance)[3]) {
    stop('The given "cr" value is too high for N or Abundance.')}
  if (fdr > dim(N)[5] | fdr > dim(Abundance)[4]) {
    stop('The given "fdr" value is too high for N or Abundance.')}
  if (length(W) != length(S) || length(W) != length(Mat)) {
    stop('W, S, or Mat has the wrong number of age classes.')}

  ##############################################################################

  # range of ages
  ages <- Rec_age:Max_age
  num <- length(ages)

  # Calculate fishing mortality
  if (Fishing == T) {
    FM[, , t, cr, fdr] <- f_mortality(t, cr, fdr, FM, A, Fb, E, S)
  } else { FM[, , t, cr, fdr] <- 0 }

  ##### Step population foward in time

  for (nm in 1:NM) {

    # Calculate recruitment and add recruits to population
    N[1, , t, cr, fdr, nm] <- R[1:A, nm]

    # natural mortality arrays
    m <- ifelse(cr < 3, Nat_mortality[1], Nat_mortality[ceiling(cr / 2)])

    # Ages rec_age + 1 to max_age - 1
    for (i in 2:(num - 1)) {
      N[i, , t, cr, fdr, ] <- N[i - 1, , t - 1, cr, fdr, ] *
        exp(-1 * (FM[i - 1, , t - 1, cr, fdr] + m))
    }

    # Final age bin
    N[num, , t, cr, fdr, ] <- N[num - 1, , t - 1, cr, fdr, ] *
      exp(-1 * (FM[num - 1, , t - 1, cr, fdr] + m)) +
      N[num, , t - 1, cr, fdr, ] * exp(-1 * (FM[num, , t - 1, cr, fdr] + m))

    # Calculate biomass of all fish
    if (A > 1) {

      # Calculate biomass of all fish
      Biomass[, t, cr, fdr] <- colSums(N[, , t, cr, fdr, 1] * W)
      # Calculate spawning stock biomass
      SSB[, t, cr, fdr, nm] <- colSums(N[, , t, cr, fdr, nm]*W*Mat)

      # Abundance of all / mature individuals
      Abundance[, t, cr, fdr, 1, nm] <- colSums(N[, , t, cr, fdr, nm])
      if (Ind_sampled == 'mature' || is.null(Ind_sampled)) {
        Abundance[, t, cr, fdr, 2, nm] <- colSums(N[A50_mat:num, , t, cr, fdr, nm])
      }

    } else if (A == 1) {

      # Calculate biomass of all fish
      Biomass[1, t, cr, fdr] <- sum(N[, 1, t, cr, fdr, 1] * W)

      # Calculate spawning stock biomass
      SSB[1, t, cr, fdr, 1] <- sum(N[, 1, t, cr, fdr, 1]*W*Mat)

      # Abundance of all / mature individuals
      Abundance[1, t, cr, fdr, 1, 1] <- sum(N[, 1, t, cr, fdr, 1])
      if (Ind_sampled == 'mature' || is.null(Ind_sampled)) {
        Abundance[1, t, cr, fdr, 2, 1] <- sum(N[A50_mat:num, 1, t, cr, fdr, 1])
      }

    }
  }

  output <- list(FM[, , t, cr, fdr],
                 N[, , t, cr, fdr, ],
                 Biomass[, t, cr, fdr],
                 SSB[, t, cr, fdr, ],
                 Abundance[, t, cr, fdr, , ])

  return(output)

}
