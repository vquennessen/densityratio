#' Initialize arrays
#'
#' \code{initialize_arrays} returns arrays necessary for other functions in the
#'    base model: TimeT (total time the model runs), L (length at age),
#'    W (weight at age), S (selectivity at age), Mat (maturity at age), A50_mat
#'    (the first age at which 50\% or more of individuals are expected to be
#'    mature), CR (number of control rules to be compared), Nat_mortality
#'    (values of estimated natural mortality), NM (number of estimates of
#'    natural mortality), N (numbers at age array), SSB (spawning stock biomass
#'    array), Abundance (abundance of all/mature individuals array),
#'    Eps (recruitment error), B0 (unfished biomass), Count (estimated counts
#'    based on sampling array), Sigma_S (sampling random variable standard
#'    deviation), NuS (sampling random variable), Delta (constant of
#'    proportionality), Gamma, FM (fishing mortality), E (effort), Catch, Yield,
#'    Rel_Biomass (relative biomass after reserve implementation),
#'    Rel_yield (relative yield after reserve implementation), Rel_SSB (relative
#'    spawning stock biomass after reserve implementation), Density_ratio.
#'
#' @param A numeric value, the number of total areas in the model. Default value
#'    is 5.
#' @param MPA numeric value, the area number to be designated as a marine
#'    reserve at Time1. Default value is 3.
#' @param Final_DRs numeric vector, the final target density ratios.
#' @param Time1 numeric value, the number of years to run the model before a
#'    marine reserve is implemented. Default value is 50.
#' @param Time2 numeric value, the number of years to run the model after a
#'    marine reserve is implemented. Default value is 20.
#' @param R0 numeric value, the unfished recruitment, set arbitrarily. Default
#'    value is 1e+5.
#' @param Rec_age numeric value, the age at recruitment, in years.
#' @param Max_age numeric value, the maximum age of fish or total lifespan, in
#'    years.
#' @param A1 numeric value, age 1 of female fish, in years.
#' @param L1 numeric value, the length of females at age 1, in cm.
#' @param A2 numeric value, age 2 of female fish, in years.
#' @param L2 numeric value, the length of females at age 2, in cm.
#' @param K numeric value, the von Bertalanffy growth parameter for females.
#' @param WA numeric value, the coefficient in the weight at length equation.
#' @param WB numeric value, the exponent in the weight at length equation.
#' @param K_mat numeric value, the slope of the maturity curve.
#' @param Fb numeric value, the historical fishing effort for the fished species
#'    on the interval (0, 1).
#' @param L50 numeric value, the length at 50\% maturity, in cm.
#' @param Sigma_R numeric value, the recruitment standard deviation.
#' @param Rho_R numeric value, the recruitment autocorrelation on the interval
#'    (-1, 1). Default value is 0.
#' @param Fleets character vector, the names of all fishing fleets that
#'    contribute to selectivity at age on the interval (0, 1).
#' @param Alpha numeric vector, the alpha values for each fishing fleet in
#'    Fleets on the interval (0, 1), gives the relative slope of the upcurve.
#' @param A50_up numeric vector, the age values at which selectivity is equal to
#'    0.5 for the upcurve each fishing fleet in Fleets on the interval (0, 1).
#' @param A50_down numeric vector, the age values at which selectivity is equal
#'    to 0.5 for the downcurve each fishing fleet in Fleets on the interval
#'    (0, 1).
#' @param F_fin numeric vector, the final selectivity values for fleets on the
#'    interval (0, 1).
#' @param Beta numeric vector, the beta values for each fishing fleet in Fleets
#'    on the interval (0, 1), gives the relative slope of the downcurve if
#'    dome-shaped.
#' @param Cf numeric vector, the fraction of the whole fishery represented by
#'    each fleet.
#' @param P numeric value, the proportion of positive transects during sampling.
#' @param X numeric value, the average value of individuals seen during positive
#'    transects.
#' @param SP numeric value, the standard deviation of individuals seen during
#'    positive transects.
#' @param M numeric value, the natural mortality on the interval (0, 1).
#' @param Control_rules numeric vector, the control rules to be compared.
#' @param Phi numeric value, the unfished recruits per spawner.
#' @param Stochasticity logical vector, does recruitment contain a stochastic
#'    component? Default value is TRUE.
#' @param D numeric value, the current depletion of the stock, on the interval
#'    (0, 1).
#' @param Transects numerical value, the number of sampling transects conducted
#'    in each area to estimate density ratio. Default value is 24.
#' @param H numeric value, the steepness of the stock-recruitment curve.
#' @param Surveys logical value, are surveys being conducted? Default value is
#'    TRUE.
#' @param Fishing logical value, is fishing occurring? Default value is TRUE.
#' @param M_Error numeric value, the error between estimated and correct natural
#'    mortality.
#' @param Sampling_Error logical value, is there any error in sampling? Default
#'    value is TRUE.
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
#' @param Ind_sampled character value, the individuals to be sampled to
#'    calculate density ratio. Values can be:
#'    'all' - sample all individuals.
#'    'mature' - sample only mature individuals.
#'    Default value is 'all'.
#'
#' @return initalizes arrays necessary for other functions in the base model,
#'    including Inside, Outside, FDR, TimeT, L, W, S, Mat, A50_mat, CR,
#'    Nat_mortality, NM, N, SSB, Biomass, Eps, B0, Count, Sigma_S, NuS, Delta,
#'    Gamma, FM, E, Catch, Yield, Density_ratio, ENM, Abundance
#' @export
#'
#' @importFrom stats rnorm
#'
#' @examples
#' initialize_arrays(A = 5,  MPA = 3, Final_DRs = c(0.2), Time1 = 50,
#'    Time2 = 20, R0 = 1e+5, Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21,
#'    A2 = 15, L2 = 47.95, K = 0.2022, WA = 1.68e-5, WB = 3, K_mat = -0.4103,
#'    Fb = 0.2, L50 = 39.53, Sigma_R = 0.5, Rho_R = 0,
#'    Fleets = c('sport', 'hook', 'trawl'), Alpha = c(0.33, 0.6, 0.64),
#'    A50_up = c(2, 5, 10), A50_down = c(6, 16, 35), F_fin = c(0.25, 0.06, 1),
#'    Beta = c(1.2, 0.6, 0), Cf = c(0.71, 0.28, 0.01), P = 0.77, X = 15.42,
#'    SP = 16.97, M = 0.14, Control_rules= c(1:6), Phi = 1.1,
#'    Stochasticity = TRUE, D = 0.488, Transects = 24, H = 0.65, Surveys = TRUE,
#'    Fishing = TRUE, M_Error = 0.05, Sampling_Error = TRUE,
#'    Recruitment_mode = 'pool', LDP = 0.1, Ind_sampled = 'all')
initialize_arrays <- function(A = 5, MPA = 3, Final_DRs, Time1 = 50, Time2 = 20,
                              R0 = 1e+5, Rec_age, Max_age, A1, L1, A2, L2, K,
                              WA, WB, K_mat, Fb, L50, Sigma_R, Rho_R = 0,
                              Fleets, Alpha, A50_up, A50_down, F_fin, Beta, Cf,
                              P, X, SP, M, Control_rules, Phi,
                              Stochasticity = TRUE, D, Transects = 24, H,
                              Surveys = TRUE, Fishing = TRUE, M_Error,
                              Sampling_Error = TRUE, Recruitment_mode = 'pool',
                              LDP = 0.1, Ind_sampled = 'all') {

  ###### Error handling ########################################################

  # classes of variables
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (MPA %% 1 != 0) {stop('MPA must be an integer value.')}
  if (!is.numeric(Final_DRs)) {stop('Final_DRs must be a numeric vector.')}
  if (Time1 %% 1 != 0) {stop('Time1 must be an integer value.')}
  if (Time2 %% 1 != 0) {stop('Time2 must be an integer value.')}
  if (R0 %% 1 != 0) {stop('R0 must be an integer value.')}
  if (Rec_age %% 1 != 0) {stop('Rec_age must be an integer value.')}
  if (Max_age %% 1 != 0) {stop('Max_age must be an integer value.')}
  if (A1 %% 1 != 0) {stop('A1 must be an integer value.')}
  if (!is.numeric(L1)) {stop('L1 must be a numeric value.')}
  if (A2 %% 1 != 0) {stop('A2 must be an integer value.')}
  if (!is.numeric(L2)) {stop('L2 must be a numeric value.')}
  if (!is.numeric(K)) {stop('K must be a numeric value.')}
  if (!is.numeric(WA)) {stop('WA must be a numeric value.')}
  if (!is.numeric(WB)) {stop('WB must be a numeric value.')}
  if (!is.numeric(K_mat)) {stop('K_mat must be a numeric value.')}
  if (!is.numeric(Fb)) {stop('Fb must be a numeric value.')}
  if (!is.numeric(L50)) {stop('L50 must be a numeric value.')}
  if (!is.numeric(Sigma_R)) {stop('Sigma_R must be a numeric value.')}
  if (!is.numeric(Rho_R)) {stop('Rho_R must be a numeric value.')}
  if (!is.character(Fleets)) {stop('Fleets must be a character vector.')}
  if (sum(A50_up %% 1 != 0) != 0) {stop('A50_up must be a vector of integers.')}
  if (sum(A50_down %% 1 != 0) != 0) {stop('A50_down must be a vector of integers.')}
  if (!is.numeric(Alpha)) {stop('Alpha must be a numeric vector.')}
  if (!is.numeric(F_fin)) {stop('F_fin must be a numeric vector.')}
  if (!is.numeric(Beta)) {stop('Beta must be a numeric vector.')}
  if (!is.numeric(Cf)) {stop('Cf must be a numeric vector.')}
  if (!is.numeric(P)) {stop('P must be a numeric value.')}
  if (!is.numeric(X)) {stop('X must be a numeric value.')}
  if (!is.numeric(SP)) {stop('SP must be a numeric value.')}
  if (!is.numeric(M)) {stop('M must be a numeric value.')}
  if (sum(Control_rules %% 1 != 0) != 0) {
    stop('Control_rules must be a vector of integers.')}
  if (!is.numeric(Phi)) {stop('Phi must be a numeric value.')}
  if (!is.logical(Stochasticity)) {
    stop('Stochasticity must be a logical value.')}
  if (!is.numeric(D)) {stop('D must be a numeric value.')}
  if (Transects %% 1 != 0) {stop('Transects must be an integer value.')}
  if (!is.numeric(H)) {stop('H must be a numeric value.')}
  if (!is.logical(Surveys)) {stop('Surveys must be a logical value.')}
  if (!is.logical(Fishing)) {stop('Fishing must be a logical value.')}
  if (!is.numeric(M_Error)) {stop('M_Error must be a numeric value.')}
  if (!is.logical(Sampling_Error)) {
    stop('Sampling_Error must be a logical value.')}
  if (!is.character(Recruitment_mode)) {
    stop('Recruitment mode must be a character value.')}
  if (!is.numeric(LDP)) {stop('LDP must be a numeric value.')}

  # acceptable values
  if (A <= 0) {stop('A must be greater than 0.')}
  if (MPA <= 0) {stop('MPA must be greater than 0.')}
  if (sum(Final_DRs <= 0) > 0) {
    stop('All values in Final_DRs must be greater than 0.')}
  if (Time1 <= 1) {stop('Time1 must be greater than 1.')}
  if (Time2 <= 1) {stop('Time2 must be greater than 1.')}
  if (R0 <= 0) {stop('R0 must be greater than 0.')}
  if (Rec_age <= 0) {stop('Rec_age must be greater than 0.')}
  if (A1 <= 0) {stop('A1 must be greater than 0.')}
  if (L1 <= 0) {stop('L1 must be greater than 0.')}
  if (K <= 0) {stop('K must be greater than 0.')}
  if (WA <= 0) {stop('WA must be greater than 0.')}
  if (WB <= 0) {stop('WB must be greater than 0.')}
  if (K_mat >= 0) {stop('K_mat must be less than 0.')}
  if (Fb < 0) {stop('Fb must be greater than or equal to 0.')}
  if (L50 <= 0) {stop('L50 must be greater than 0.')}
  if (Sigma_R <= 0) {stop('Sigma_R must be greater than 0.')}
  if (Rho_R < -1 || Rho_R > 1) {stop('Rho_R must be between -1 and 1.')}
  if (sum(A50_up <= 0) > 0) {stop('All values in A50_up must be greater than 0.')}
  if (sum(A50_down < 0) > 0) {
    stop('All values in A50_down must be greater than 0.')}
  if (sum(Alpha < 0) > 0) {
    stop('All values in Alpha must be greater than 0.')}
  if (sum(Beta < 0) > 0) {
    stop('All values in Beta must be greater than or equal to 0.')}
  if (sum(Cf <= 0) > 0) {stop('All values in Cf must be greater than 0.')}
  if (P < -1 || P > 1) {stop('P must be between -1 and 1.')}
  if (X < 0) {stop('X must be greater than or equal to 0.')}
  if (SP < 0) {stop('SP must be greater than or equal to 0.')}
  if (M <= 0 || M > 1) {stop('M must be between 0 and 1.')}
  if (sum(Control_rules < 1) > 0) {
    stop('All values in Control_rules must be greater than or equal to 1.')}
  if (Phi <= 0) {stop('Phi must be greater than 0.')}
  if (D <= 0 || D > 1) {stop('D must be between 0 and 1.')}
  if (Transects <= 0) {stop('Transects must be greater than 0.')}
  if (H <= 0 || H > 1) {stop('H must be between 0 and 1.')}
  if (M_Error < 0) {stop('M_Error must be greater than or equal to 0.')}
  if (Recruitment_mode != 'pool' && Recruitment_mode != 'closed' &&
      Recruitment_mode != 'regional_DD' && Recruitment_mode != 'local_DD') {
    stop('Recruitment_mode must be either "pool", "closed", "regional_DD", or
         "local_DD".')}
  if (LDP < 0) {stop('LDP must be greater than or equal to 0.')}

  # relational values
  if (MPA > A) {stop('MPA must be less than or equal to A.')}
  if (Rec_age >= Max_age) {stop('Rec_age must be less than Max_age.')}
  if (A1 >= A2) {stop('A1 must be less than A2.')}
  if (L1 >= L2) {stop('L1 must be less than L2.')}
  if (length(A50_up) != length(Fleets)) {
    stop('A50_up must have the same number of elements as Fleets.')}
  if (length(A50_down) != length(Fleets)) {
    stop('A50_down must have the same number of elements as Fleets.')}
  if (length(Alpha) != length(Fleets)) {
    stop('Alpha must have the same number of elements as Fleets.')}
  if (length(F_fin) != length(Fleets)) {
    stop('F_fin must have the same number of elements as Fleets.')}
  if (length(Beta) != length(Fleets)) {
    stop('Beta must have the same number of elements as Fleets')}
  if (length(Cf) != length(Fleets)) {
    stop('Cf must have the same number of elements as Fleets.')}

  ##############################################################################

  # set areas in and out of marine reserves
  areas <- 1:A
  Inside <- areas[MPA]
  Outside <- areas[-MPA]

  # final density ratios
  FDR <- length(Final_DRs)

  # total amount of timesteps (years)
  TimeT <- Time1 + Time2

  # ages for which fish have recruited
  ages <- Rec_age:Max_age
  num <- length(ages)

  # Length at age
  # Dimensions = 1 * age
  L <- length_age(Rec_age, Max_age, A1, L1, A2, L2, K, All_ages = F)

  # Weight at age
  # Dimensions = 1 * age
  W <- weight(L, WA, WB)

  # Maturity at age
  # Dimensions = 1 * age
  Mat <- maturity(Rec_age, Max_age, K_mat, L, L50)

  # Selectivity at age (updated)
  # Dimensions = 1 * age
  S <- selectivity(Rec_age, Max_age, A1, L1, A2, L2, K, Fleets, A50_up,
                   A50_down, Alpha, F_fin, Beta, Cf)

  # Cutoff for maturity
  A50_mat <- ages[min(which(Mat > 0.5))]

  # Number of control rules
  CR <- length(Control_rules)

  # Range of natural mortalities (low, correct, and high) if error =/= 0
  if (M_Error != 0) {
    Nat_mortality <- c(M - M_Error, M, M + M_Error)
  } else { Nat_mortality <- M}
  NM <- length(Nat_mortality)

  # Initialize age-structured population size matrix
  # Dimensions = age * area * time * CR * M * FDR values (3)
  N <- array(rep(0, num*A*TimeT*CR*NM*FDR), c(num, A, TimeT, CR, NM, FDR))

  # Initialize spawning stock biomass array
  # Dimensions = area * time * cr * M * FDR values (3)
  SSB <- array(rep(0, A*TimeT*CR*NM*FDR), c(A, TimeT, CR, NM, FDR))

  # Initialize abundance arrays
  # Dimensions = area * time * CR * M * FDR values (3)
  dimension <- ifelse(Ind_sampled == 'all', 1, 2)
  Abundance <- array(rep(0, A*TimeT*CR*NM*FDR*dimension),
                     c(A, TimeT, CR, NM, FDR, dimension))

  # Initialize biomass array
  # Dimensions = area * time * CR * M * FDR values (3)
  Biomass <- array(rep(0, A*TimeT*CR*NM*FDR), c(A, TimeT, CR, NM, FDR))

  # Recruitment normal variable
  # Dimensions = area * timeT * CR * M * FDR values (3)
  if (Stochasticity == T) {
    NuR <- array(rnorm(A*TimeT*CR*NM*FDR, 0, Sigma_R), c(A, TimeT, CR, NM, FDR))
  } else if (Stochasticity == F) {
    NuR <- array(rep(0, A*TimeT*CR*NM*FDR), c(A, TimeT, CR, NM, FDR))
  }

  # Recruitment error
  # Dimensions = area * timeT * CR * M * FDR values (3)
  Eps <- epsilon(A, TimeT, CR, NM, FDR, NuR, Rho_R)

  # Unfished spawning stock biomass
  B0 <- R0 / Phi

  # Initialize count array
  # Dimensions = area * time * transects * 2 * CR * M * FDR values (3)
 Count <- array(rep(0, A*TimeT*Transects*dimension*CR*NM*FDR),
                   c(A, TimeT, Transects, dimension, CR, NM, FDR))

 if (Sampling_Error == TRUE) {

    # Calculate standard deviation of normal variable for sampling
    # Based on Babcock & MacCall (2011): Eq. (15)
    Sigma_S <- sqrt(log(1 + (SP / X)^2))

    # Sampling normal variable
    # Dimensions = area * timeT * CR * M * FDR values (3)
    if (Stochasticity == T) {
      NuS <- array(rnorm(A*TimeT*CR*NM*FDR*dimension, 0, Sigma_S),
                   c(A, TimeT, CR, NM, FDR, dimension))
    } else if (Stochasticity == F) {
      NuS <- array(rep(0, A*TimeT*CR*NM*FDR*dimension),
                   c(A, TimeT, CR, NM, FDR, dimension))
    }

    # Calculate delta - constant of proportionality
    # Based on Babcock & MacCall (2011): Eq. (13)
    Delta <- P / D

    # Calculate gamma
    # Based on Babcock & MacCall (2011): Eq. (16)
    Gamma <- X / D

  }

  # Initialize fishing mortality rate
  # Dimensions = age * area * time * CR * M * FDR values (3)
  FM <- array(rep(0, num*A*TimeT*CR*NM*FDR), c(num, A, TimeT, CR, NM, FDR))

  # Initialize fishing effort in each area
  # Dimensions = area * time * CR * M * FDR values (3)
  E <- array(rep(0, A*TimeT*CR*NM*FDR), c(A, TimeT, CR, NM, FDR))

  # Initialize catch-at-age matrix
  # Dimensions = age * area * time * CR * M * FDR values (3)
  Catch <- array(rep(0, num*A*TimeT*CR*NM*FDR), c(num, A, TimeT, CR, NM, FDR))

  # Initialize yield matrix
  # Dimensions = area * time * CR * M * FDR values (3)
  Yield <- array(rep(0, A*TimeT*CR*NM*FDR), c(A, TimeT, CR, NM, FDR))

  # set constant fishing effort for first 50 years
  if (Fishing == T) {

    # Initial fishing effort
    E[, 1:Time1, , , ] <- rep(1/A, A*CR*Time1*NM*FDR)

  }

  # Stable age distribution, derived from equilibrium conditions with Fb = 0
  # Dimensions age
  SAD <- stable_AD(Rec_age, Max_age, W, R0, Mat, H, B0, Sigma_R, Fb, S, M,
                   eq_time = 150, A50_mat, Stochasticity = FALSE, Rho_R,
                   Recruitment_mode, LDP)

  # initialize density ratio matrix
  # Dimensions = timeT * CR * FDR
  Density_ratio <- array(rep(0, TimeT*CR*FDR), c(TimeT, CR, FDR))

  # ENM value - the nm value that represents the 'true' population
  ENM <- ifelse(M_Error == 0, 1, 2)

  # Enter N, abundance, biomasses, and E for time = 1 to rec_age
  # Dimensions = age * area * time * CR
  for (t in 1:Rec_age) {
    for (cr in 1:CR) {
      for (fdr in 1:FDR) {
        for (nm in 1:NM) {
          N[, , t, cr, nm, fdr] <- array(rep(SAD, A), c(num, A))
          FM[, , t, cr, nm, fdr] <- f_mortality(t, cr, nm, fdr, FM, A,
                                                Fb, E, S)
          Biomass[, t, cr, nm, fdr] <- colSums(N[, , t, cr, nm, fdr] * W)
          SSB[, t, cr, nm, fdr] <- colSums(N[, , t, cr, nm, fdr]*W*Mat)
          E[, t, cr, nm, fdr] <- rep(0.2, A)
          Catch[, , t, cr, nm, fdr] <- catch(t, cr, nm, fdr, FM,
                                             Nat_mortality, N, A, Fb, E, Catch)
          Yield[, t, cr, nm, fdr] <- colSums(Catch[, , t, cr, nm, fdr]*W)
          Abundance[, t, cr, nm, fdr, 1] <- colSums(N[, , t, cr, nm, fdr])

          # Abundance
          if (Ind_sampled == 'mature' || is.null(Ind_sampled)) {
            Abundance[, t, cr, nm, fdr, 2] <- colSums(N[A50_mat:(Max_age-Rec_age + 1),
                                                        , t, cr, nm, fdr])}
          if (t > 1 && Sampling_Error == TRUE) {
            Count[, t, , , cr, nm, fdr] <- sampling(t, cr, nm, fdr, Delta,
                                                    Gamma, Abundance, Transects,
                                                    X, Count, NuS, A, Ind_sampled)

          }

        }

        Density_ratio[t, cr, fdr] <- true_DR(t, cr, fdr, Abundance, Inside,
                                             Outside, Density_ratio, ENM)
      }
    }
  }

  if (Sampling_Error == TRUE) {

    output <- list(Inside, Outside, FDR, TimeT, L, W, S, Mat, A50_mat, CR,
                   Nat_mortality, NM, N, SSB, Biomass, Eps, B0, FM, E, Catch,
                   Yield, Density_ratio, ENM, Abundance, Count, Sigma_S, NuS,
                   Delta, Gamma)
  } else { output <- list(Inside, Outside, FDR, TimeT, L, W, S, Mat, A50_mat,
                          CR, Nat_mortality, NM, N, SSB, Biomass, Eps, B0, FM,
                          E, Catch, Yield, Density_ratio, ENM, Abundance, Count)
  }

  return(output)

}
