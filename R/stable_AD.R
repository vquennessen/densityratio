#' Stable Age Distribution
#'
#' \code{stableAD} produces the stable age distribution, or the proportional
#'    numbers at age from recruitment age to maximum age at equilibrium, after
#'    enough years of no fishing for proportional numbers at age to stabilize
#'
#' @param Rec_age numeric value, the age at recruitment, in years.
#' @param Max_age numeric value, the maximum age of fish or total lifespan, in
#'    years.
#' @param W numeric vector, the estimated weight at age from age at recruitment
#'    to maximum age, in kg.
#' @param R0 numeric value, set arbitrarily, the unfished recruitment.
#' @param Mat numeric vector, the estimated fraction of individuals mature at
#'    each age, from age at recruitment to maximum age, on the interval (0, 1).
#' @param H numeric value, the steepness of the stock-recruitment curve.
#' @param B0 numeric value, the unfished biomass, in kg.
#' @param Sigma_R numeric value, the recruitment standard deviation.
#' @param Fb numeric value, the historical fishing effort for the fished species.
#' @param S numeric vector, the selectivities at age from age at recruitment to
#'    maximum age, on the interval (0, 1).
#' @param M numeric value, the natural mortality on the interval (0, 1).
#' @param eq_time numeric value, the number of years to run the function to
#'    determine the stable age distribution.
#' @param A50_mat numeric value, the first age at which 50\% or more individuals
#'    are estimated to be mature, on the interval (Rec_age, Max_age).
#' @param Stochasticity logical vector, does recruitment contain a stochastic
#'    component? Default value is TRUE.
#' @param Rho_R numeric value, the recruitment autocorrelation on the interval
#'    (-1, 1). Default value is 0.
#' @param Nat_mortality numeric vector, the estimates of natural mortality.
#' @param Recruitment_mode character value, values can be:
#'    'closed' - the recruits in each area originate from adults in that area.
#'    'pool' - the recruits in each area come from a pool of larvae produced by
#'       adults in all areas.
#'    Default value is 'pool'.
#' @param A numeric value, the number of total areas in the model. Default
#'    value is 5.
#'
#' @return a numeric vector of numbers at age after enough years with no fishing
#'    that proportions remain stable over time.
#' @export
#'
#' @importFrom stats rnorm
#'
#' @examples
#' L <- length_age(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#'    L2 = 47.95, K = 0.2022, All_ages = FALSE)
#' W <- weight(L, WA = 1.68e-5, WB = 3)
#' Mat <- maturity(Rec_age = 2, Max_age = 35, K_mat = -0.4103, L, L50 = 39.53)
#' S <- selectivity(Rec_age = 2, Max_age = 35, A1 = 5, L1 = 32.21, A2 = 15,
#'    L2 = 47.95, K = 0.2022, Fleets = c('sport', 'hook', 'trawl'),
#'    A50_up = c(2, 5, 10), A50_down = c(6, 16, 35), Alpha = c(0.33, 0.6, 0.64),
#'    F_fin = c(0.25, 0.06, 1), Beta = c(1.2, 0.6, 0), Cf = c(0.71, 0.28, 0.01))
#' stable_AD(Rec_age = 2, Max_age = 35, W, R0 = 1e+5, Mat, H = 0.65, B0 = 1e+5/1.1,
#'    Sigma_R = 0.5, Fb = 0.2, S, M = 0.14, eq_time = 150, A50_mat = 8,
#'    Stochasticity = TRUE, Rho_R = 0, Nat_mortality = 0.14,
#'    Recruitment_mode = 'pool', A = 5)
stable_AD <- function(Rec_age, Max_age, W, R0, Mat, H, B0, Sigma_R, Fb, S, M,
                     eq_time, A50_mat, Stochasticity, Rho_R, Nat_mortality,
                     Recruitment_mode, A) {

  # range of ages
  ages <- Rec_age:Max_age
  num <- length(ages)

  # Initialize population size and fishing mortality arrays
  # Dimensions = age * 1 * time * 1 * 1
  N2 <- array(rep(0, num*eq_time), c(num, 1, eq_time, 1, 1))
  FM2 <- array(rep(0, num*eq_time), c(num, 1, eq_time, 1, 1))

  # Initialize effort, biomass, and SSB arrays
  # Dimensions = 1 * time * 1 * 1
  E2 <- array(rep(1, eq_time), c(1, eq_time, 1, 1))
  biomass2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
  SSB2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
  abundance_all2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
  abundance_mature2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))

  # Recruitment normal variable
  # Dimensions = area * timeT * CR * 1
  nuR2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))

  # Recruitment error
  # Dimensions = area * timeT * CR * 1
  Eps2 <- epsilon(A = 1, eq_time, CR = 1, NM = 1, nuR2, Rho_R)

  # Start each age class with 10 individuals
  # Enter FM, N, abundance, and biomasses for time = 1 to Rec_age
  # Dimensions = age * area * time * CR * M values (3)
  start_age <- A50_mat - Rec_age + 1
  for (t in 1:(Rec_age + 1)) {
    N2[, 1, t, 1, 1] <- rep(100, num)
    biomass2[1, t, 1, 1] <- sum(N2[, 1, t, 1, 1] * W)
    SSB2[1] <- sum(N2[, 1, t, 1, 1]*W*Mat)
    abundance_all2[1, t, 1, 1] <- sum(N2[, 1, t, 1, 1])
    abundance_mature2[1, t, 1, 1] <- sum(N2[start_age:num, 1, t, 1, 1])
  }

  # Step population forward in time with set fishing level
  for (t in (Rec_age + 1):(eq_time - 1)) {

    # biology
    PD <- pop_dynamics(a = 1, t, cr = 1, nm = 1, Rec_age, Max_age, SSB2,
                       N2, W, Mat, A, R0, H, B0, Eps2, Sigma_R, Fb, E2, S,
                       NM = 1, FM2, A50_mat, abundance_all2, abundance_mature2,
                       biomass2, Fishing = F, Nat_mortality = c(M),
                       Recruitment_mode)

    FM2[, 1, t, 1, 1]               <- PD[[1]]
    N2[, 1, t, 1, 1]                <- PD[[2]]
    abundance_all2[1, t, 1, 1]      <- PD[[3]]
    abundance_mature2[1, t, 1, 1]   <- PD[[4]]
    biomass2[1, t, 1, 1]            <- PD[[5]]
    SSB2[1, t, 1, 1]                <- PD[[6]]

  }

  # plotting for troubleshooting
  # plot(1:eq_time, N2[1, 1, 1:eq_time, 1], type = 'l', ylim = c(0, 2e4), col = 'green')
  # for (x in 2:(Num - 1)) {
  #   lines(1:eq_time, N2[x, 1, 1:eq_time, 1], col = 'red')
  # }
  # lines(1:eq_time, N2[Num, 1, 1:eq_time, 1], col = 'blue')

  SAD <- N2[, 1, eq_time - 1, 1, 1]

  return(SAD)

}
