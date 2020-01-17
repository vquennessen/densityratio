#' Population Dynamics
#'
#' \code{pop_dynamics} returns updated values for fishing mortality, numbers at
#'    age, abundance of all and mature individuals, biomass, and spawning stock
#'    biomass in a specific area, at a specific time step, under a specific
#'    control rule, and with a specific estimate of natural mortality
#'
#' @param a temporary numeric value, the current area
#' @param t temporary numeric value, the current time step
#' @param cr temporary numeric value, the current control rule
#' @param nm temporary numeric value, the current natural mortality estimate
#' @param Rec_age numeric value, gives age at which fish has entered the fishery
#'    in years
#' @param Max_age numeric value, gives maximum age of fish or total lifespan in
#'    years
#' @param SSB numeric matrix, gives the spawning stock biomass of the whole
#'    stock for each area, at each timestep, under each control rule, and for
#'    each estimate of natural mortality
#' @param N numeric matrix, gives the number of individuals at each age, in each
#'    area, at each timestep, under each control rule, and for each estimate of
#'    natural mortality
#' @param W numeric vector, gives the estimated weight at age from age at
#'    recruitment to maximum age
#' @param Mat numeric vector, gives the estimated fraction of individuals mature
#'    at each age, from age at recruitment to maximum age
#' @param A numeric value, the number of total areas in the model. The default
#'    value is 5.
#' @param R0 numeric value, set arbitrarily, the unfished recruitment
#' @param H numeric value, the steepness of the stock-recruitment curve
#' @param B0 numeric value, set arbitrarily, the unfished biomass
#' @param Eps numeric matrix, recruitment error terms
#' @param Sigma_R numeric value, recruitment standard deviation
#' @param Fb numeric value, the historical fishing effort for the fished species
#' @param E numeric matrix, the relative fishing effort displayed in each area,
#'    at each time step, under each control rule, and for each natural mortality
#'    estimate
#' @param S numeric vector, the selectivities at age from age at recruitment to
#'    maximum age, on the interval (0, 1)
#' @param NM numeric value, the total number of estimated values of natural
#'    mortality. The default value is 3.
#' @param FM numeric matrix that corresponds to the fishing mortality at each
#'    age in each area, at each timestep, under all control rules, with all
#'    estimates of natural mortality
#' @param A50_mat numeric value, the first age at which 50\% or more individuals
#'    are estimated to be mature
#' @param Abundance_all numeric matrix, the total number of individuals in each
#'    area, at each timestep, under all control rules, with all estimates of
#'    natural mortality
#' @param Abundance_mature numeric matrix, the number of mature individuals in
#'    each area, at each timestep, under all control rules, with all estimates
#'    of natural mortality
#' @param Biomass numeric matrix, the total biomass in each area, at each time
#'    step, under all control rules, with all estimates of natural mortality,
#'    in kg
#' @param Fishing logical value, is fishing occurring? The default value is TRUE
#' @param Nat_mortality numeric vector, the estimates of natural mortality
#' @param Recruitment_mode character value, values can be 'closed' (if the
#'    recruits in each area originate from adults in that area) or 'pool' (if
#'    the recruits in each area come from a pool of larvae produced by all
#'    reproducing individuals in all areas) - the default value is 'pool'
#'
#' @return numeric matrix, with updated values of fishing mortality (FM),
#'    numbers at age (N), Abundance_all, Abundance_mature, Biomass, and spawning
#'    stock biomass (SSB)
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
#' S <- selectivity(Rec_age = 2, Max_age = 35, L,
#'    Fleets = c('sport', 'hook', 'trawl'), A50_up = c(2, 5, 10),
#'    A50_down = c(6, 16, 35), Alpha = c(0.33, 0.6, 0.64),
#'    F_fin = c(0.25, 0.06, 1), Beta = c(1.2, 0.6, 0), Cf = c(0.71, 0.28, 0.01))
#' pop_dynamics(a = 1, t = 3, cr = 1, nm = 2, Rec_age = 2, Max_age = 35, SSB, N,
#'    W, Mat, A = 5, R0 = 1e+5, H = 0.7, B0 = 1e+5/1.1, Eps, Sigma_R = 0.5,
#'    Fb = 0.2, E, S, NM = 3, FM, A50_mat = 8, Abundance_all, Abundance_mature,
#'    Biomass, Fishing = TRUE, Nat_mortality = c(0.09, 0.14, 0.19),
#'    Recruitment_mode = 'pool')
pop_dynamics <- function(a, t, cr, nm, Rec_age, Max_age, SSB, N, W, Mat, A = 5,
                         R0 = 1e+5, H, B0, Eps, Sigma_R, Fb, E, S, NM = 3, FM,
                         A50_mat, Abundance_all, Abundance_mature, Biomass,
                         Fishing = T, Nat_mortality, Recruitment_mode = 'pool') {

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
  N[num, a, t, cr, nm] <- N[num - 1, a, t - 1, cr, nm] * exp(-1 * (FM[num - 1, a, t - 1, cr, nm] +
                                                                     Nat_mortality[nm])) +
    N[num, a, t - 1, cr, nm] * exp(-1 * (FM[num, a, t - 1, cr, nm] + Nat_mortality[nm]))

  # Calculate abundance of all fish
  Abundance_all[a, t, cr, nm] <- sum(N[, a, t, cr, nm])

  # Calculate abundance of mature fish
  Abundance_mature[a, t, cr, nm] <- sum(N[A50_mat:Max_age - 1, a, t, cr, nm])

  # Calculate biomass of all fish
  Biomass[a, t, cr, nm] <- sum(N[, a, t, cr, nm] * W)

  # Calculate spawning stock biomass
  SSB[a, t, cr, nm] <- sum(N[, a, t, cr, nm]*W*Mat)

  output <- list(FM[, a, t, cr, nm], N[, a, t, cr, nm],
                 Abundance_all[a, t, cr, nm], Abundance_mature[a, t, cr, nm],
                 Biomass[a, t, cr, nm], SSB[a, t, cr, nm])

  return(output)

}
