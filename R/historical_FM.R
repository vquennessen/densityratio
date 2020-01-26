#' Historical fishing mortality to cause estimated depletion
#'
#' \code{historical_FM} gives the historical fishing mortality that is most
#'    likely to have resulted in the currently estimated depletion value for the
#'    fished stock.
#'
#' @param Species character value, the species to analyze. Species names are in
#'    the format species abbreviation.state or stock.year assessed. Example:
#'    BR.OR.2003 is the Oregon black rockfish stock assessed in 2003.
#' @param eq_time numeric value, the number of years to run the function to
#'    determine the historic fishing mortality. Default value is 150.
#' @param R0 numeric value, set arbitrarily, the unfished recruitment. Default
#'    value is 1e+5.
#' @param A numeric value, the number of total areas in the model. Default
#'    value is 5.
#' @param Stochasticity logical vector, does recruitment contain a stochastic
#'    component? Default value is FALSE.
#' @param Recruitment_mode character value, values can be:
#'    'closed' - the recruits in each area originate from adults in that area.
#'    'pool' - the recruits in each area come from a pool of larvae produced by
#'       adults in all areas.
#'    Default value is 'pool'.
#' @param Nat_mortality numeric vector, the estimates of natural mortality.
#'
#' @return a numeric value, the historical value of fishing mortality that would
#'    result in the estimated depletion value, on the interval (0, 1).
#' @export
#'
#' @examples
#' historical_FM(Species = 'BR.OR.2003', eq_time = 150, R0 = 1e+5, A = 5,
#'    Stochasticity = FALSE, Recruitment_mode = 'pool',
#'    Nat_mortality = c(0.14))
historical_FM <- function(Species, eq_time = 150, R0 = 1e+5, A = 5,
                          Stochasticity = FALSE, Recruitment_mode = 'pool',
                          Nat_mortality) {

##### load species parameters #####
par <- parameters(Species)

Max_age                <- par[[1]]        # maximum age
M                      <- par[[2]]        # natural mortality
Rec_age                <- par[[3]]        # age at recruitment
WA  <- par[[4]];  WB   <- par[[5]]        # weight at length parameters (f)
A1  <- par[[6]];  L1   <- par[[7]]        # growth parameters (f)
A2  <- par[[8]];  L2   <- par[[9]]
K   <- par[[10]]
L50                    <- par[[11]]       # length at 50% maturity
K_mat                  <- par[[12]]       # slope of maturity curve
H                      <- par[[14]]       # steepness
Phi                    <- par[[15]]       # unfished recruits per spawner
Sigma_R                <- par[[16]]       # recruitment standard deviation
Rho_R                  <- par[[17]]       # recruitment autocorrelation
                                          #       in PISCO monitoring data
D                      <- par[[19]]
SP                     <- par[[23]]       # std of positive transects

####### selectivity parameters #######
Fleets                 <- par[[26]]       # fishery fleet names
Alpha                  <- par[[27]]       # slope for upcurve
Beta                   <- par[[28]]       # slope for downcurve
F_fin                  <- par[[29]]       # F_fin for fishery, 0 if asymptotic
A50_up                 <- par[[30]]       # L50 for upcurve
A50_down               <- par[[31]]       # L50 for downcurve
Cf                     <- par[[32]]       # fraction of fishery caught / fleet

##### Calculate set values #####
ages <- Rec_age:Max_age                            # applicable ages
num <- length(ages)                                # number of age bins
# length at age
L   <- length_age(Rec_age, Max_age, A1, L1, A2, L2, K, All_ages = F)
# weight at age
W <- weight(L, WA, WB)
# fraction mature at age
Mat <- maturity(Rec_age, Max_age, K_mat, L, L50)
# age at 50% mature
A50_mat <- ages[min(which(Mat > 0.5))]
# unfished biomass
B0 <- R0/Phi
# selectivity at age
S <- selectivity(Rec_age, Max_age, L, Fleets, A50_up, A50_down, Alpha, F_fin,
                 Beta, Cf)
# NM - number of estimates of natural mortality
NM <- length(Nat_mortality)

# Recruitment error = 0 without stochasticity
Eps2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))

##### Initialize arrays #####

# Set historical FM to 0
Fb <- 0

# Fishing effort stays constant
E2 <- array(rep(1, eq_time), c(1, eq_time, 1, 1))

# Initialize FM and depletion levels
FM_values <- seq(from = 0, to = 1, by = 0.01)
fn <- length(FM_values)

# Initialize depletion vector
dep <- rep(0, fn)

# Initialize population size and catch arrays
# Dimensions = age * 1 * time * 1
N2 <- catch2 <- array(rep(0, num*eq_time), c(num, 1, eq_time, 1, 1))

# Initialize biomass, SSB, and recruitment error
# Dimensions = 1 * time * 1
SSB2 <- biomass2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
abundance_all2 <- abundance_mature2 <- array(rep(0, eq_time),
                                             c(1, eq_time, 1, 1))

rec_biomass <- array(rep(NA, fn*eq_time), c(fn, eq_time))

SAD <- stable_AD(Rec_age, Max_age, W, R0, Mat, H, B0, Sigma_R, Fb, S, M,
                 eq_time = 150, A50_mat, Stochasticity, Rho_R, Nat_mortality,
                 Recruitment_mode, A)

# Enter N, abundance, catch, and biomasses for time = 1 to rec_age
start_age <- A50_mat - Rec_age + 1
for (t in 1:Rec_age) {
  N2[, 1, t, 1, 1] <- SAD
  biomass2[1, t, 1, 1] <- sum(N2[, 1, t, 1, 1] * W)
  catch2[, 1, t, 1, 1] <- rep(0, num)
  SSB2[1, t, 1, 1] <- sum(N2[, 1, t - Rec_age, 1, 1]*W*Mat)
  abundance_all2[1, t, 1, 1] <- sum(N2[, 1, t, 1, 1])
  abundance_mature2[1, t, 1, 1] <- sum(N2[start_age:Max_age - 1, 1, t, 1, 1])
}

# Initialize FM matrix
FM2 <- array(rep(0, num*eq_time), c(num, 1, eq_time, 1, 1))

# Step population forward in time with set fishing level
for (t in (Rec_age + 1):eq_time) {

  PD <- pop_dynamics(a = 1, t, cr = 1, nm = 1, Rec_age, Max_age, SSB2,
                     N2, W, Mat, A = 1, R0, H, B0, Eps2, Sigma_R, Fb = 0, E2,
                     S, NM, FM2, A50_mat, abundance_all2, abundance_mature2,
                     biomass2, Fishing = T, Nat_mortality = M, Recruitment_mode)

  FM2[, 1, t, 1, 1]               <- rep(0, num)
  N2[, 1, t, 1, 1]                <- PD[[2]]
  abundance_all2[1, t, 1, 1]      <- PD[[3]]
  abundance_mature2[1, t, 1, 1]   <- PD[[4]]
  biomass2[1, t, 1, 1]            <- PD[[5]]
  SSB2[1, t, 1, 1]                <- PD[[6]]

}

# Calculate final biomass given zero fishing
FM0_biomass <- biomass2[1, eq_time, 1, 1]

# Substitute in values for Fb to get depletion level
for (i in 2:fn) {

  FM2 <- array(rep(FM_values[i], num*eq_time), c(num, 1, eq_time, 1, 1))

  # Step population forward in time with set fishing level
  for (t in (Rec_age + 1):eq_time) {

    PD <- pop_dynamics(a = 1, t, cr = 1, nm = 1, Rec_age, Max_age, SSB2,
                       N2, W, Mat, A = 1, R0, H, B0, Eps2, Sigma_R, Fb = 0, E2,
                       S, NM, FM2, A50_mat, abundance_all2, abundance_mature2,
                       biomass2, Fishing = T, Nat_mortality = M,
                       Recruitment_mode)

    FM2[, 1, t, 1, 1]               <- rep(FM_values[i], num)
    N2[, 1, t, 1, 1]                <- PD[[2]]
    abundance_all2[1, t, 1, 1]      <- PD[[3]]
    abundance_mature2[1, t, 1, 1]   <- PD[[4]]
    biomass2[1, t, 1, 1]            <- PD[[5]]
    SSB2[1, t, 1, 1]                <- PD[[6]]

  }

  dep[i] <- 1 - (biomass2[1, eq_time, 1, 1] / FM0_biomass)

}

closest_FM <- FM_values[which.min(abs(dep - D))]

graphics::plot(FM_values, dep)
graphics::abline(v = closest_FM, col = 'red')
graphics::abline(h = D, col = 'green')

return(closest_FM)

}