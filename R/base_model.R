#' Base Model
#'
#' \code{base_model} runs one simulation of the density ratio base model, based
#'    on Babcock & MacCall (2011).
#'
#' @param Species character value, the species to analyze. Species names are in
#'    the format species abbreviation.state or stock.year assessed. Example:
#'    BR.CA.2003 is the California black rockfish stock assessed in 2003.
#' @param R0 numeric value, the unfished recruitment, set arbitrarily. Default
#'    value is 1e+5.
#' @param A numeric value, the number of total areas in the model. Default value
#'    is 5.
#' @param MPA numeric value, the area number to be designated as a marine
#'    reserve at Time1. Default value is 3.
#' @param Time1 numeric value, the number of years to run the model before a
#'    marine reserve is implemented. Default value is 50.
#' @param Time2 numeric value, the number of years to run the model after a
#'    marine reserve is implemented. Default value is 20.
#' @param Recruitment_mode character value, values can be:
#'    'closed' - the recruits in each area originate from adults in that area.
#'    'pool' - the recruits in each area come from a pool of larvae produced by
#'       adults in all areas.
#'    Default value is 'pool'.
#' @param Error numeric value, the error between estimated and correct natural
#'    mortality.
#' @param Stochasticity logical vector, does recruitment contain a stochastic
#'    component? Default value is TRUE.
#' @param Surveys logical value, are surveys being conducted? Default value is
#'    TRUE.
#' @param Fishery_management logical value, is the fishery being managed using
#'    the control rules? Default value is TRUE.
#' @param Fishing logical value, is fishing occurring? Default value is TRUE.
#' @param Transects numerical value, the number of sampling transects conducted
#'    in each area to estimate density ratio. Default value is 24.
#' @param Adult_movement logical value, do adults move between areas? Default
#'    value is TRUE.
#' @param Plotting logical value, should individual runs be plotted? Default
#'    value is FALSE.
#' @param Final_DRs numeric vector, the final target density ratios.
#' @param Years_sampled numeric value, the number of years of sampling upon
#'    which to base the estimate of density ratio. Default value is 1.
#' @param Areas_sampled character value, the areas to be sampled to calculate
#'    density ratio. Values can be:
#'    'all' - sample all areas.
#'    'far' - sample only the first and last areas (assuming the reserve is
#'       somewhere in the middle.
#'    Default value is 'all'.
#' @param Ind_sampled character value, the individuals to be sampled to
#'    calculate density ratio. Values can be:
#'    'all' - sample all individuals.
#'    'mature' - sample only mature individuals.
#'    Default value is 'all'.
#' @param Allocation character value, how effort is to be allocated. Values can
#'    be:
#'    'IFD' - the ideal free distribution. Effort is allocated proportional to
#'       the yield caught in each area in the previous timestep.
#'    'equal' - Effort is allocated equally between all areas.
#'    Default value is 'IFD'.
#' @param BM logical value, are the control rules from Babcock and MacCall 2011?
#'    Default value is FALSE.
#' @param Control_rules numeric vector, the control rules to be compared.
#'    Default value is c(1:6).
#' @param Output.FM logical value, should the output include FM? Default value
#'    is FALSE.
#' @param Output.N logical value, should the output include N? Default value is
#'    FALSE.
#' @param Output.Abundance.All logical value, should the output include
#'    Abundance_all? Default value is FALSE.
#' @param Output.Abundance.Mature logical value, should the output include
#'    Abundance_mature? Default value is FALSE.
#' @param Output.Biomass logical value, should the output include biomass?
#'    Default value is FALSE.
#' @param Output.SSB logical value, should the output include spawning stock
#'    biomass? Default value is FALSE.
#' @param Output.Catch logical value, should the output include catch? Default
#'    value is FALSE.
#' @param Output.Yield logical value, should the output include yield? Default
#'    value is FALSE.
#' @param Output.Effort logical value, should the output include effort? Default
#'    value is FALSE.
#' @param Output.Density.Ratio logical value, should the output include density
#'    ratios? Default value is FALSE.
#'
#' @return a list, containing numerical arrays of the relative biomass, yield,
#'    and spawning stock biomass, as well as the true density ratios over time
#'    and the transient target density ratios over time
#' @export
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics axis layout legend lines plot
#'
#' @examples
#' base_model(Species = 'BR_CA_2003', R0 = 1e+5, A = 5, MPA = 3, Time1 = 50,
#'    Time2 = 20, Recruitment_mode = 'pool', Error = 0.05, Stochasticity = TRUE,
#'    Surveys = TRUE, Fishery_management = TRUE, Fishing = TRUE, Transects = 24,
#'    Adult_movement = TRUE, Plotting = TRUE, Final_DRs = c(0.2, 0.4, 0.6, 0.8),
#'    Years_sampled = 1, Areas_sampled = 'all', Ind_sampled = 'all',
#'    Allocation = 'IFD', BM = FALSE, Control_rules = c(1:6), Output.FM = FALSE,
#'    Output.N = TRUE, Output.Abundance.All = FALSE,
#'    Output.Abundance.Mature = FALSE, Output.Biomass = TRUE, Output.SSB = TRUE,
#'    Output.Catch = FALSE, Output.Yield = TRUE, Output.Effort = TRUE,
#'    Output.Density.Ratio = TRUE)
base_model <- function(Species, R0 = 1e+5, A = 5, MPA = 3, Time1 = 50,
                       Time2 = 20, Recruitment_mode = 'pool', Error = 0.05,
                       Stochasticity = TRUE, Surveys = TRUE,
                       Fishery_management = TRUE, Fishing = TRUE,
                       Transects = 24, Adult_movement = TRUE, Plotting = FALSE,
                       Final_DRs, Years_sampled = 1, Areas_sampled = 'all',
                       Ind_sampled = 'all', Allocation = 'IFD', BM = FALSE,
                       Control_rules = c(1:6), Output.FM = FALSE,
                       Output.N = FALSE, Output.Abundance.All = FALSE,
                       Output.Abundance.Mature = FALSE, Output.Biomass = FALSE,
                       Output.SSB = FALSE, Output.Catch = FALSE,
                       Output.Yield = FALSE, Output.Effort = FALSE,
                       Output.Density.Ratio = FALSE) {

  ###### Error handling ########################################################

  # classes of variables
  if (!is.character(Species)) {
    stop('Study species must be a character string.')}
  if (R0 %% 1 != 0) {stop('R0 must be an integer value.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (MPA %% 1 != 0) {stop('MPA must be an integer value.')}
  if (Time1 %% 1 != 0) {stop('Time1 must be an integer value.')}
  if (Time2 %% 1 != 0) {stop('Time2 must be an integer value.')}
  if (!is.character(Recruitment_mode)) {
    stop('Recruitment mode must be a character value.')}
  if (!is.numeric(Error)) {stop('Error must be a numeric value.')}
  if (!is.logical(Stochasticity)) {
    stop('Stochasticity must be a logical value.')}
  if (!is.logical(Surveys)) {stop('Surveys must be a logical value.')}
  if (!is.logical(Fishery_management)) {
    stop('Fishery_management must be a logical value.')}
  if (!is.logical(Fishing)) {stop('Fishing must be a logical value.')}
  if (Transects %% 1 != 0) {stop('Transects must be an integer value.')}
  if (!is.logical(Adult_movement)) {
    stop('Adult_movement must be a logical value.')}
  if (!is.logical(Plotting)) {stop('Plotting must be a logical value.')}
  if (!is.numeric(Final_DRs)) {stop('Final_DRs must be a numeric vector.')}
  if (Years_sampled %% 1 != 0 && !is.null(Years_sampled)) {
    stop('Years_sampled must be an integer value or NULL.')}
  if (!is.character(Areas_sampled) && !is.null(Areas_sampled)) {
    stop('Areas_sampled must be a character value or NULL.')}
  if (!is.character(Ind_sampled) && !is.null(Ind_sampled)) {
    stop('Ind_sampled must be a character value or NULL.')}
  if (!is.character(Allocation)) {stop('Allocation must be a character value.')}
  if (!is.logical(BM)) {stop('BM must be a logical value.')}
  if (sum(Control_rules %% 1 != 0) != 0) {
    stop('Control_rules must be a vector of integers.')}
  if (!is.logical(Output.FM)) {stop('Output.FM must be a logical value.')}
  if (!is.logical(Output.N)) {stop('Output.N must be a logical value.')}
  if (!is.logical(Output.Abundance.All)) {
    stop('Output.Abundance.All must be a logical value.')}
  if (!is.logical(Output.Abundance.Mature)) {
    stop('Output.Abundance.Mature must be a logical value.')}
  if (!is.logical(Output.Biomass)) {
    stop('Output.Biomass must be a logical value.')}
  if (!is.logical(Output.SSB)) {stop('Output.SSB must be a logical value.')}
  if (!is.logical(Output.Catch)) {stop('Output.Catch must be a logical value.')}
  if (!is.logical(Output.Yield)) {stop('Output.Yield must be a logical value.')}
  if (!is.logical(Output.Effort)) {
    stop('Output.Effort must be a logical value.')}
  if (!is.logical(Output.Density.Ratio)) {
    stop('Output.Density.Ratio must be a logical value.')}

  # acceptable values
  if (R0 <= 0) {stop('R0 must be greater than 0.')}
  if (A <= 0) {stop('A must be greater than 0.')}
  if (MPA <= 0) {stop('MPA must be greater than 0.')}
  if (Time1 <= 0) {stop('Time1 must be greater than 0.')}
  if (Time2 <= 0) {stop('Time2 must be greater than 0.')}
  if (Recruitment_mode != 'pool' && Recruitment_mode != 'closed') {
    stop('Recruitment_mode must be either "pool" or "closed".')}
  if (Error < 0) {stop('Error must be greater than or equal to 0.')}
  if (Transects <= 0) {stop('Transects must be greater than 0.')}
  if (sum(Final_DRs <= 0) > 0) {
    stop('All values in Final_DRs must be greater than 0.')}
  if (Years_sampled <= 0 && !is.null(Years_sampled)) {
    stop('Years_sampled must be greater than 0 or NULL.')}
  if (is.numeric(Years_sampled) && Years_sampled <= 0) {
    stop('Years_sampled must be greater than 0 or NULL.')}
  if (is.character(Areas_sampled) && Areas_sampled != 'far' &&
      Areas_sampled != 'all' ) {
    stop('Areas_sampled must be either "far" or "all" or NULL.')}
  if (is.character(Ind_sampled) && Ind_sampled != 'mature' &&
      Ind_sampled != 'all') {
    stop('Ind_sampled must be either "mature" or "all" or NULL.')}
  if (sum(Control_rules <= 0) > 0) {
    stop('All values in Control_rules must be greater than 0.')}
  if (is.null(Years_sampled) && BM == FALSE) {
    stop('BM must be TRUE for Years_sampled to be NULL.')}
  if (is.null(Areas_sampled) && BM == FALSE) {
    stop('BM must be TRUE for Areas_sampled to be NULL.')}
  if (is.null(Ind_sampled) && BM == FALSE) {
    stop('BM must be TRUE for Ind_sampled to be NULL.')}

  # relationtional values
  if (MPA > A) {stop('MPA must be less than or equal to A.')}

  ##############################################################################

  ##### Load life history characteristics for species ##########################

  par <- parameters(Species)

  Max_age                <- par[[1]]        # maximum age
  M                      <- par[[2]]        # natural mortality
  Rec_age                <- par[[3]]        # age at recruitment
  WA <- par[[4]];  WB    <- par[[5]]        # weight at length parameters (f)
  A1 <- par[[6]];  L1    <- par[[7]]        # growth parameters (f)
  A2 <- par[[8]];  L2    <- par[[9]]
  K                      <- par[[10]]
  L50                    <- par[[11]]       # length at 50% maturity
  K_mat                  <- par[[12]]       # slope of maturity curve
  LDP                    <- par[[13]]       # larval drift proportion
  H                      <- par[[14]]       # steepness
  Phi                    <- par[[15]]       # unfished recruits per spawner
  Sigma_R                <- par[[16]]       # recruitment standard deviation
  Rho_R                  <- par[[17]]       # recruitment autocorrelation
  AMP                    <- par[[18]]       # adult movement proportion
  D                      <- par[[19]]       # depletion
  Fb                     <- par[[20]]       # fishing mortality to cause D
  P                      <- par[[21]]       # proportion of positive transects
  X                      <- par[[22]]       # mean of positive transects
  SP                     <- par[[23]]       # std of positive transects
  Fleets                 <- par[[24]]       # fishery fleet names
  Alpha                  <- par[[25]]       # slope for upcurve
  Beta                   <- par[[26]]       # slope for downcurve
  F_fin                  <- par[[27]]       # F_fin for fishery, 0 if asymptotic
  A50_up                 <- par[[28]]       # A50 for upcurve
  A50_down               <- par[[29]]       # A50 for downcurve
  Cf                     <- par[[30]]       # fraction of fishery caught / fleet

  ##### Population Dynamics - Non-Time Varying #################################

  # Initialize arrays for time-varying dynamics
  IA <- initialize_arrays(A, MPA, Final_DRs, Time1, Time2, R0, Rec_age, Max_age,
                          A1, L1, A2, L2, K, WA, WB, K_mat, Fb, L50, Sigma_R,
                          Rho_R, Fleets, Alpha, A50_up, A50_down, F_fin, Beta,
                          Cf, P, X, SP, M, Control_rules, Phi, Stochasticity, D,
                          Transects, H, Surveys, Fishing, Error,
                          Recruitment_mode)

  Inside           <- IA[[1]]     # Area(s) in the marine reserve
  Outside          <- IA[[2]]     # Areas not in the marine reserve
  FDR              <- IA[[3]]
  TimeT            <- IA[[4]]     # total amount of timesteps (years)
  L                <- IA[[5]]     # Length at age, dim = 1*age
  W                <- IA[[6]]     # Weight at age, dim = 1*age
  S                <- IA[[7]]     # Selectivity at age
  Mat              <- IA[[8]]     # Fraction mature at age, dim = 1*age
  A50_mat          <- IA[[9]]     # Age at which fraction mature > 0.5
  CR               <- IA[[10]]    # Number of control rules
  Nat_mortality    <- IA[[11]]    # Range of potential natural mortality values
  NM               <- IA[[12]]    # Number of potential natural mortality values
  N                <- IA[[13]]    # Population size, dim = age*area*time
  SSB              <- IA[[14]]    # Spawning stock biomass, dim = area*time
  Abundance_all    <- IA[[15]]    # Abundance, dim = area*time
  Abundance_mature <- IA[[16]]    # Abundance, dim = area*time
  Biomass          <- IA[[17]]    # Biomass, dim = area*time
  Eps              <- IA[[18]]    # Epsilon vector, dim = area*time*CR
  B0               <- IA[[19]]    # Unfished spawning stock biomass
  Count            <- IA[[20]]    # Species count when sampling, dim = area*time
  Sigma_S          <- IA[[21]]    # Sampling normal standard deviation
  NuS              <- IA[[22]]    # Sampling normal variable, dim = area*time*CR
  Delta            <- IA[[23]]    # Constant of proportionality
  Gamma            <- IA[[24]]    # Gamma
  FM               <- IA[[25]]    # Fishing mortality rate, dim = age*area*time
  E                <- IA[[26]]    # nominal fishing effort in each area
  Catch            <- IA[[27]]    # Catch at age
  Yield            <- IA[[28]]    # Yield per area
  Rel_biomass      <- IA[[29]]    # Relative biomass after reserve implementation
  Rel_yield        <- IA[[30]]    # Relative yield after reserve implementation
  Rel_SSB          <- IA[[31]]    # Relative SSB after reserve implementation
  Density_ratio    <- IA[[32]]    # Density ratios
  ENM              <- IA[[33]]    # nm value that represents "true" M

  ##### Population Dynamics - Time Varying #####################################

  for (fdr in 1:FDR) {

    for (t in (Rec_age + 1):(Time1 - 1)) {

      for (cr in 1:CR) {

        for (nm in 1:NM) {

          # effort allocation
          E <- effort_allocation(t, cr, nm, fdr, Allocation, E, Yield, Time1,
                                 Inside, Outside)

          # If there is adult movement, add movement
          if (Adult_movement == TRUE) {N <- movement(t, cr, nm, fdr, N, A, AMP)}

          # Recruitment / larval movement (if applicable)
          R <- recruitment(t, cr, nm, fdr, SSB, A, R0, H, B0, Eps, Sigma_R,
                           Rec_age, Recruitment_mode, LDP)

          for (a in 1:A) {

            # biology
            PD <- pop_dynamics(a, t, cr, nm, fdr, Rec_age, Max_age, SSB,
                               N, W, Mat, A, Fb, E, S, NM, FM, A50_mat,
                               Abundance_all, Abundance_mature, Biomass, Fishing,
                               Nat_mortality, R)

            FM[, a, t, cr, nm, fdr]               <- PD[[1]]
            N[, a, t, cr, nm, fdr]                <- PD[[2]]
            Abundance_all[a, t, cr, nm, fdr]      <- PD[[3]]
            Abundance_mature[a, t, cr, nm, fdr]   <- PD[[4]]
            Biomass[a, t, cr, nm, fdr]            <- PD[[5]]
            SSB[a, t, cr, nm, fdr]                <- PD[[6]]

            # sampling
            if (Surveys == TRUE) {
              Count[a, t, , , cr, nm, fdr] <- sampling(a, t, cr, nm, fdr, Delta,
                                                       Gamma, Abundance_all,
                                                       Abundance_mature,
                                                       Transects, X, Count, NuS,
                                                       A)
            }

            # fishing
            if (Fishing == TRUE) {
              Catch[, a, t, cr, nm, fdr] <- catch(a, t, cr, nm, fdr, FM,
                                                  Nat_mortality, N, A, Fb, E,
                                                  Catch)
              Yield[a, t, cr, nm, fdr] <- sum(Catch[, a, t, cr, nm, fdr]*W)
            }
          }
        }
      }
    }
  }

  ##### Implement Reserve, and apply control rules #############################

  for (fdr in 1:FDR) {

    for (t in Time1:TimeT) {

      for (cr in 1:CR) {

        for (nm in 1:NM) {

          # effort allocation
          E <- effort_allocation(t, cr, nm, fdr, Allocation, E, Yield, Time1,
                                 Inside, Outside)

          # If there is adult movement, add movement
          if (Adult_movement == TRUE) {N <- movement(t, cr, nm, fdr, N, A, AMP)}

          # Recruitment / larval movement (if applicable)
          R <- recruitment(t, cr, nm, fdr, SSB, A, R0, H, B0, Eps, Sigma_R,
                           Rec_age, Recruitment_mode, LDP)

          for (a in 1:A) {

            # biology
            PD <- pop_dynamics(a, t, cr, nm, fdr, Rec_age, Max_age, SSB,
                               N, W, Mat, A, Fb, E, S, NM, FM, A50_mat,
                               Abundance_all, Abundance_mature, Biomass, Fishing,
                               Nat_mortality, R)

            FM[, a, t, cr, nm, fdr]               <- PD[[1]]
            N[, a, t, cr, nm, fdr]                <- PD[[2]]
            Abundance_all[a, t, cr, nm, fdr]      <- PD[[3]]
            Abundance_mature[a, t, cr, nm, fdr]   <- PD[[4]]
            Biomass[a, t, cr, nm, fdr]            <- PD[[5]]
            SSB[a, t, cr, nm, fdr]                <- PD[[6]]

            # sampling
            if (Surveys == TRUE) {
              Count[a, t, , , cr, nm, fdr] <- sampling(a, t, cr, nm, fdr, Delta,
                                                       Gamma, Abundance_all,
                                                       Abundance_mature,
                                                       Transects, X, Count, NuS,
                                                       A)
            }

            # fishing
            if (Fishing == TRUE) {
              Catch[, a, t, cr, nm, fdr] <- catch(a, t, cr, nm, fdr, FM,
                                                  Nat_mortality, N, A, Fb, E,
                                                  Catch)
              Yield[a, t, cr, nm, fdr] <- sum(Catch[, a, t, cr, nm, fdr]*W)
            }
          }
        }

        # management
        if (Fishery_management == TRUE && t < TimeT) {
          E[, t, cr, , fdr] <- control_rule(t, cr, nm, fdr, A, E, Count, Time1,
                                            TimeT, Transects, Nat_mortality,
                                            Final_DRs, Inside, Outside,
                                            Areas_sampled, Ind_sampled,
                                            Years_sampled, BM)
        }

        # calculate true density ratio
        Density_ratio <- true_DR(t, cr, fdr, Abundance_all, Inside, Outside,
                                 Density_ratio, Time1, ENM)

      }
    }
  }

  # initialize output list
  output <- list()

  # add output depending on arguments passed to base_model.R
  if (Output.FM == TRUE) { output$FM <- FM[, 1:(MPA - 1), , , ENM, ] }
  if (Output.N == TRUE) { output$N <- N[, 1:MPA, , , ENM, ] }
  if (Output.Abundance.All == TRUE) {
    output$Abundance_all <- Abundance_all[1:MPA, , , ENM, ] }
  if (Output.Abundance.Mature == TRUE) {
    output$Abundance_mature <- Abundance_mature[1:MPA, , , ENM, ] }
  if (Output.Biomass == TRUE) { output$Biomass <- Biomass[1:MPA, , , ENM, ] }
  if (Output.SSB == TRUE) { output$SSB <- SSB[1:MPA, , , ENM, ] }
  if (Output.Catch == TRUE) { output$Catch <- Catch[, 1:MPA, , , ENM, ] }
  if (Output.Yield == TRUE) { output$Yield <- Yield[1:(MPA - 1), , , ENM, ] }
  if (Output.Effort == TRUE) { output$Effort <- E[1:(MPA - 1), , , ENM, ] }
  if (Output.Density.Ratio == TRUE) {output$Density_ratio <- Density_ratio }

  return(output)

}
