#' Base Model
#'
#' \code{base_model} runs one simulation of the density ratio base model, based
#'    on Babcock & MacCall (2011).
#'
#' @param Species character value, the species to analyze. Species names are in
#'    the format species abbreviation.state or stock.year assessed. Example:
#'    BR_CA_2003 is the California black rockfish stock assessed in 2003.
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
#'    'regional_DD' - larvae experience regional density dependence before
#'       settling evenly across all areas
#'    'local_DD' - larvae experience local density dependence before settling
#'       evely across all areas
#'    Default value is 'regional_DD'.
#' @param M_Error numeric value, the error between estimated and correct natural
#'    mortality.
#' @param Sampling_Var logical value, is there any variability in sampling?
#'    Default value is TRUE.
#' @param Recruitment_Var logical vector, is there any variability in
#'    recruitment? Default value is TRUE.
#' @param Surveys logical value, are surveys being conducted? Default value is
#'    TRUE.
#' @param Fishery_management logical value, is the fishery being managed using
#'    the control rules? Default value is TRUE.
#' @param Fishing logical value, is fishing occurring? Default value is TRUE.
#' @param Transects numerical value, the number of sampling transects conducted
#'    in each area to estimate density ratio. Default value is 24.
#' @param Adult_movement logical value, do adults move between areas? Default
#'    value is TRUE.
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
#' @param Floor_DR numeric value, the DR value under which effort will be
#'    reduced to 10\% of its starting value. Default value is 0.2.
#' @param Allocation character value, how effort is to be allocated. Values can
#'    be:
#'    'IFD' - the ideal free distribution. Effort is allocated proportional to
#'       the yield caught in each area in the previous timestep.
#'    'equal' - Effort is allocated equally between all areas.
#'    Default value is 'IFD'.
#' @param BM logical value, are the control rules from Babcock and MacCall 2011?
#'    Default value is FALSE.
#' @param LDP numeric value, the proportion of larvae that drift to adjacent
#'    areas. Default value is 0.1.
#' @param Output.N logical value, should the output include N? Default value is
#'    FALSE.
#' @param Output.Abundance logical value, should the output include Abundance?
#'    Default value is FALSE.
#' @param Output.Biomass logical value, should the output include biomass?
#'    Default value is FALSE.
#' @param Output.SSB logical value, should the output include spawning stock
#'    biomass? Default value is FALSE.
#' @param Output.Yield logical value, should the output include yield? Default
#'    value is FALSE.
#' @param Output.Effort logical value, should the output include effort? Default
#'    value is FALSE.
#' @param Output.Density.Ratio logical value, should the output include density
#'    ratios? Default value is TRUE.
#'
#' @return a list, containing numerical arrays of the relative biomass, yield,
#'    and spawning stock biomass, as well as the true density ratios over time
#'    and the transient target density ratios over time
#' @export
#'
#' @examples
#' base_model(Species = 'BR_OR_2015', R0 = 1e+5, A = 5, MPA = 3, Time1 = 50,
#'    Time2 = 20, Recruitment_mode = 'regional_DD', M_Error = 0.05,
#'    Sampling_Var = TRUE, Recruitment_Var = TRUE, Surveys = TRUE,
#'    Fishery_management = TRUE, Fishing = TRUE, Transects = 24,
#'    Adult_movement = TRUE, Final_DRs = c(0.6, 0.8), Years_sampled = 1,
#'    Areas_sampled = 'all', Ind_sampled = 'all', Floor_DR = 0.2,
#'    Allocation = 'IFD', BM = FALSE, LDP = 0.1, Output.N = FALSE,
#'    Output.Abundance = FALSE, Output.Biomass = FALSE, Output.SSB = FALSE,
#'    Output.Yield = FALSE, Output.Effort = FALSE, Output.Density.Ratio = TRUE)
base_model <- function(Species, R0 = 1e+5, A = 5, MPA = 3, Time1 = 50,
                       Time2 = 20, Recruitment_mode = 'regional_DD',
                       M_Error = 0.05, Sampling_Var = TRUE,
                       Recruitment_Var = TRUE, Surveys = TRUE,
                       Fishery_management = TRUE, Fishing = TRUE,
                       Transects = 24, Adult_movement = TRUE, Final_DRs,
                       Years_sampled = 1, Areas_sampled = 'all',
                       Ind_sampled = 'all', Floor_DR = 0.2, Allocation = 'IFD',
                       BM = FALSE, LDP = 0.1, Output.N = FALSE,
                       Output.Abundance = FALSE, Output.Biomass = FALSE,
                       Output.SSB = FALSE, Output.Yield = FALSE,
                       Output.Effort = FALSE, Output.Density.Ratio = TRUE) {

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
  if (!is.numeric(M_Error)) {stop('M_Error must be a numeric value.')}
  if (!is.logical(Sampling_Var)) {
    stop('Sampling_Var must be a logical value.')}
  if (!is.logical(Recruitment_Var)) {
    stop('Recruitment_Var must be a logical value.')}
  if (!is.logical(Surveys)) {stop('Surveys must be a logical value.')}
  if (!is.logical(Fishery_management)) {
    stop('Fishery_management must be a logical value.')}
  if (!is.logical(Fishing)) {stop('Fishing must be a logical value.')}
  if (Transects %% 1 != 0) {stop('Transects must be an integer value.')}
  if (!is.logical(Adult_movement)) {
    stop('Adult_movement must be a logical value.')}
  if (!is.numeric(Final_DRs)) {stop('Final_DRs must be a numeric vector.')}
  if (Years_sampled %% 1 != 0 && !is.null(Years_sampled)) {
    stop('Years_sampled must be an integer value or NULL.')}
  if (!is.character(Areas_sampled) && !is.null(Areas_sampled)) {
    stop('Areas_sampled must be a character value or NULL.')}
  if (!is.character(Ind_sampled) && !is.null(Ind_sampled)) {
    stop('Ind_sampled must be a character value or NULL.')}
  if (!is.numeric(Floor_DR)) {stop('Floor_DR must be a numeric vector.')}
  if (!is.character(Allocation)) {stop('Allocation must be a character value.')}
  if (!is.logical(BM)) {stop('BM must be a logical value.')}
  if (!is.numeric(LDP)) {stop('LDP must be a numeric value.')}
  if (!is.logical(Output.N)) {stop('Output.N must be a logical value.')}
  if (!is.logical(Output.Abundance)) {
    stop('Output.Abundance must be a logical value.')}
  if (!is.logical(Output.Biomass)) {
    stop('Output.Biomass must be a logical value.')}
  if (!is.logical(Output.SSB)) {stop('Output.SSB must be a logical value.')}
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
  if (Recruitment_mode != 'pool' && Recruitment_mode != 'closed' &&
      Recruitment_mode != 'regional_DD' && Recruitment_mode != 'local_DD') {
    stop('Recruitment_mode must be either "pool", "closed", "regional_DD", or
         "local_DD".')}
  if (M_Error < 0) {stop('M_Error must be greater than or equal to 0.')}
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
  if (Floor_DR <= 0) {stop('Floor_DR must be greater than 0.')}
  if (LDP < 0 || LDP > 1) {stop('LDP must be between 0 and 1.')}
  if (is.null(Years_sampled) && BM == FALSE) {
    stop('BM must be TRUE for Years_sampled to be NULL.')}
  if (is.null(Areas_sampled) && BM == FALSE) {
    stop('BM must be TRUE for Areas_sampled to be NULL.')}
  if (is.null(Ind_sampled) && BM == FALSE) {
    stop('BM must be TRUE for Ind_sampled to be NULL.')}

  # relationtional values
  if (MPA > A) {stop('MPA must be less than or equal to A.')}
  if (Floor_DR > min(Final_DRs)) {
    stop('Floor_DR must be less than or equal to the minimum final density
         ratio.')}

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
  H                      <- par[[13]]       # steepness
  Sigma_R                <- par[[14]]       # recruitment standard deviation
  Rho_R                  <- par[[15]]       # recruitment autocorrelation
  AMP                    <- par[[16]]       # adult movement proportion
  D                      <- par[[17]]       # depletion
  Fb                     <- par[[18]]       # fishing mortality to cause D
  P                      <- par[[19]]       # proportion of positive transects
  X                      <- par[[20]]       # mean of positive transects
  SP                     <- par[[21]]       # std of positive transects
  Fleets                 <- par[[22]]       # fishery fleet names
  Alpha                  <- par[[23]]       # slope for upcurve
  Beta                   <- par[[24]]       # slope for downcurve
  F_fin                  <- par[[25]]       # F_fin for fishery, 0 if asymptotic
  A50_up                 <- par[[26]]       # A50 for upcurve
  A50_down               <- par[[27]]       # A50 for downcurve
  Cf                     <- par[[28]]       # fraction of fishery caught / fleet

  ##### Population Dynamics - Non-Time Varying #################################

  # Initialize arrays for time-varying dynamics
  IA <- initialize_arrays(A, MPA, Final_DRs, Time1, Time2, R0, Rec_age, Max_age,
                          A1, L1, A2, L2, K, WA, WB, K_mat, Fb, L50, Sigma_R,
                          Rho_R, Fleets, Alpha, A50_up, A50_down, F_fin, Beta,
                          Cf, P, X, SP, M, Recruitment_Var, D, Transects, H,
                          Surveys, Fishing, M_Error, Sampling_Var,
                          Recruitment_mode, LDP, Ind_sampled, BM)

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
  N                <- IA[[12]]    # Population size, dim = age*area*time
  SSB              <- IA[[13]]    # Spawning stock biomass, dim = area*time
  Biomass          <- IA[[14]]    # Biomass, dim = area*time
  Eps              <- IA[[15]]    # Epsilon vector, dim = area*time*CR
  B0               <- IA[[16]]    # Unfished spawning stock biomass
  FM               <- IA[[17]]    # Fishing mortality rate, dim = age*area*time
  E                <- IA[[18]]    # nominal fishing effort in each area
  Catch            <- IA[[19]]    # Catch at age
  Yield            <- IA[[20]]    # Yield per area
  Density_ratio    <- IA[[21]]    # Density ratios
  Abundance        <- IA[[22]]    # Abundance of all and/or mature individuals
  Transects        <- IA[[23]]    # Species count when sampling, dim = area*time
  Count            <- IA[[24]]    # Species count when sampling, dim = area*time

  if (Sampling_Var == TRUE) {
    Sigma_S          <- IA[[25]]  # Sampling normal standard deviation
    NuS              <- IA[[26]]  # Sampling normal variable, dim = area*time*CR
    Delta            <- IA[[27]]  # Constant of proportionality
    Gamma            <- IA[[28]]  # Gamma
  }

  ##### Population Dynamics - Time Varying #####################################
  for (fdr in 1:FDR) {

    for (cr in 1:CR) {

      for (t in (Rec_age + 1):TimeT) {

        # effort allocation
        E[, t, cr, fdr] <- effort_allocation(t, cr, fdr, Allocation, E, Yield,
                                             Time1, Inside, Outside)

        # If there is adult movement, add movement
        if (Adult_movement == TRUE) {
          N[, , t, cr, fdr] <- movement(t, cr, fdr, N, A, AMP)
        }

        # Recruitment / larval movement (if applicable)
        R <- recruitment(t, cr, fdr, SSB, A, R0, H, B0, Eps, Sigma_R, Rec_age,
                         Recruitment_mode, LDP)

        # biology
        PD <- pop_dynamics(t, cr, fdr, Rec_age, Max_age, SSB, N, W, Mat, A, Fb,
                           E, S, FM, A50_mat, Biomass, Abundance, Fishing,
                           Nat_mortality, R, Ind_sampled)

        FM[, , t, cr, fdr]               <- PD[[1]]
        N[, , t, cr, fdr]                <- PD[[2]]
        Biomass[, t, cr, fdr]            <- PD[[3]]
        SSB[, t, cr, fdr]                <- PD[[4]]
        Abundance[, t, cr, fdr, ]        <- PD[[5]]

        # fishing
        if (Fishing == TRUE) {
          Catch[, , t, cr, fdr] <- catch(t, cr, fdr, FM, Nat_mortality, N, A,
                                         Fb, E, Catch)
          Yield[, t, cr, fdr] <- colSums(Catch[, , t, cr, fdr]*W)
        }

        # sampling
        if (Surveys == TRUE & Sampling_Var == TRUE) {
          Count[, t, , , cr, fdr] <- sampling(t, cr, fdr, Delta, Gamma,
                                              Abundance, Transects, X, Count,
                                              NuS, A, Ind_sampled)
        }

        # calculate true density ratio
        Density_ratio[t, cr, fdr] <- true_DR(t, cr, fdr, Abundance, Inside,
                                             Outside, Density_ratio)
        # management
        if (Fishery_management == TRUE && t >= Time1 && t < TimeT) {
          E[, t + 1, cr, fdr] <- control_rule(t, cr, fdr, A, E, Count, Time1,
                                              TimeT, Transects, Nat_mortality,
                                              Final_DRs, Inside, Outside,
                                              Years_sampled, Areas_sampled,
                                              Ind_sampled, Floor_DR, BM,
                                              Sampling_Var, Density_ratio,
                                              Abundance)
        }

      }

    }

  }

  ######## initialize output list ##############################################
  output <- list()

  if (Output.N == TRUE) { output$N <- N[, , Time1:TimeT, , ] }
  if (Output.Abundance == TRUE) {
    output$Abundance <- Abundance[, Time1:TimeT, , , ] }
  if (Output.Biomass == TRUE) {
    output$Biomass <- Biomass[, Time1:TimeT, , ] }
  if (Output.SSB == TRUE) { output$SSB <- SSB[, Time1:TimeT, , ] }
  if (Output.Yield == TRUE) {
    output$Yield <- colSums(Yield[, Time1:TimeT, , ]) }
  if (Output.Effort == TRUE) { output$Effort <- colSums(E[, , , ]) }
  if (Output.Density.Ratio == TRUE) {
    output$Density_ratio <- Density_ratio[Time1:TimeT, , ] }

  return(output)

}
