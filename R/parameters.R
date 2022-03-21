#' Life history parameters
#'
#' \code{parameters} loads the life history parameters for a specific fished
#'    stock.
#'
#' @param Species character value, the species to analyze. Species names are in
#'    the format species abbreviation_state or stock_year assessed. Example:
#'    BR.CA.2003 is the California black rockfish stock assessed in 2003.
#'
#' @return a list of numeric values and vectors required to run the base model
#'    and supporting functions.
#' @export
#'
#' @examples
#' parameters(Species = 'BR_CA_2003')

parameters = function(Species) {

  species_list <- c('BR_CA_2003', 'BR_OR_2015', 'CAB_CA_2005', 'CAB_OR_2019',
                    'LING_OW_2017', 'CR_OR_2015', 'China_OR_2015',
                    'BR_OR_2015_overfished', 'CAB_OR_2019_overfished',
                    'LING_OW_2017_overfished', 'CR_OR_2015_overfished')

  ###### Error handling ########################################################

  if (!is.character(Species)) {
    stop('Study species must be a character string.')
  }

  if (!is.element(Species, species_list)) {
    stop('Desired stock parameters are not known. Try a different species,
         region, and/or year.')
  }


  ##### Black Rockfish (OR) 2015 assessment, overfished #####
  # source: Cope et al. 2016
  if (Species == 'BR_OR_2015_overfished') { #####
    Max_age  <- 40                          # maximum age
    M        <- 0.17                        # natural mortality
    Rec_age  <- 3                           # age at recruitment
    WA       <- 2.6e-5;   WB <- 2.88        # weight at length parameters (f)
    A1       <- 1;        L1 <- 20.32       # growth parameters (f)
    A2       <- 40;       L2 <- 49.67
    K        <- 0.21
    L50      <- 43.69                       # length at 50% maturity
    K_mat    <- -0.66                       # slope of maturity curve
    H        <- 0.77                        # steepness
    Sigma_R  <- 0.5                         # recruitment standard deviation
    Rho_R    <- 0                           # recruitment autocorrelation
    AMP      <- 0.1                         # adult movement proportion
    D        <- 0.2                         # depletion
    Fb       <- 0.16                        # fishing mortality to cause D
    P        <- 0.77                        # Proportion of positive transects
                                            #       in PISCO monitoring data
    X        <- 15.42                       # mean of positive transects
    SP       <- 16.97                       # std of positive transects
    Fleets   <- c('trawl', 'live', 'dead',  # names of fleets
                  'ocean', 'shore')
    Alpha    <- c(0.325, 0.4, 0.35,
                  0.65, 0.425)              # slope of upcurve per fleet
    Beta     <- c(0.25, 0.5, 0.4, 1.1, 0.5) # slope of downcurve per fleet
    F_fin    <- c(0.325, 0.05, -0.11,
                  -0.025, 0.135)            # final selectivity if dome-shaped
    A50_up   <- c(7, 5, 5, 5, 3)            # A50 value for upcurve
    A50_down <- c(15, 13, 13, 12, 6)        # A50 value for downcurve
    Cf       <- c(0.0001, 0.1679, 0.0982,   # fraction of fishery caught / fleet
                  0.6979, 0.0358)
  }

  ##### Cabezon (OR) 2019 assessment, overfished #####
  # source: Cope et al. 2019
  if (Species == 'CAB_OR_2019_overfished') {
    Max_age  <- 20                          # maximum age
    M        <- 0.26                        # natural mortality
    Rec_age  <- 4                           # age at recruitment
    WA       <- 1.90e-5;  WB  <- 2.99       # weight at length parameters (f)
    A1       <- 4;        L1  <- 44.30      # growth parameters (f)
    A2       <- 20;       L2  <- 63.35
    K <- 0.225
    L50      <- 43                          # length at 50% maturity
    K_mat    <- -0.7                        # slope of maturity curve
    H        <- 0.7                         # steepness
    Sigma_R  <- 0.5                         # recruitment standard deviation
    Rho_R    <- 0                           # recruitment autocorrelation
    AMP      <- 0.1                         # adult movement proportion
    D        <- 0.44                        # depletion
    Fb       <- 0.24                        # fishing mortality to cause D
    P        <- 0.77                        # Proportion of positive transects
                                            #       in PISCO monitoring data
    X        <- 15.42                       # mean of positive transects
    SP       <- 16.97                       # std of positive transects
    Fleets   <- c('live', 'dead', 'ocean',  # names of fleets
                  'shore')
    Alpha    <- c(0.4, 0.33, 0.35, 0.9)     # slope of upcurve per fleet
    Beta     <- c(0.35, 0, 0, 0.2)          # slope of downcurve per fleet
    F_fin    <- c(0.7, 1, 1, 0.07)          # final select. if dome-shaped
    A50_up   <- c(3, 4, 2, 1)               # A50 value for upcurve
    A50_down <- c(17, 1, 1, 3)              # A50 value for downcurve
    Cf       <- c(0.6033, 0.0415, 0.3423,   # fraction of fishery
                  0.0130)
  }

  ##### Lingcod (OR and WA) 2017 assessment, overfished #####
  # source: Haltuch et al. 2018
  if (Species == 'LING_OW_2017_overfished') {
    Max_age  <- 25                          # maximum age
    M        <- 0.28                        # natural mortality
    Rec_age  <- 3                           # age at recruitment
    WA       <- 2.76e-6;  WB <- 3.28        # weight at length parameters (f)
    A1       <- 1;        L1 <- 17.28;      # growth parameters (f)
    A2       <- 20;       L2 <- 120;
    K        <- 0.128
    L50      <- 56.7                        # length at 50% maturity
    K_mat    <- -0.27                       # slope of maturity curve
    H        <- 0.7                         # steepness
    Sigma_R  <- 0.55                        # recruitment standard deviation
    Rho_R    <- 0                           # recruitment autocorrelation
    AMP      <- 0.1                         # adult movement proportion
    D        <- 0.37                        # depletion
    Fb       <- 0.16                        # fishing mortality to cause D
    P        <- 0.77                        # Proportion of positive transects
                                            #       in PISCO monitoring data
    X        <- 15.42                       # mean of positive transects
    SP       <- 16.97                       # std of positive transects
    Fleets   <- c('trawl', 'fixed_gear',    # names of fleets
                  'WArec', 'ORrec')
    Alpha    <- c(0.25, 0.25, 0.55, 1)      # slope of upcurve per fleet
    Beta     <- c(0.09, 0.3, 0.17, 0.15)    # slope of downcurve per fleet
    F_fin    <- c(0.07, 0, 0, 0)            # final select. if dome-shaped
    A50_up   <- c(3, 5, 5, 3)               # A50 value for upcurve
    A50_down <- c(15, 12, 10, 9)            # A50 value for downcurve
    Cf       <- c(0.2872, 0.1379, 0.3253,   # fraction of fishery
                  0.2496)
  }

  ##### Canary Rockfish (OR) 2015 assessment, overfished #####
  # source: Thorson & Wetzel 2015
  if (Species == 'CR_OR_2015') {
    Max_age  <- 84                          # maximum age
    M        <- 0.0521                      # natural mortality
    Rec_age  <- 3                           # age at recruitment
    WA       <- 1.18e-5;   WB   <- 3.094    # weight at length parameters (f)
    A1       <- 1;         L1   <- 9.05;    # growth parameters (f)
    A2       <- 20;        L2   <- 60.05;
    K        <- 0.129
    L50      <- 42                          # length at 50% maturity
    K_mat    <- -0.25                       # slope of maturity curve
    H        <- 0.773                       # steepness
    Sigma_R  <- 0.5                         # recruitment standard deviation
    Rho_R    <- 0                           # recruitment autocorrelation
    AMP      <- 0.1                         # adult movement proportion
    D        <- 0.4                         # depletion
    Fb       <- 0.04                        # fishing mortality to cause D
    P        <- 0.77                        # Proportion of positive transects
    #       in PISCO monitoring data
    X        <- 15.42                       # mean of positive transects
    SP       <- 16.97                       # std of positive transects
    Fleets   <- c('trawl', 'non-trawl',     # names of fleets
                  'rec', 'hake', 'research')
    Alpha    <- c(0.3, 0.6, 1, 1, 1)        # slope of upcurve per fleet
    Beta     <- c(1, 0, 1, 1, 0.08)         # slope of downcurve per fleet
    F_fin    <- c(0.36, 1, 0.175, 0.65, 0.8)# final selectivity if dome-shaped
    A50_up   <- c(5, 5, 4, 8, 1)            # A50 value for upcurve
    A50_down <- c(10, 50, 7, 11, 30)        # A50 value for downcurve
    Cf       <- c(0.3908, 0.3122, 0.2246,   # fraction of fishery caught / fleet
                  0.0295, 0.0429)           #       from upcurve to 1
  }

  ##### China Rockfsh (Northern OR) 2015 assessment #####
  # source: Cope et al. 2015
  if (Species == 'China_OR_2015') {
    Max_age  <- 79                          # maximum age
    M        <- 0.07                        # natural mortality
    Rec_age  <- 5                           # age at recruitment
    WA       <- 7.79E-06;   WB   <- 3.177   # weight at length parameters (f) #
    A1       <- 10;         L1   <- 30.5;   # growth parameters (f)
    A2       <- 30;         L2   <- 36.85;
    K        <- 0.147
    L50      <- 27                          # length at 50% maturity
    K_mat    <- -0.467                      # slope of maturity curve
    H        <- 0.773                       # steepness
    Sigma_R  <- 0.5                         # recruitment standard deviation
    Rho_R    <- 0                           # recruitment autocorrelation
    AMP      <- 0.2                         # adult movement proportion
    D        <- 0.6149                      # depletion
    Fb       <- 0.07                        # fishing mortality to cause D
    P        <- 0.77                        # Proportion of positive transects
    #       in PISCO monitoring data
    X        <- 15.42                       # mean of positive transects
    SP       <- 16.97                       # std of positive transects
    Fleets   <- c('Commercial', 'Rec_PC',
                  'Rec_PR')                 # names of fleets
    Alpha    <- c(1.05, 0.45, 0.45)         # slope of upcurve per fleet
    Beta     <- c(1, 1, 1)                  # slope of downcurve per fleet
    F_fin    <- c(1, 1, 1)                  # final selectivity if dome-shaped
    A50_up   <- c(10, 12, 12)               # A50 value for upcurve
    A50_down <- c(79, 79, 79)               # A50 value for downcurve
    Cf       <- c(0.02, 0.54, 0.44)         # fraction of fishery caught / fleet

  }

  ##### Black Rockfish (OR) 2015 assessment, actual #####
  # source: Cope et al. 2016
  if (Species == 'BR_OR_2015') { #####
    Max_age  <- 40                          # maximum age
    M        <- 0.17                        # natural mortality
    Rec_age  <- 3                           # age at recruitment
    WA       <- 2.6e-5;   WB <- 2.88        # weight at length parameters (f)
    A1       <- 1;        L1 <- 20.32       # growth parameters (f)
    A2       <- 40;       L2 <- 49.67
    K        <- 0.21
    L50      <- 43.69                       # length at 50% maturity
    K_mat    <- -0.66                       # slope of maturity curve
    H        <- 0.77                        # steepness
    Sigma_R  <- 0.5                         # recruitment standard deviation
    Rho_R    <- 0                           # recruitment autocorrelation
    AMP      <- 0.1                         # adult movement proportion
    D        <- 0.604                       # depletion
    Fb       <- 0.05                        # fishing mortality to cause
                                            #       overfishing (D = 0.)
    P        <- 0.77                        # Proportion of positive transects
                                            #       in PISCO monitoring data
    X        <- 15.42                       # mean of positive transects
    SP       <- 16.97                       # std of positive transects
    Fleets   <- c('trawl', 'live', 'dead',  # names of fleets
                  'ocean', 'shore')
    Alpha    <- c(0.325, 0.4, 0.35,
                  0.65, 0.425)              # slope of upcurve per fleet
    Beta     <- c(0.25, 0.5, 0.4, 1.1, 0.5) # slope of downcurve per fleet
    F_fin    <- c(0.325, 0.05, -0.11,
                  -0.025, 0.135)            # final selectivity if dome-shaped
    A50_up   <- c(7, 5, 5, 5, 3)            # A50 value for upcurve
    A50_down <- c(15, 13, 13, 12, 6)        # A50 value for downcurve
    Cf       <- c(0.0001, 0.1679, 0.0982,   # fraction of fishery caught / fleet
                  0.6979, 0.0358)
  }

  output <- list(Max_age, M, Rec_age, WA, WB, A1, L1, A2, L2, K, L50, K_mat, H,
                 Sigma_R, Rho_R, AMP, D, Fb, P, X, SP, Fleets, Alpha, Beta,
                 F_fin, A50_up, A50_down, Cf)

  return(output)

  ##### Cabezon (OR) 2019 assessment #####
  # source: Cope et al. 2019
  if (Species == 'CAB_OR_2019') {
    Max_age  <- 20                          # maximum age
    M        <- 0.26                        # natural mortality
    Rec_age  <- 4                           # age at recruitment
    WA       <- 1.90e-5;  WB  <- 2.99       # weight at length parameters (f)
    A1       <- 4;        L1  <- 44.30      # growth parameters (f)
    A2       <- 20;       L2  <- 63.35
    K <- 0.225
    L50      <- 43                          # length at 50% maturity
    K_mat    <- -0.7                        # slope of maturity curve
    H        <- 0.7                         # steepness
    Sigma_R  <- 0.5                         # recruitment standard deviation
    Rho_R    <- 0                           # recruitment autocorrelation
    AMP      <- 0.1                         # adult movement proportion
    D        <- 0.528                       # depletion
    Fb       <- 0.17                        # fishing mortality to cause D
    P        <- 0.77                        # Proportion of positive transects
    #       in PISCO monitoring data
    X        <- 15.42                       # mean of positive transects
    SP       <- 16.97                       # std of positive transects
    Fleets   <- c('live', 'dead', 'ocean',  # names of fleets
                  'shore')
    Alpha    <- c(0.4, 0.33, 0.35, 0.9)     # slope of upcurve per fleet
    Beta     <- c(0.35, 0, 0, 0.2)          # slope of downcurve per fleet
    F_fin    <- c(0.7, 1, 1, 0.07)          # final select. if dome-shaped
    A50_up   <- c(3, 4, 2, 1)               # A50 value for upcurve
    A50_down <- c(17, 1, 1, 3)              # A50 value for downcurve
    Cf       <- c(0.6033, 0.0415, 0.3423,   # fraction of fishery
                  0.0130)
  }


  ##### Lingcod (OR and WA) 2017 assessment #####
  # source: Haltuch et al. 2018
  if (Species == 'LING_OW_2017') {
    Max_age  <- 25                          # maximum age
    M        <- 0.28                        # natural mortality
    Rec_age  <- 3                           # age at recruitment
    WA       <- 2.76e-6;  WB <- 3.28        # weight at length parameters (f)
    A1       <- 1;        L1 <- 17.28;      # growth parameters (f)
    A2       <- 20;       L2 <- 120;
    K        <- 0.128
    L50      <- 56.7                        # length at 50% maturity
    K_mat    <- -0.27                       # slope of maturity curve
    H        <- 0.7                         # steepness
    Sigma_R  <- 0.55                        # recruitment standard deviation
    Rho_R    <- 0                           # recruitment autocorrelation
    AMP      <- 0.1                         # adult movement proportion
    D        <- 0.579                       # depletion
    Fb       <- 0.08                        # fishing mortality to cause D
    P        <- 0.77                        # Proportion of positive transects
    #       in PISCO monitoring data
    X        <- 15.42                       # mean of positive transects
    SP       <- 16.97                       # std of positive transects
    Fleets   <- c('trawl', 'fixed_gear',    # names of fleets
                  'WArec', 'ORrec')
    Alpha    <- c(0.25, 0.25, 0.55, 1)      # slope of upcurve per fleet
    Beta     <- c(0.09, 0.3, 0.17, 0.15)    # slope of downcurve per fleet
    F_fin    <- c(0.07, 0, 0, 0)            # final select. if dome-shaped
    A50_up   <- c(3, 5, 5, 3)               # A50 value for upcurve
    A50_down <- c(15, 12, 10, 9)            # A50 value for downcurve
    Cf       <- c(0.2872, 0.1379, 0.3253,   # fraction of fishery
                  0.2496)
  }


  ##### Canary Rockfish (OR) 2015 assessment #####
  # source: Thorson & Wetzel 2015
  if (Species == 'CR_OR_2015') {
    Max_age  <- 84                          # maximum age
    M        <- 0.0521                      # natural mortality
    Rec_age  <- 3                           # age at recruitment
    WA       <- 1.18e-5;   WB   <- 3.094    # weight at length parameters (f)
    A1       <- 1;         L1   <- 9.05;    # growth parameters (f)
    A2       <- 20;        L2   <- 60.05;
    K        <- 0.129
    L50      <- 42                          # length at 50% maturity
    K_mat    <- -0.25                       # slope of maturity curve
    H        <- 0.773                       # steepness
    Sigma_R  <- 0.5                         # recruitment standard deviation
    Rho_R    <- 0                           # recruitment autocorrelation
    AMP      <- 0.1                         # adult movement proportion
    D        <- 0.555                       # depletion
    Fb       <- 0.02                        # fishing mortality to cause D
    P        <- 0.77                        # Proportion of positive transects
    #       in PISCO monitoring data
    X        <- 15.42                       # mean of positive transects
    SP       <- 16.97                       # std of positive transects
    Fleets   <- c('trawl', 'non-trawl',     # names of fleets
                  'rec', 'hake', 'research')
    Alpha    <- c(0.3, 0.6, 1, 1, 1)        # slope of upcurve per fleet
    Beta     <- c(1, 0, 1, 1, 0.08)         # slope of downcurve per fleet
    F_fin    <- c(0.36, 1, 0.175, 0.65, 0.8)# final selectivity if dome-shaped
    A50_up   <- c(5, 5, 4, 8, 1)            # A50 value for upcurve
    A50_down <- c(10, 50, 7, 11, 30)        # A50 value for downcurve
    Cf       <- c(0.3908, 0.3122, 0.2246,   # fraction of fishery caught / fleet
                  0.0295, 0.0429)           #       from upcurve to 1
  }


}
