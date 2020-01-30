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
#' @param MPAs numeric vector, the area numbers to be designated as marine
#'    reserves at Time1. Default value is c(3).
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
#' @param Final_DR numeric value, the final target density ratio.
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
#' @param Control_rules numeric vector, the control rules to be compared.
#'    Default value is c(1:6).
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
#' base_model(Species = 'BR.CA.2003', R0 = 1e+5, A = 5, MPAs = c(3), Time1 = 50,
#'    Time2 = 20, Recruitment_mode = 'pool', Error = 0.05, Stochasticity = TRUE,
#'    Surveys = TRUE, Fishery_management = TRUE, Fishing = TRUE, Transects = 24,
#'    Adult_movement = TRUE, Plotting = TRUE, Final_DR = 0.6, Years_sampled = 1,
#'    Areas_sampled = 'all', Ind_sampled = 'all', Allocation = 'IFD',
#'    Control_rules = c(1:6))
base_model <- function(Species, R0 = 1e+5, A = 5, MPAs = c(3), Time1 = 50,
                       Time2 = 20, Recruitment_mode = 'pool', Error,
                       Stochasticity = TRUE, Surveys = TRUE,
                       Fishery_management = TRUE, Fishing = TRUE,
                       Transects = 24, Adult_movement = TRUE, Plotting = FALSE,
                       Final_DR, Years_sampled = 1, Areas_sampled = 'all',
                       Ind_sampled = 'all', Allocation = 'IFD',
                       Control_rules = c(1:6)) {

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
  R                      <- par[[21]]       # proportion of positive transects
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
  IA <- initialize_arrays(A, MPAs, Time1, Time2, R0, Rec_age, Max_age, A1, L1,
                          A2, L2, K, WA, WB, K_mat, Fb, L50, Sigma_R, Rho_R,
                          Fleets, Alpha, A50_up, A50_down, F_fin, Beta, Cf, R,
                          X, SP, M, Control_rules, Phi, Stochasticity, D,
                          Transects, H, Surveys, Fishing, Error,
                          Recruitment_mode)

  Inside           <- IA[[1]]     # Area(s) in the marine reserve
  Outside          <- IA[[2]]     # Areas not in the marine reserve
  TimeT            <- IA[[3]]     # total amount of timesteps (years)
  L                <- IA[[4]]     # Length at age, dim = 1*age
  W                <- IA[[5]]     # Weight at age, dim = 1*age
  S                <- IA[[6]]     # Selectivity at age
  Mat              <- IA[[7]]     # Fraction mature at age, dim = 1*age
  A50_mat          <- IA[[8]]     # Age at which fraction mature > 0.5
  CR               <- IA[[9]]     # Number of control rules
  Nat_mortality    <- IA[[10]]    # Range of potential natural mortality values
  NM               <- IA[[11]]    # Number of potential natural mortality values
  N                <- IA[[12]]    # Population size, dim = age*area*time
  SSB              <- IA[[13]]    # Spawning stock biomass, dim = area*time
  Abundance_all    <- IA[[14]]    # Abundance, dim = area*time
  Abundance_mature <- IA[[15]]    # Abundance, dim = area*time
  Biomass          <- IA[[16]]    # Biomass, dim = area*time
  Eps              <- IA[[17]]    # Epsilon vector, dim = area*time*CR
  B0               <- IA[[18]]    # Unfished spawning stock biomass
  Count            <- IA[[19]]    # Species count when sampling, dim = area*time
  Sigma_S          <- IA[[20]]    # Sampling normal standard deviation
  NuS              <- IA[[21]]    # Sampling normal variable, dim = area*time*CR
  Delta            <- IA[[22]]    # Constant of proportionality
  Gamma            <- IA[[23]]    # Gamma
  FM               <- IA[[24]]     # Fishing mortality rate, dim = age*area*time
  E                <- IA[[25]]     # nominal fishing effort in each area
  Catch            <- IA[[26]]    # Catch at age
  Yield            <- IA[[27]]    # Yield per area
  Rel_biomass      <- IA[[28]]    # Relative biomass after reserve implementation
  Rel_yield        <- IA[[29]]    # Relative yield after reserve implementation
  Rel_SSB          <- IA[[30]]    # Relative SSB after reserve implementation
  Density_ratio    <- IA[[31]]    # Density ratios

  ##### Population Dynamics - Time Varying #####################################

  for (t in 3:(Time1 - 1)) {

    for (cr in 1:CR) {

      for (nm in 1:NM) {

        # If there is adult movement, add movement
        if (Adult_movement == TRUE) { N <- movement(t, cr, nm, N, A, AMP) }

        for (a in 1:A) {

          # biology
          PD <- pop_dynamics(a, t, cr, nm, Rec_age, Max_age, SSB, N, W, Mat, A,
                             R0, H, B0, Eps, Sigma_R, Fb, E, S, NM, FM, A50_mat,
                             Abundance_all, Abundance_mature, Biomass, Fishing,
                             Nat_mortality, Recruitment_mode)

          FM[, a, t, cr, nm]               <- PD[[1]]
          N[, a, t, cr, nm]                <- PD[[2]]
          Abundance_all[a, t, cr, nm]      <- PD[[3]]
          Abundance_mature[a, t, cr, nm]   <- PD[[4]]
          Biomass[a, t, cr, nm]            <- PD[[5]]
          SSB[a, t, cr, nm]                <- PD[[6]]

          # sampling
          if (Surveys == TRUE) {
              Count[a, t, , , cr, nm] <- sampling(a, t, cr, nm, Delta, Gamma,
                                                  Abundance_all, Abundance_mature,
                                                  Transects, X, Count, NuS, A)
          }

          # fishing
          if (Fishing == TRUE) {
            Catch[, a, t, cr, nm] <- catch(a, t, cr, nm, FM, Nat_mortality, N,
                                           A, Fb, E, Catch)
            Yield[a, t, cr, nm] <- sum(Catch[, a, t, cr, nm]*W)
          }
        }
      }
    }
  }

  ##### Implement Reserve, and apply control rules #############################

  for (t in Time1:TimeT) {

    for (cr in 1:CR) {

      for (nm in 1:NM) {

        # effort allocation
        E <- effort_allocation(t, cr, nm, Allocation, E, Yield, Time1, Inside,
                               Outside)

        # If there is adult movement, add movement
        if (Adult_movement == TRUE) { N <- movement(t, cr, nm, N, A, AMP) }

        for (a in 1:A) {

          # biology
          PD <- pop_dynamics(a, t, cr, nm, Rec_age, Max_age, SSB, N, W, Mat, A,
                             R0, H, B0, Eps, Sigma_R, Fb, E, S, NM, FM, A50_mat,
                             Abundance_all, Abundance_mature, Biomass, Fishing,
                             Nat_mortality, Recruitment_mode)

          FM[, a, t, cr, nm]               <- PD[[1]]
          N[, a, t, cr, nm]                <- PD[[2]]
          Abundance_all[a, t, cr, nm]      <- PD[[3]]
          Abundance_mature[a, t, cr, nm]   <- PD[[4]]
          Biomass[a, t, cr, nm]            <- PD[[5]]
          SSB[a, t, cr, nm]                <- PD[[6]]

          # sampling
          if (Surveys == TRUE) {
            Count[a, t, , , cr, nm] <- sampling(a, t, cr, nm, Delta, Gamma,
                                            Abundance_all, Abundance_mature,
                                            Transects, X, Count, NuS, A)
          }

          # fishing
          if (Fishing == TRUE) {
            Catch[, a, t, cr, nm] <- catch(a, t, cr, nm, FM, Nat_mortality, N,
                                           A, Fb, E, Catch)
            Yield[a, t, cr, nm] <- sum(Catch[, a, t, cr, nm]*W)
          }

        }

      }

      # management
        if (Fishery_management == TRUE && t < TimeT) {
          E[, t, cr, ] <- control_rule(t, cr, A, E, Count, Time1, TimeT,
                                       Transects, Nat_mortality, Final_DR,
                                       Inside, Outside, Areas_sampled,
                                       Ind_sampled, Years_sampled)
        }

      # calculate true density ratio
      Density_ratio <- true_DR(t, cr, Abundance_all, Inside, Outside,
                               Density_ratio, Time1)

    }

  }

  # Troubleshooting plots
  # plot(1:TimeT, N[1, 1, 1:TimeT, 1], type = 'l', col = 'green', ylim = c(0, 2e4))
  # for (x in 2:(Num - 1)) {
  #   lines(1:TimeT, N[x, 1, 1:TimeT, 1], col = 'red')
  # }
  # lines(1:TimeT, N[Num, 1, 1:TimeT, 1], col = 'blue')

  ##### Calculate relative biomass, yield, and SSB ############################

  # calculate relative biomass since reserve implementation
  for (a in 1:A) {
    for (cr in 1:CR) {
      Rel_biomass[a, , cr] <- Biomass[a, Time1:TimeT, cr, 2]/Biomass[a, Time1, cr, 2]
    }
  }

  # calculate relative biomass since reserve implementation
  for (a in 1:A) {
    for (cr in 1:CR) {
      Rel_yield[a, , cr] <- Yield[a, Time1:TimeT, cr, 2]/Yield[a, Time1, cr, 2]
    }
  }

  # calculate relative biomass since reserve implementation
  for (a in 1:A) {
    for (cr in 1:CR) {
      Rel_SSB[a, , cr] <- SSB[a, Time1:TimeT, cr, 2]/SSB[a, Time1, cr, 2]
    }
  }

  ##### Plotting ######

  if (Plotting == T) {

    # use red-blue color palette
    palette <- colorRampPalette(c('red', 'blue'))
    color <- palette(CR)

    # set line types - solid for correct M, dashed for high M, dotted for low M
    line_type <- c(2, 1, 3, 2, 1, 3)

    # set layout matrix for all plots
    layout_m <- matrix(c(1, 3, 2, 3), nrow = 2, ncol = 2, byrow = T)

    # set legend title and text and position
    legend_title <- expression(bold('Control Rule'))
    legend_text  <- c("\n Static \n Low M", "\n Static \n Correct M",
                     "\n Static \n High M", "\n Transient \n Low M",
                     "\n Transient \n Correct M", "\n Transient \n High M")
    position <- 'left'

    # transient DR for population with correct M
    y_DR <- transient_DR(Time1, TimeT, Final_DR, Nat_mortality, nm = 2)

    ##### Plot relative biomass over time after reserve implementation #########

    # y-axis limits
    y1 <- 0
    y2 <- 4
    y_by <- (y2 - y1)/2

    # x-axis limits
    x1 <- 0
    x2 <- Time2
    x_by <- x2/4

    # DR y-axis limits
    y1_dr <- 0
    y2_dr <- 2
    by_dr <- y2_dr/2

    for (a in 1:3) {

      # set plotting layout
      layout(mat = layout_m,
             widths = c(2, 0.4),                 # Widths of the 2 columns
             heights = c(4, 2))                  # Heights of the 2 rows

      area <- ifelse(a < 2, 'far from', ifelse(a == 3, 'in', 'near'))
      title <- sprintf("Relative biomass: %s reserve", area)

      # plot the relative biomass
      par(mar = c(0.1, 4.5, 3.1, 0.1))
      plot(1, type = 'l',                        # make an empty line graph
           main = title,                         # title of plot
           ylab = 'Relative Biomass',            # axis labels
           xaxt = 'n',
           yaxt = 'n',                           # get rid of y-axis
           xlim = c(x1, x2),                     # set x-axis limits
           ylim = c(y1, y2),
           cex.lab = 1.5, cex.main = 1.5)

      # set specific y-axis
      ytick <- seq(y1, y2, by = y_by)            # set y axis tick marks
      axis(side = 2,                             # specify y axis
           at = ytick,                           # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal

      for (cr in 1:CR) {
        lines(x1:x2, Rel_biomass[a, , cr],
              col = color[cr],                   # use pre-defined color palette
              lwd = 2,                           # set line width
              lty = (cr %% 3) + 1)               # set line type
      }

      # add a gray dotted line at y = 1
      lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3)

      # plot the density ratios over time
      par(mar = c(4.1, 4.5, 3.1, 0.1))
      plot(1, type = 'l',                        # make an empty line graph
           main = 'Density Ratios Over Time',    # title of plot
           ylab = 'Density Ratio',               # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n',
           yaxt = 'n',                           # get rid of y-axis
           xlim = c(0, Time2),                   # set x-axis limits
           ylim = c(0, y2_dr),
           cex.lab = 1.5, cex.main = 1.5)

      # set specific y-axis
      dr_ytick <- seq(y1_dr, y2_dr, by_dr)       # set y axis tick marks
      axis(side = 2,                             # specify y axis
           at = dr_ytick,                        # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal

      # set specific x-axis
      xtick <- seq(x1, x2, by = x_by)            # set x axis tick marks
      axis(side = 1,                             # specify x axis
           at = xtick,                           # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal

      for (cr in 1:CR) {
        lines(x1:x2, Density_ratio[x1:x2 + 1, cr],
              col = color[cr],                   # use pre-defined color palette
              lwd = 2,                           # set line width
              lty = (cr %% 3) + 1)}              # set line type

      # add a gray dotted line at target_DR over time
      lines(0:Time2, y_DR, col = 'gray', lty = 3)

      # add a legend
      par(mar = c(0.1, 0.1, 0.1, 0.1))
      plot(1, type = 'n', axes = F, xlab = '', ylab = '')
      legend(x = position, inset = 0, horiz = F, # position
             col = color,                        # apply color palette
             lwd = 2,                            # apply line thicknesses
             lty = line_type,                    # apply line patterns
             title = legend_title,               # add legend title
             legend = legend_text,               # add legend labels
             seg.len = 3,                        # adjust length of lines
             cex = 1.1,                          # adjust legend text size
             bty = 'n')

    }

    ##### Plot relative yield over time after reserve implementation ###########

    # y-axis limits
    yy1 <- 0
    yy2 <- 3
    yy_by <- (yy2 - yy1)/2

    for (a in 1:2) {

      # set plotting layout
      layout(mat = layout_m,
             widths = c(2, 0.4),                 # Widths of the 2 columns
             heights = c(4, 2))                  # Heights of the 2 rows

      area <- ifelse(a == 1, 'far from', 'near')
      title <- sprintf("Relative yield: %s reserve", area)

      # plot the relative yield
      par(mar = c(0.1, 4.5, 3.1, 0.1))
      plot(1, type = 'l',                        # make an empty line graph
           main = title,                         # title of plot
           ylab = 'Relative Yield',              # axis labels
           xaxt = 'n',
           yaxt = 'n',                           # get rid of y-axis
           xlim = c(x1, x2),                     # set x-axis limits
           ylim = c(yy1, yy2),
           cex.lab = 1.5, cex.main = 1.5)

      # set specific y-axis
      yytick <- seq(yy1, yy2, by = yy_by)        # set yaxis tick marks
      axis(side = 2,                             # specify y axis
           at = yytick,                          # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal

      for (cr in 1:CR) {
        lines(x1:x2, Rel_yield[a, , cr],
              col = color[cr],                   # use pre-defined color palette
              lwd = 2,                           # set line width
              lty = (cr %% 3) + 1)               # set line type
      }

      # add a gray dotted line at y = 1
      lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3)

      # plot the density ratio over time
      par(mar = c(4.1, 4.5, 3.1, 0.1))
      plot(1, type = 'l',                        # make an empty line graph
           main = 'Density Ratios Over Time',    # title of plot
           ylab = 'Density Ratio',               # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n',
           yaxt = 'n',                           # get rid of y-axis
           xlim = c(0, Time2),                   # set x-axis limits
           ylim = c(0, y2_dr),
           cex.lab = 1.5, cex.main = 1.5)

      # set specific y-axis
      dr_ytick <- seq(y1_dr, y2_dr, by_dr)       # set y axis tick marks
      axis(side = 2,                             # specify y axis
           at = dr_ytick,                        # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal

      # set specific x-axis
      xtick <- seq(x1, x2, by = x_by)            # set x axis tick marks
      axis(side = 1,                             # specify x axis
           at = xtick,                           # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal

      for (cr in 1:CR) {
        lines(x1:x2, Density_ratio[x1:x2 + 1, cr],
              col = color[cr],                   # use pre-defined color palette
              lwd = 2,                           # set line width
              lty = (cr %% 3) + 1)}              # set line type

      # add a gray dotted line at target_DR over time
      lines(0:Time2, y_DR, col = 'gray', lty = 3)

      # add a legend
      par(mar = c(0.1, 0.1, 0.1, 0.1))
      plot(1, type = 'n', axes = F, xlab = '', ylab = '')
      legend(x = position, inset = 0, horiz = F, # position
             col = color,                        # apply color palette
             lwd = 2,                            # apply line thicknesses
             lty = line_type,                    # apply line patterns
             title = legend_title,               # add legend title
             legend = legend_text,               # add legend labels
             seg.len = 3,                        # adjust length of lines
             cex = 1.1,                          # adjust legend text size
             bty = 'n')

    }

    ###### Plot relative SSB over time after reserve implementation ##############

    # y-axis limits
    yyy1 <- 0
    yyy2 <- 2
    yyy_by <- (yyy2 - yyy1)/2

    for (a in 1:3) {

      # set plotting layout
      layout(mat = layout_m,
             widths = c(2, 0.4),                 # Widths of the 2 columns
             heights = c(4, 2))                  # Heights of the 2 rows

      area <- ifelse(a < 2, 'far from', ifelse(a == 3, 'in', 'near'))
      title <- sprintf("Relative biomass: %s reserve", area)

      # plot the relative SSB
      par(mar = c(0.1, 4.5, 3.1, 0.1))
      plot(1, type = 'l',                        # make an empty line graph
           main = title,                         # title of plot
           ylab = 'Relative SSB',                # axis labels
           xaxt = 'n',
           yaxt = 'n',                           # get rid of y-axis
           xlim = c(x1, x2),                     # set x-axis limits
           ylim = c(yyy1, yyy2),
           cex.lab = 1.5, cex.main = 1.5)

      # set specific y-axis
      yyytick <- seq(yyy1, yyy2, by = yyy_by)    # set yaxis tick marks
      axis(side = 2,                             # specify y axis
           at = yyytick,                         # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal

      for (cr in 1:CR) {
        lines(x1:x2, Rel_SSB[a, , cr],
              col = color[cr],                   # use pre-defined color palette
              lwd = 2,                           # set line width
              lty = (cr %% 3) + 1)               # set line type
      }

      # add a gray dotted line at y = 1
      lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3)

      # plot the density ratio over time
      par(mar = c(4.1, 4.5, 3.1, 0.1))
      plot(1, type = 'l',                        # make an empty line graph
           main = 'Density Ratios Over Time',    # title of plot
           ylab = 'Density Ratio',               # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n',
           yaxt = 'n',                           # get rid of y-axis
           xlim = c(0, Time2),                   # set x-axis limits
           ylim = c(0, y2_dr),
           cex.lab = 1.5, cex.main = 1.5
      )

      # set specific y-axis
      dr_ytick <- seq(y1_dr, y2_dr, by_dr)       # set y axis tick marks
      axis(side = 2,                             # specify y axis
           at = dr_ytick,                        # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal

      # set specific x-axis
      dr_xtick <- seq(0, Time2, by = Time2/4)    # set x axis tick marks
      axis(side = 1,                             # specify x axis
           at = dr_xtick,                        # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal

      for (cr in 1:CR) {
        lines(x1:x2, Density_ratio[x1:x2 + 1, cr],
              col = color[cr],                   # use pre-defined color palette
              lwd = 2,                           # set line width
              lty = (cr %% 3) + 1)}              # set line type

      # add a gray dotted line at target_DR over time
      lines(0:Time2, y_DR, col = 'gray', lty = 3)

      # add a legend
      par(mar = c(0.1, 0.1, 0.1, 0.1))
      plot(1, type = 'n', axes = F, xlab = '', ylab = '')
      legend(x = position, inset = 0, horiz = F, # position
             col = color,                        # apply color palette
             lwd = 2,                            # apply line thicknesses
             lty = line_type,                    # apply line patterns
             title = legend_title,               # add legend title
             legend = legend_text,               # add legend labels
             seg.len = 3,                        # adjust length of lines
             cex = 1.1,                          # adjust legend text size
             bty = 'n')

      }

  }

  output <- list(Rel_yield, Rel_biomass, Rel_SSB, Density_ratio)

  return(output)

}
