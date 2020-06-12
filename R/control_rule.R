#' Control rule
#'
#' \code{control_rule} determines how to calculate the density ratio and manage
#'    the fishery (i.e. set fishing effort for the next timestep).
#'
#' @param t temporary numeric value, the current time step.
#' @param cr temporary numeric value, the current control rule.
#' @param fdr temporary numeric value, the current final target density ratio.
#' @param A numeric value, the number of total areas in the model. Default
#'    value is 5.
#' @param E numeric array, the relative fishing effort displayed in each area,
#'    at each time step, under each control rule, and for each natural mortality
#'    estimate.
#' @param Count numeric array, the number of individuals estimated to be in each
#'    area, at each timestep, under each control rule, for each estimate of
#'    natural mortality, for both all individuals and just mature individuals.
#' @param Time1 numeric value, the number of years to run the model before a
#'    marine reserve is implemented. Default value is 50.
#' @param TimeT numeric value, the number of years to run the model total.
#'    Default value is 70.
#' @param Transects numerical value, the number of sampling transects conducted
#'    in each area to estimate density ratio. Default value is 24.
#' @param Nat_mortality numeric vector, the estimates of natural mortality.
#' @param Final_DRs numeric vector, the final target density ratios.
#' @param Inside numeric vector, the area(s) inside the marine reserve. Default
#'    is c(3).
#' @param Outside numeric vector, the area(s) outside the marine reserve.
#'    Default is c(1, 2, 4, 5).
#' @param Years_sampled numeric value, the number of years of sampling upon
#'    which to base the estimate of density ratio. Default value is 1.
#' @param Areas_sampled character value, the areas to be sampled to calculate
#'    density ratio. Values can be:
#'    'all' - sample all areas.
#'    'far' - sample only areas farthest from the marine reserve.
#'    Default value is 'all'.
#' @param Ind_sampled character value, the individuals to be sampled to
#'    calculate density ratio. Values can be:
#'    'all' - sample all individuals.
#'    'mature' - sample only mature individuals.
#'    Default value is 'all'.
#' @param Floor_DR numeric value, the DR value under which effort will be
#'    reduced to 10\% of its starting value. Default value is 0.2.
#' @param BM logical value, are the control rules from Babcock and MacCall 2011?
#'    Default value is FALSE.
#' @param Sampling_Error logical value, is there any error in sampling? Default
#'    value is TRUE.
#' @param Density_ratio numeric array, the true density ratio at each time step
#'    under each control rule and for each final density ratio.
#' @return a numeric vector of fishing effort for the next timestep, under the
#'    specific control rule, with a specific estimate of natural mortality.
#' @param Abundance numeric array, the total number of individuals in each area,
#'    at each timestep, under all control rules, with all estimates of natural
#'    mortality.
#' @export
#'
#' @examples
#' A = 5; TimeT = 70; CR = 6; NM = 2; FDR = 4; Transects = 24
#' E <- array(rep(1, A*TimeT*CR*NM*FDR), c(A, TimeT, CR, NM, FDR))
#' Count <- array(rep(5, A*TimeT*Transects*2*CR*NM*FDR),
#'    c(A, TimeT, Transects, 2, CR, NM, FDR))
#' Density_Ratio <- array(rep(0.5, TimeT*CR*FDR), c(TimeT, CR, FDR))
#' Abundance <- array(rep(3400, A*TimeT*CR*NM*FDR*1),
#'    c(A, TimeT, CR, NM, FDR, 1))
#' control_rule(t = 51, cr = 1, fdr = 1, A = 5, E, Count, Time1 = 50,
#'    TimeT = 70, Transects = 24, Nat_mortality = c(0.14, 0.09, 0.19),
#'    Final_DRs = c(0.2, 0.4, 0.6, 0.8), Inside = 3, Outside = c(1, 2, 4, 5),
#'    Years_sampled = 1, Areas_sampled = 'all', Ind_sampled = 'all',
#'    Floor_DR = 0.2, BM = FALSE, Sampling_Error = TRUE, Density_Ratio,
#'    Abundance)
control_rule <- function(t, cr, fdr, A = 5, E, Count, Time1 = 50,
                         TimeT = 70, Transects = 24, Nat_mortality, Final_DRs,
                         Inside = 3, Outside = c(1, 2, 4, 5), Years_sampled = 1,
                         Areas_sampled = 'all', Ind_sampled = 'all',
                         Floor_DR = 0.2, BM = FALSE, Sampling_Error = TRUE,
                         Density_ratio, Abundance) {

  ###### Error handling ########################################################

  # classes of variables
  if (t %% 1 != 0) {stop('t must be an integer value.')}
  if (cr %% 1 != 0) {stop('cr must be an integer value.')}
  if (fdr %% 1 != 0) {stop('fdr must be an integer value.')}
  if (A %% 1 != 0) {stop('A must be an integer value.')}
  if (!is.numeric(E)) {stop('E must be a numeric array.')}
  if (!is.numeric(Count)) {stop('Count must be a numeric array.')}
  if (Time1 %% 1 != 0) {stop('Time1 must be an integer value.')}
  if (TimeT %% 1 != 0) {stop('TimeT must be an integer value.')}
  if (Transects %% 1 != 0) {stop('Transects must be an integer value.')}
  if (!is.numeric(Nat_mortality)) {
    stop('Nat_mortality must be a numeric vector.')}
  if (!is.numeric(Final_DRs)) {stop('Final_DRs must be a numeric vector.')}
  if (sum(Inside %% 1 != 0) != 0) {stop('Inside must be a vector of integers.')}
  if (sum(Outside %% 1 != 0) != 0) {
    stop('Outside must be a vector of integers.')}
  if (Years_sampled %% 1 != 0 && !is.null(Years_sampled)) {
    stop('Years_sampled must be an integer value or NULL.')}
  if (!is.character(Areas_sampled) && !is.null(Areas_sampled)) {
    stop('Areas_sampled must be a character value or NULL.')}
  if (!is.character(Ind_sampled) && !is.null(Ind_sampled)) {
    stop('Ind_sampled must be a character value or NULL.')}
  if (!is.numeric(Floor_DR)) {stop('Floor_DR must be a numeric vector.')}
  if (!is.logical(BM)) {stop('BM must be a logical value.')}
  if (!is.logical(Sampling_Error)) {
    stop('Sampling_Error must be a logical value.')}
  if (!is.numeric(Density_ratio)) {stop('Density_ratio must be a numeric array.')}
  if (!is.numeric(Abundance)) {stop('Abundance must be a numeric array.')}

  # acceptable values
  if (t <= 0) {stop('t must be greater than 0.')}
  if (cr <= 0) {stop('cr must be greater than 0.')}
  if (fdr <= 0) {stop('fdr must be greater than 0.')}
  if (A <= 0) {stop('A must be greater than 0.')}
  if (sum(E < 0) > 0) {
    stop('All values in E must be greater than or equal to 0.')}
  if (sum(Count < 0) > 0) {
    stop('All values in Count must be greater than or equal to 0.')}
  if (Time1 <= 0) {stop('Time1 must be greater than 0.')}
  if (TimeT <= 0) {stop('TimeT must be greater than 0.')}
  if (Transects <= 0) {stop('Transects must be greater than 0.')}
  if (sum(Nat_mortality <= 0) > 0 || sum(Nat_mortality > 1) > 0) {
    stop('All values in Nat_mortality must be between 0 and 1.')}
  if (sum(Final_DRs <= 0) > 0) {
    stop('All values in Final_DRs must be greater than 0.')}
  if (sum(Inside < 0) > 0) {
    stop('All values in Inside must be greater than or equal to 0.')}
  if (sum(Outside < 0) > 0) {
    stop('All values in Outside must be greater than or equal to 0.')}
  if (is.numeric(Years_sampled) && Years_sampled <= 0) {
    stop('Years_sampled must be greater than 0 or NULL.')}
  if (is.character(Areas_sampled) && Areas_sampled != 'far' &&
      Areas_sampled != 'all' ) {
    stop('Areas_sampled must be either "far" or "all" or NULL.')}
  if (is.character(Ind_sampled) && Ind_sampled != 'mature' &&
      Ind_sampled != 'all') {
    stop('Ind_sampled must be either "mature" or "all" or NULL.')}
  if (Floor_DR <= 0) {stop('Floor_DR must be greater than 0.')}
  if (sum(Density_ratio < 0) > 0) {
    stop('All values in Density_ratio must be greater than or equal to 0.')}
  if (sum(Abundance < 0) > 0) {
    stop('All values in Abundance must be greater than or equal to 0.')}

  # relational values
  if (sum(Inside > A) > 0) {
    stop('All values in Inside must be less than or equal to A.')}
  if (sum(Outside > A) > 0) {
    stop('All values in Outside must be less than or equal to A.')}
  if (sum(intersect(Inside, Outside)) > 0) {
    stop('Areas cannot both be inside and outside the marine reserve.')}
  if (Time1 >= TimeT) {stop('TimeT must be greater than Time1.')}
  if(dim(E)[1] != dim(Count)[1] | dim(E)[1] != A | dim(Abundance)[1] != A) {
    stop('E, Count, or Abundance has an incorrect number of areas.')}
  if(dim(E)[2] != dim(Count)[2] | dim(E)[2] != TimeT ||
     dim(Density_ratio)[1] != TimeT | dim(Abundance)[2] != TimeT) {
    stop('E, Count, Density_ratio, or Abundance has an incorrect number of time
         steps.')}
  if(dim(Count)[3] != Transects) {
    stop('Count has the wrong number of transects.')}
  if(dim(E)[3] != dim(Count)[5]  | dim(E)[3] != dim(Density_ratio)[2] |
     dim(Abundance)[3] != dim(E)[3]) {
    stop('E, Count, Density_ratio, or Abundance has an incorrect number of
         control rules.')}
  if(dim(E)[4] != dim(Count)[6] | dim(E)[4] != dim(Abundance)[4]) {
    stop('E, Count, Density_ratio, or Abundance has an incorrect number of
         natural mortality estimates.')}
  if(dim(E)[5] != dim(Count)[7] | dim(E)[5] != dim(Density_ratio)[3] |
     dim(E)[5] != dim(Abundance)[5]) {
    stop('E, Count, Density_ratio, or Abundance has an incorrect number of final
         density ratios.')}
  if (t > dim(E)[2]) {stop('The given "t" value is too high for E.')}
  if (cr > dim(E)[3]) {stop('The given "cr" value is too high for E.')}
  if (fdr > dim(E)[5]) {stop('The given "fdr" value is too high for E.')}
  if (Floor_DR > min(Final_DRs)) {
    stop('Floor_DR must be less than or equal to the minimum final density
         ratio.')}

  ##############################################################################

  # if there is no sampling error, use true DR to dictate effort for next time
  # step
  True_DR <- Density_ratio[t, cr, fdr]

  if (BM == FALSE) {

    nm <- ifelse(cr < 3, 1, 2)
    j <- ceiling(cr / 2)

    if (Sampling_Error == TRUE) {
        DR <- density_ratio(t, cr, nm, fdr, A, Count, Years_sampled,
                          Areas_sampled, Ind_sampled, Transects, Inside,
                          Outside)
      } else { DR <- True_DR }

    # static control rules, with constant target DRs (cr = 1, 3, 5)
    if (cr %% 2 == 1) {

      # calculate effort at the next timestep
      E[, t + 1, cr, , fdr] <- management(t, cr, fdr, E, DR,
                                          target_DR = Final_DRs[fdr],
                                          floor_DR = Floor_DR,
                                          effort_inc_allowed = 0.1, Time1)

      # transient control rules with shifting target DRs (cr = 2, 4, 6)
    } else if (cr %% 2 == 0) {

      target <- transient_DR(Time1, TimeT, Final_DRs, Nat_mortality, nm = j,
                             fdr)

      # calculate effort at the next timestep
      E[, t + 1, cr, , fdr] <- management(t, cr, fdr, E, DR,
                                          target_DR = target[t - Time1 + 1],
                                          floor_DR = Floor_DR,
                                          effort_inc_allowed = 0.1, Time1)
    }

  } else if (BM == TRUE) {

    if (cr == 1) {

      # calculate effort at the next timestep
      E[, t + 1, cr, 1, fdr] <- 1.1*E[, t, cr, 1, fdr]

    } else if (cr == 2) {

      DR <- true_DR(t, cr, fdr = 1, Abundance, Inside, Outside, Density_ratio,
                    BM = TRUE, Years_sampled = 3, Areas_sampled = 'all',
                    Ind_sampled = 'all', A = 5)

        # DR <- density_ratio(t, cr, nm = 1, fdr = 1, A, Count, Years_sampled = 3,
        #                   Areas_sampled = 'all', Ind_sampled = 'all',
        #                   Transects, Inside, Outside)

      # calculate effort at the next timestep
      E[, t + 1, cr, 1, fdr] <- management(t, cr, fdr = 1, E, DR,
                                          target_DR = 0.6, floor_DR = 0.2,
                                          effort_inc_allowed = 0.1, Time1)

    } else if (cr == 3) {

      DR <- true_DR(t, cr, fdr = 1, Abundance, Inside, Outside, Density_ratio,
                    BM = TRUE, Years_sampled = 1, Areas_sampled = 'all',
                    Ind_sampled = 'all', A = 5)

        # DR <- density_ratio(t, cr, nm = 1, fdr = 1, A, Count, Years_sampled = 1,
        #                   Areas_sampled = 'all', Ind_sampled = 'all',
        #                   Transects, Inside, Outside)

      # calculate effort at the next timestep
      E[, t + 1, cr, 1, fdr] <- management(t, cr, fdr = 1, E, DR,
                                          target_DR = 0.6, floor_DR = 0.2,
                                          effort_inc_allowed = 0.1, Time1)

    } else if (cr == 4) {

      DR <- true_DR(t, cr, fdr = 1, Abundance, Inside, Outside, Density_ratio,
                    BM = TRUE, Years_sampled = 1, Areas_sampled = 'far',
                    Ind_sampled = 'all', A = 5)

        # DR <- density_ratio(t, cr, nm = 1, fdr = 1, A, Count, Years_sampled = 1,
        #                   Areas_sampled = 'far', Ind_sampled = 'all',
        #                   Transects, Inside, Outside)

      # calculate effort at the next timestep
      E[, t + 1, cr, 1, fdr] <- management(t, cr, fdr = 1, E, DR,
                                          target_DR = 0.6, floor_DR = 0.2,
                                          effort_inc_allowed = 0.1, Time1)

    } else if (cr == 5) {

      DR <- true_DR(t, cr, fdr = 1, Abundance, Inside, Outside, Density_ratio,
                    BM = TRUE, Years_sampled = 1, Areas_sampled = 'all',
                    Ind_sampled = 'mature', A = 5)

        # DR <- density_ratio(t, cr, nm = 1, fdr = 1, A, Count, Years_sampled = 1,
        #                   Areas_sampled = 'all', Ind_sampled = 'mature',
        #                   Transects, Inside, Outside)

      # calculate effort at the next timestep
      E[, t + 1, cr, 1, fdr] <- management(t, cr, fdr = 1, E, DR,
                                          target_DR = 0.6, floor_DR = 0.2,
                                          effort_inc_allowed = 0.1, Time1)

    } else if (cr == 6) {

      DR <- true_DR(t, cr, fdr = 1, Abundance, Inside, Outside, Density_ratio,
                    BM = TRUE, Years_sampled = 1, Areas_sampled = 'all',
                    Ind_sampled = 'all', A = 5)

        # DR <- density_ratio(t, cr, nm = 1, fdr = 1, A, Count, Years_sampled = 1,
        #                   Areas_sampled = 'all', Ind_sampled = 'all',
        #                   Transects, Inside, Outside)

      # calculate effort at the next timestep
      E[, t + 1, cr, 1, fdr] <- management(t, cr, fdr = 1, E, DR,
                                          target_DR = 0.8, floor_DR = 0.2,
                                          effort_inc_allowed = 0.1, Time1)

    } else if (cr == 7) {

      # calculate effort at the next timestep
      E[, t + 1, cr, 1, fdr] <- E[, t, cr, 1, fdr]

    } else if (cr == 8) {

      DR <- true_DR(t, cr, fdr = 1, Abundance, Inside, Outside, Density_ratio,
                    BM = TRUE, Years_sampled = 1, Areas_sampled = 'all',
                    Ind_sampled = 'all', A = 5)

        # DR <- density_ratio(t, cr, nm = 1, fdr = 1, A, Count, Years_sampled = 1,
        #                   Areas_sampled = 'all', Ind_sampled = 'all',
        #                   Transects, Inside, Outside)

      # calculate effort at the next timestep
      E[, t + 1, cr, 1, fdr] <- management(t, cr, fdr = 1, E, DR,
                                          target_DR = 0.6, floor_DR = 0.2,
                                          effort_inc_allowed = 0, Time1)
    }
  }

  return(E[, t + 1, cr, , fdr])

}
