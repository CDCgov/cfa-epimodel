#' @title Module for clinical aspects of M. genitalium infection
#'
#' @description This module handles the clinical aspects of M. genitalium infection,
#' including screening, testing, diagnosis, and treatment outcomes
#' depending on the specified scenario.
#'
#' @inheritParams vitals
#'
#' @export

mod_mgen_clinical <- function(dat, at) {
  # Get nodal attributes
  active <- get_attr(dat, "active")
  female <- get_attr(dat, "female") # biological sex, where 1 = female and 0 = male
  status <- get_attr(dat, "status") # infection status (s = susceptible, i = infected)
  sympt <- get_attr(dat, "sympt") # symptom status (1 = symptomatic, 0 = asymptomatic)
  amr_m <- get_attr(dat, "amr_m") # macrolide resistance status, where 0 = susceptible and 1 = resistant
  amr_q <- get_attr(dat, "amr_q") # quinolone resistance status, where 0 = susceptible and 1 = resistant
  curr_tx <- get_attr(dat, "curr_tx") # current treatment status (NA = no treatment, 1 = doxycycline, 2 = moxifloxacin)
  tx_end_day <- get_attr(dat, "tx_end_day") # day of treatment evaluation (NA if not treated)

  # Get parameters
  rec_state <- get_param(dat, "rec_state") # recovery state (e.g., "s" for susceptible)

  # Get indices of infected individuals by symptom status and sex
  ids_sympt_m <- which(active == 1 & status == "i" & sympt == 1 & female == 0) # indices of symptomatic infected, male
  ids_sympt_f <- which(active == 1 & status == "i" & sympt == 1 & female == 1) # indices of symptomatic infected, female
  ids_asympt_m <- which(active == 1 & status == "i" & sympt == 0 & female == 0) # indices of asymptomatic infected, male
  ids_asympt_f <- which(active == 1 & status == "i" & sympt == 0 & female == 1) # indices of asymptomatic infected, female

  # Establish empty vectors for clinical flow logic and tracking of clinical outcomes
  # Not all vectors will be used in all scenarios
  ids_naat_elig_m <- NULL # indices of infected eligible for NAAT testing, male
  ids_naat_elig_f <- NULL # indices of infected eligible for NAAT testing, female
  #ids_tested_m <- NULL # indices of infected who get tested, male
  #ids_tested_f <- NULL # indices of infected who get tested, female

  ids_tx_doxy_m <- NULL # indices of infected who get treated with doxycycline, male
  ids_tx_doxy_f <- NULL # indices of infected who get treated with doxycycline, female
  ids_tx_moxi_m <- NULL # indices of infected who get treated with moxifloxacin, male
  ids_tx_moxi_f <- NULL # indices of infected who get treated with moxifloxacin, female

  # Get scenario number for additional parameters and clinical flow logic
  scenario <- get_param(dat, "scenario")

  if (scenario == 1) {
    # Scenario 1: Baseline scenario with standard treatment
    # No asymptomatic screening, only symptomatic individuals seek care
    # Given course of doxycycline treatment
    # Most infections fail to clear and may persist with or without AMR
    # No resistance-guided therapy
    # NAAT test for M.gen only upon failure of doxycycline
    # If NAAT confirms M.gen infection, treat with moxifloxacin
    # Infections either clear or persist w/ AMR after moxifloxacin treatment
    # Men and women experience same clinical flow
    p_seek_care_m <- get_param(dat, "p_seek_care_m")
    p_seek_care_f <- get_param(dat, "p_seek_care_f")
    duration_doxy_tx <- get_param(dat, "duration_doxy_tx")
    duration_moxi_tx <- get_param(dat, "duration_moxi_tx")
    naat_sensitivity <- get_param(dat, "naat_sensitivity")
    p_doxy_failure <- get_param(dat, "p_doxy_failure")
    p_moxi_failure <- get_param(dat, "p_moxi_failure")

    # Get ids of symptomatic infected individuals eligible for first-line treatment (doxycycline)
    ids_tx_doxy_m <- get_successful_ids_binom(ids_sympt_m, p_seek_care_m)
    ids_tx_doxy_f <- get_successful_ids_binom(ids_sympt_f, p_seek_care_f)
    all_tx_doxy_ids <- c(ids_tx_doxy_m, ids_tx_doxy_f)

    # Update treatment status and day of treatment evaluation for those treated with doxycycline
    if (length(all_tx_doxy_ids) > 0) {
      curr_tx[all_tx_doxy_ids] <- 1 # 1 = doxycycline
      tx_end_day[all_tx_doxy_ids] <- at + duration_doxy_tx
    }

    # Get ids of those treated with doxycycline who fail treatment at eval day
    ids_doxy_eval <- which(
      active == 1 &
        sympt == 1 &
        status == "i" &
        curr_tx == 1 &
        tx_end_day == at
    )
    ids_doxy_failure <- get_successful_ids_binom(ids_doxy_eval, p_doxy_failure)
    ids_doxy_success <- setdiff(ids_doxy_eval, ids_doxy_failure)

    if (length(ids_doxy_success) > 0) {
      # Clear infection for those with successful doxycycline treatment
      status[ids_doxy_success] <- rec_state
      curr_tx[ids_doxy_success] <- NA
      tx_end_day[ids_doxy_success] <- NA
    }

    if (length(ids_doxy_failure) > 0) {
      # Get ids of those with doxycycline failure who get NAAT test
      ids_naat_elig_m <- which(
        active == 1 &
          sympt == 1 &
          status == "i" &
          curr_tx == 1 &
          tx_end_day == at
      )
    }
  }
  # Update Epi Trackers
  dat <- set_epi(dat, "n_tx_doxy_m", at, length(ids_tx_doxy_m))
  dat <- set_epi(dat, "n_tx_doxy_f", at, length(ids_tx_doxy_f))

  # Update Attributes
  dat <- set_attr(dat, "curr_tx", curr_tx)
  dat <- set_attr(dat, "tx_end_day", tx_end_day)

  # Return dat object
  dat
}
