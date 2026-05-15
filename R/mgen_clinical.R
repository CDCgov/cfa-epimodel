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
  inf_time <- get_attr(dat, "inf_time") # time of infection
  rec_time <- get_attr(dat, "rec_time") # time of recovery
  amr_m <- get_attr(dat, "amr_m") # macrolide resistance status, where 0 = susceptible and 1 = resistant
  amr_q <- get_attr(dat, "amr_q") # quinolone resistance status, where 0 = susceptible and 1 = resistant
  curr_tx <- get_attr(dat, "curr_tx") # current treatment status (NA = no treatment, 1 = doxycycline, 2 = moxifloxacin)
  tx_end_day <- get_attr(dat, "tx_end_day") # day of treatment evaluation (NA if not treated)

  # Get parameters
  # NEED TO TRACK RECOVERIES HERE AND IN RECOVERY MODULE TO ENSURE PROPER CLINICAL FLOW
  # IN BOTH PLACES check if is.flow epi value assigned, then add to it
  rec_state <- get_param(dat, "rec_state") # recovery state (e.g., "s" for susceptible)
  flow_name <- paste0("i", rec_state, "_flow") # name of epi variable to track flow into recovery state

  # Establish empty vectors for clinical flow logic and tracking of clinical outcomes
  # These will be updated at each time step regardless of scenario, assign default vals here
  # Intermediate vectors not used for tracking but necessary for clinical flow logic
  # do NOT need to be initialized here
  n_recovered_m <- 0 # n infected who recover from doxycycline tx, male
  n_recovered_f <- 0 # n infected who recover from doxycycline tx, female
  n_tx_doxy_m <- 0 # n infected who begin treatment with doxycycline, male
  n_tx_doxy_f <- 0 # n infected who begin treatment with doxycycline, female
  n_tx_moxi_m <- 0 # n infected who begin treatment with moxifloxacin, male
  n_tx_moxi_f <- 0 # n infected who begin treatment with moxifloxacin, female
  n_tx_az_m <- 0 # n infected who begin treatment with azithromycin, male
  n_tx_az_f <- 0 # n infected who begin treatment with azithromycin, female

  # get current number of recoveries, some may already naturally cleared infection
  n_recovered_m <- get_epi(
    dat,
    paste0(flow_name, "_m"),
    at,
    override.null.error = TRUE
  )
  n_recovered_f <- get_epi(
    dat,
    paste0(flow_name, "_f"),
    at,
    override.null.error = TRUE
  )
  if (is.NULL(n_recovered_m)) {
    n_recovered_m <- 0
  }
  if (is.NULL(n_recovered_f)) {
    n_recovered_f <- 0
  }

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

    # Step 1: New Patients Seek Care and Get Treated with Doxycycline ---------
    ## Get ids of symptomatic infected individuals eligible for
    ## first-line treatment (doxycycline)
    ids_doxy_elig_m <- which(
      active == 1 &
        female == 0 &
        status == "i" &
        sympt == 1 &
        curr_tx == 0 &
        is.na(tx_end_day)
    )
    ids_doxy_elig_f <- which(
      active == 1 &
        female == 1 &
        status == "i" &
        sympt == 1 &
        curr_tx == 0 &
        is.na(tx_end_day)
    )
    ## Determine which eligible individuals seek care and receive doxycycline treatment
    ids_tx_doxy_m <- get_successful_ids_binom(ids_doxy_elig_m, p_seek_care_m)
    ids_tx_doxy_f <- get_successful_ids_binom(ids_doxy_elig_f, p_seek_care_f)
    all_tx_doxy_ids <- c(ids_tx_doxy_m, ids_tx_doxy_f)

    ## Update treatment status and day of treatment evaluation for those treated with doxycycline
    curr_tx[all_tx_doxy_ids] <- 1 # 1 = doxycycline
    tx_end_day[all_tx_doxy_ids] <- at + duration_doxy_tx

    ## Update tracker for individuals who begin doxycycline treatment
    n_tx_doxy_m <- length(ids_tx_doxy_m)
    n_tx_doxy_f <- length(ids_tx_doxy_f)

    # Step 2 - Sucess/Failure of Doxycycline Treatment ---------
    ## Get ids of those treated with doxycycline who finish tx course today
    ids_doxy_eval <- which(
      active == 1 &
        status == "i" &
        sympt == 1 &
        curr_tx == 1 &
        tx_end_day == at
    )
    ids_doxy_failure <- get_successful_ids_binom(ids_doxy_eval, p_doxy_failure)
    ids_doxy_success <- setdiff(ids_doxy_eval, ids_doxy_failure)

    ## Clear infection for those with successful doxycycline treatment
    status[ids_doxy_success] <- rec_state
    sympt[ids_doxy_success] <- NA # reset symptom status upon recovery
    rec_time[ids_doxy_success] <- at
    inf_time[ids_doxy_success] <- NA
    curr_tx[ids_doxy_success] <- NA
    tx_end_day[ids_doxy_success] <- NA

    ## Add recovery counts to tracker
    n_recovered_m <- n_recovered_m + sum(female[ids_doxy_success] == 0)
    n_recovered_f <- n_recovered_f + sum(female[ids_doxy_success] == 1)

    # Step 3 - NAAT for Doxycycline Failures & Moxi tx---------
    ## Get ids of those with doxycycline failure who get NAAT test
    ids_naat_elig_m <- which(
      active == 1 &
        sympt == 1 &
        status == "i" &
        curr_tx == 1 &
        female == 0 &
        tx_end_day < at # allow for some delay in care-seeking after treatment failure
    )

    ids_naat_elig_f <- which(
      active == 1 &
        sympt == 1 &
        status == "i" &
        curr_tx == 1 &
        female == 1 &
        tx_end_day < at # allow for some delay in care-seeking after treatment failure
    )

    ## Determine which eligible ids return and have sucessful NAAT test
    ## This setup means NAAT failures remain infected and may return for another test in future
    ids_naat_pos_m <- get_successful_ids_binom(
      ids_naat_elig_m,
      p_seek_care_m * naat_sensitivity
    )
    ids_naat_pos_f <- get_successful_ids_binom(
      ids_naat_elig_f,
      p_seek_care_f * naat_sensitivity
    )
    all_naat_pos_ids <- c(ids_naat_pos_m, ids_naat_pos_f)

    ## Update attributes & tracker value to reflect moxifloxacin treatment
    curr_tx[all_naat_pos_ids] <- 2 # 2 = moxifloxacin
    tx_end_day[all_naat_pos_ids] <- at + duration_moxi_tx

    ## Update tracker for individuals who begin moxifloxacin treatment
    n_tx_moxi_m <- length(ids_naat_pos_m)
    n_tx_moxi_f <- length(ids_naat_pos_f)

    # Step 4 - Sucess/Failure of Moxifloxacin Treatment ---------
    ## Get ids of those treated with moxifloxacin who finish tx course today
    ids_moxi_eval <- which(
      active == 1 &
        status == "i" &
        sympt == 1 &
        curr_tx == 2 &
        tx_end_day == at
    )
    ids_moxi_failure <- get_successful_ids_binom(ids_moxi_eval, p_moxi_failure)
    ids_moxi_success <- setdiff(ids_moxi_eval, ids_moxi_failure)

    ## Clear infection for those with successful moxifloxacin treatment
    status[ids_moxi_success] <- rec_state
    sympt[ids_moxi_success] <- NA # reset symptom status upon recovery
    rec_time[ids_moxi_success] <- at
    inf_time[ids_moxi_success] <- NA
    curr_tx[ids_moxi_success] <- NA
    tx_end_day[ids_moxi_success] <- NA

    ## Add recovery counts to tracker
    n_recovered_m <- n_recovered_m + sum(female[ids_moxi_success] == 0)
    n_recovered_f <- n_recovered_f + sum(female[ids_moxi_success] == 1)

    ## Update quinolone AMR status for those with moxifloxacin failure
    amr_q[ids_moxi_failure] <- 1 # set to resistant
    ## WHAT DO WE DO WITH THESE FOLKS????
  }
  # Update Epi Trackers
  dat <- set_epi(dat, "n_tx_doxy_m", at, n_tx_doxy_m)
  dat <- set_epi(dat, "n_tx_doxy_f", at, n_tx_doxy_f)
  dat <- set_epi(dat, "n_tx_moxi_m", at, n_tx_moxi_m)
  dat <- set_epi(dat, "n_tx_moxi_f", at, n_tx_moxi_f)
  dat <- set_epi(dat, "n_tx_az_m", at, n_tx_az_m)
  dat <- set_epi(dat, "n_tx_az_f", at, n_tx_az_f)
  dat <- set_epi(dat, "is_flow", at, sum(n_recovered_m, n_recovered_f))
  dat <- set_epi(dat, "is_flow_m", at, n_recovered_m)
  dat <- set_epi(dat, "is_flow_f", at, n_recovered_f)

  # Update Attributes
  dat <- set_attr(dat, "curr_tx", curr_tx)
  dat <- set_attr(dat, "tx_end_day", tx_end_day)
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "inf_time", inf_time)
  dat <- set_attr(dat, "rec_time", rec_time)
  dat <- set_attr(dat, "sympt", sympt)
  dat <- set_attr(dat, "amr_q", amr_q)
  dat <- set_attr(dat, "amr_m", amr_m)

  # Return dat object
  dat
}
