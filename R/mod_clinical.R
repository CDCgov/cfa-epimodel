#' @title Module for clinical aspects of M. genitalium infection
#'
#' @description This module handles the clinical aspects of M. genitalium infection,
#' including screening, testing, diagnosis, and treatment outcomes
#' depending on the specified scenario. If treatment is sucessful, patients are flagged
#' for recovery but tracking and attribute updating occurs in the recovery module.
#'
#' @inheritParams vitals
#'
#' @export

mod_clinical_mgen <- function(dat, at) {
  # Get nodal attributes
  # NOTE: these attributes are used to determine clinical eligibility at beginning of each time step
  # DO NOT modify these attributes directly within this module, use helper functions instead
  active <- get_attr(dat, "active")
  female <- get_attr(dat, "female") # biological sex, where 1 = female and 0 = male
  status <- get_attr(dat, "status") # infection status (s = susceptible, i = infected, e = exposed)
  sympt <- get_attr(dat, "sympt") # symptom status (1 = symptomatic, 0 = asymptomatic)
  amr_m <- get_attr(dat, "amr_m") # macrolide resistance status, where 0 = susceptible and 1 = resistant
  amr_q <- get_attr(dat, "amr_q") # quinolone resistance status, where 0 = susceptible and 1 = resistant
  seek_tx_day <- get_attr(dat, "seek_tx_day")
  curr_tx <- get_attr(dat, "curr_tx") # current treatment status (NA = no treatment, 1 = doxycycline, 2 = moxifloxacin, 3 = azithromycin, 4 = other)
  tx_end_day <- get_attr(dat, "tx_end_day") # day of treatment evaluation (NA if not treated)
  tx_success <- get_attr(dat, "tx_success") # treatment success status (NA = not treated, 1 = successful & will recover this time step)
  tx_return_day <- get_attr(dat, "tx_return_day") # day the patient is expected to return to clinic after treatment failure

  # Get scenario number for additional parameters and clinical flow logic
  scenario <- get_param(dat, "scenario")
  mean_days_to_clinic_visit <- get_param(
    dat,
    "mean_days_to_clinic_visit"
  )

  # Establish empty vectors for clinical flow logic and tracking of clinical outcomes
  # These will be updated at each time step regardless of scenario, assign default vals here
  # Intermediate vectors not used for tracking but necessary for clinical flow logic
  # do NOT need to be initialized here
  n_naat_tests <- 0
  n_amr_m_tests <- 0
  n_amr_q_tests <- 0
  n_tx_doxy_m <- 0 # n infected who begin treatment with doxycycline, male
  n_tx_doxy_f <- 0 # n infected who begin treatment with doxycycline, female
  n_tx_moxi_m <- 0 # n infected who begin treatment with moxifloxacin, male
  n_tx_moxi_f <- 0 # n infected who begin treatment with moxifloxacin, female
  n_tx_az_m <- 0 # n infected who begin treatment with azithromycin, male
  n_tx_az_f <- 0 # n infected who begin treatment with azithromycin, female
  n_tx_other_m <- 0 # n infected who begin treatment with other, male
  n_tx_other_f <- 0 # n infected who begin treatment with other, female

  if (scenario == 1) {
    # Scenario 1: Baseline scenario with standard treatment
    # No asymptomatic screening, only symptomatic individuals seek care
    # Given course of doxycycline treatment
    # Most infections fail to clear and may persist with or without AMR
    # No resistance-guided therapy
    # Assume that GC/CT test usually administered w/ concurrent doxy tx is negative
    # If doxy fails to treat MG, proceed to moxifloxacin
    # Infections either clear or persist w/ AMR after moxifloxacin treatment
    # Men and women experience same clinical flow
    p_doxy_success <- get_param(dat, "p_doxy_success")
    p_doxy_failure <- 1 - p_doxy_success
    p_moxi_success <- get_param(dat, "p_moxi_success")
    p_moxi_failure <- 1 - p_moxi_success
    p_thirdline_success <- get_param(dat, "p_thirdline_success")
    p_thirdline_failure <- 1 - p_thirdline_success

    # Step 1: New Patients Seek Care and Get Treated with Doxycycline ---------
    ## Get ids of symptomatic infected individuals eligible for
    ## first-line treatment (doxycycline)
    ids_doxy_elig <- which(
      active == 1 &
        sympt == 1 &
        status == "i" &
        is.na(curr_tx) &
        seek_tx_day == at
    )

    ## Update treatment status and day of treatment evaluation for those treated with doxycycline
    dat <- update_attrs_for_new_treatment(
      dat,
      at,
      ids_doxy_elig,
      "doxy"
    )

    ## Update tracker for individuals who begin doxycycline treatment
    ids_doxy_elig_m <- ids_doxy_elig[female[ids_doxy_elig] == 0]
    ids_doxy_elig_f <- ids_doxy_elig[female[ids_doxy_elig] == 1]
    n_tx_doxy_m <- length(ids_doxy_elig_m)
    n_tx_doxy_f <- length(ids_doxy_elig_f)

    # Step 2 - Success/Failure of Doxycycline Treatment ------------------------------
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

    ## Flag outcome of doxycycline treatment
    dat <- flag_ids_for_recovery(dat, at, ids_doxy_success)
    dat <- flag_ids_for_treatment_failure(dat, at, ids_doxy_failure)

    ## Calculate return day for those who failed doxycycline treatment
    dat <- set_event_time_rnorm(
      dat,
      at,
      ids_doxy_failure,
      "tx_return_day",
      mean_days_to_clinic_visit
    )

    # Step 3 - Moxifloxacin Treatment for Doxycycline Failures ---------
    ## Get ids of those with doxycycline failure returning to clinic
    ids_moxi_elig <- which(
      active == 1 &
        sympt == 1 &
        status == "i" &
        curr_tx == 1 & # doxy is current treatment
        tx_success == 0 & # failed doxy
        tx_return_day == at # only consider those who are due to return today
    )

    ## Update treatment status and day of treatment evaluation for those treated with moxifloxacin
    dat <- update_attrs_for_new_treatment(
      dat,
      at,
      ids_moxi_elig,
      "moxi"
    )

    ## Update tracker for individuals who begin moxifloxacin treatment
    ids_moxi_elig_m <- ids_moxi_elig[female[ids_moxi_elig] == 0]
    ids_moxi_elig_f <- ids_moxi_elig[female[ids_moxi_elig] == 1]
    n_tx_moxi_m <- length(ids_moxi_elig_m)
    n_tx_moxi_f <- length(ids_moxi_elig_f)

    # Step 4 - Sucess/Failure of Moxifloxacin Treatment -----------------------------
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

    ## Flag outcome of moxifloxacin treatment
    dat <- flag_ids_for_recovery(dat, at, ids_moxi_success)
    dat <- flag_ids_for_treatment_failure(dat, at, ids_moxi_failure)

    ## Calculate return day for those who failed moxifloxacin treatment
    dat <- set_event_time_rnorm(
      dat,
      at,
      ids_moxi_failure,
      "tx_return_day",
      mean_days_to_clinic_visit
    )

    ## Update quinolone AMR status for those with moxifloxacin failure
    dat <- update_amr_status(dat, at, ids_moxi_failure, "amr_q")

    # STEP 5 THIRD LINE TREATMENT, THEN LOSS TO FOLLOW UP IF THIRD LINE FAILS -----------------------------
    # STOPPED HERE
  }

  if (scenario == 2) {
    # Scenario 2: Resistance-guided therapy available
    # Symptomatic individuals seek care as in scenario 1
    # Given course of doxycycline treatment
    # Most infections fail to clear and may persist with or without AMR
    # NAAT test for M.gen and macrolide resistance upon failure of doxycycline
    # If NAAT confirms M.gen infection and macrolide susceptibility, treat with moxifloxacin
    # If NAAT confirms M.gen infection and macrolide resistance, treat with azithromycin
    # Infections either clear or persist w/ AMR after moxifloxacin or azithromycin treatment

    # Parameters specific to this scenario
    p_seek_care_m <- get_param(dat, "p_seek_care_m")
    p_seek_care_f <- get_param(dat, "p_seek_care_f")
    duration_doxy_tx <- get_param(dat, "duration_doxy_tx")
    duration_moxi_tx <- get_param(dat, "duration_moxi_tx")
    duration_az_tx <- get_param(dat, "duration_az_tx")
    naat_sensitivity <- get_param(dat, "naat_sensitivity")
    p_doxy_success <- get_param(dat, "p_doxy_success")

    # Step 1: New Patients Seek Care and Get Treated with Doxycycline ---------
    # Step 2: Sucess/Failure of Doxycycline Treatment ---------
    # Step 3: NAAT for Doxycycline Failures & Resistance-Guided Tx---------
    # Step 4: Sucess/Failure of Moxifloxacin or Azithromycin Treatment ---------
  }
  # Update Epi Trackers
  dat <- set_epi(dat, "n_naat_tests", at, n_naat_tests)
  dat <- set_epi(dat, "n_amr_m_tests", at, n_amr_m_tests)
  dat <- set_epi(dat, "n_amr_q_tests", at, n_amr_q_tests)
  dat <- set_epi(dat, "n_tx_doxy_m", at, n_tx_doxy_m)
  dat <- set_epi(dat, "n_tx_doxy_f", at, n_tx_doxy_f)
  dat <- set_epi(dat, "n_tx_moxi_m", at, n_tx_moxi_m)
  dat <- set_epi(dat, "n_tx_moxi_f", at, n_tx_moxi_f)
  dat <- set_epi(dat, "n_tx_az_m", at, n_tx_az_m)
  dat <- set_epi(dat, "n_tx_az_f", at, n_tx_az_f)

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

#' @title  Clinical Helper Functions for Treatment and AMR Management
#' @inheritParams vitals
#' @param ids A vector of individual IDs for whom the event or treatment is being set or updated.
#' @return The updated `dat` object with modified attributes based on the specified event or treatment.
#' @name clinical_helpers
#'
NULL

#' @rdname clinical_helpers
#' @description Sets the time for a specified event for a given set of individuals,
#' drawing from a normal distribution with the specified mean and standard deviation.
#' @param event_name The name of the event attribute to be updated for the specified individuals.
#' @param param1 The mean of the normal distribution from which the event times are drawn.
#' @param param2 The standard deviation of the normal distribution from which the event times are drawn. Defaults to 1.
#' @export
set_event_time_rnorm <- function(
  dat,
  at,
  ids,
  event_name,
  param1,
  param2 = NULL
) {
  if (length(ids) > 0) {
    event_times_from_now <- nnorm(
      length(ids),
      mean = param1,
      sd = ifelse(is.null(param2), 1, param2)
    ) # assuming a standard deviation of 1 day for return times unless otherwise specified
    event_attr <- get_attr(dat, event_name)
    event_attr[ids] <- at + event_times_from_now
    dat <- set_attr(dat, event_name, event_attr)
  }

  # Return
  dat
}

#' @rdname clinical_helpers
#' @description Flags the specified individuals as having recovered from treatment.
#' @param tx_success_attr The name of the attribute indicating treatment success. Defaults to `"tx_success"`.
#' @param tx_success_val The value to set for successful treatment. Defaults to `1`.
#' @export
flag_ids_for_recovery <- function(
  dat,
  at,
  ids,
  tx_success_attr = "tx_success",
  tx_success_val = 1
) {
  if (length(ids) > 0) {
    tx_success <- get_attr(dat, tx_success_attr)
    tx_success[ids] <- tx_success_val
    dat <- set_attr(dat, tx_success_attr, tx_success)
  }

  # Return
  dat
}

#' @rdname clinical_helpers
#' @description Flags the specified individuals as having experienced treatment failure.
#' @param tx_success_attr The name of the attribute indicating treatment success. Defaults to `"tx_success"`.
#' @param tx_success_val The value to set for treatment failure. Defaults to `0`.
#' @export
flag_ids_for_treatment_failure <- function(
  dat,
  at,
  ids,
  tx_success_attr = "tx_success",
  tx_success_val = 0
) {
  if (length(ids) > 0) {
    tx_success <- get_attr(dat, tx_success_attr)
    tx_success[ids] <- tx_success_val
    dat <- set_attr(dat, tx_success_attr, tx_success)
  }

  # Return
  dat
}

#' @rdname clinical_helpers
#' @description Updates the antimicrobial resistance (AMR) status for the specified individuals.
#' @param amr_attr The name of the attribute indicating AMR status.
#' @param amr_present_val The value to set for individuals with AMR. Defaults to `1`.
#' @export
update_amr_status <- function(dat, at, ids, amr_attr, amr_present_val = 1) {
  if (length(ids) > 0) {
    amr_status <- get_attr(dat, amr_attr)
    amr_status[ids] <- amr_present_val
    dat <- set_attr(dat, amr_attr, amr_status)
  }

  # Return
  dat
}

#' @rdname clinical_helpers
#' @description Updates the attributes for individuals receiving a new treatment, including current treatment, treatment end day, treatment success, and return day.
#' @param treatment_type The type of treatment being administered. Should be one of `"doxy"`, `"moxi"`, `"az"`, or `"other"`.
#' @param curr_tx_attr The name of the attribute indicating the current treatment. Defaults to `"curr_tx"`.
#' @param tx_end_day_attr The name of the attribute indicating the treatment end day. Defaults to `"tx_end_day"`.
#' @param tx_success_attr The name of the attribute indicating treatment success. Defaults to `"tx_success"`.
#' @param tx_return_day_attr The name of the attribute indicating the return day for treatment. Defaults to `"tx_return_day"`.
#' @export
update_attrs_for_new_treatment <- function(
  dat,
  at,
  ids,
  treatment_type,
  curr_tx_attr = "curr_tx",
  tx_end_day_attr = "tx_end_day",
  tx_success_attr = "tx_success",
  tx_return_day_attr = "tx_return_day"
) {
  if (length(ids) > 0) {
    curr_tx <- get_attr(dat, curr_tx_attr)
    tx_end_day <- get_attr(dat, tx_end_day_attr)
    tx_success <- get_attr(dat, tx_success_attr)
    tx_return_day <- get_attr(dat, tx_return_day_attr)

    tx_duration <- get_param(dat, paste0("duration_", treatment_type, "_tx"))

    treatment <- switch(
      treatment_type,
      "doxy" = 1,
      "moxi" = 2,
      "az" = 3,
      "other" = 4
    )

    curr_tx[ids] <- treatment
    tx_end_day[ids] <- at + tx_duration
    tx_success[ids] <- NA
    tx_return_day[ids] <- NA

    dat <- set_attr(dat, curr_tx_attr, curr_tx)
    dat <- set_attr(dat, tx_end_day_attr, tx_end_day)
    dat <- set_attr(dat, tx_success_attr, tx_success)
    dat <- set_attr(dat, tx_return_day_attr, tx_return_day)
  }

  # Return
  dat
}
