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
  amr_m <- get_attr(dat, "amr_m") # macrolide resistance status
  amr_q <- get_attr(dat, "amr_q") # quinolone resistance status
  curr_tx <- get_attr(dat, "curr_tx") # current treatment status (NA = no treatment, 1 = doxycycline, 2 = moxifloxacin)
  tx_day_eval <- get_attr(dat, "tx_day_eval") # day of treatment evaluation (NA if not yet evaluated)
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

  ids_tx_doxi_m <- NULL # indices of infected who get treated with doxycycline, male
  ids_tx_doxi_f <- NULL # indices of infected who get treated with doxycycline, female
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

    get_successful_ids_binom <- function(elig_ids, eval_prob) {
      l_elig_ids <- length(elig_ids)
      l_eval_prob <- length(eval_prob)
      if (length(elig_ids) > 0) {
        if (l_eval_prob == 1 | l_eval_prob == l_elig_ids) {
          elig_ids[rbinom(length(elig_ids), 1, eval_prob) == 1]
        } else {
          stop(
            "Length of eval_prob must be either 1 or equal to the number of eligible IDs."
          )
        }
      } else {
        NULL
      }
    }

    # Get ids of symptomatic infected individuals eligible for first-line treatment (doxycycline)
    ids_tx_doxi_m <- get_successful_ids_binom(ids_sympt_m, p_seek_care_m)
    ids_tx_doxi_f <- get_successful_ids_binom(ids_sympt_f, p_seek_care_f)
    all_tx_doxi_ids <- c(ids_tx_doxi_m, ids_tx_doxi_f)

    # Update treatment status and day of treatment evaluation for those treated with doxycycline
    if (length(all_tx_doxi_ids) > 0) {
      curr_tx[all_tx_doxi_ids] <- 1 # 1 = doxycycline
      tx_day_eval[all_tx_doxi_ids] <- at + duration_doxy_tx
    }
  }

  # Update Epi Trackers
  dat <- set_epi(dat, "n_tx_doxi_m", length(ids_tx_doxi_m))
  dat <- set_epi(dat, "n_tx_doxi_f", length(ids_tx_doxi_f))

  # Update Attributes
  dat <- set_attr(dat, "curr_tx", curr_tx)
  dat <- set_attr(dat, "tx_day_eval", tx_day_eval)

  # Return dat object
  dat
}
