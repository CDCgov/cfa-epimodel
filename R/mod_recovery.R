#' @title Recover from M. genitalium infection
#'
#' @description Evaluates eligible infected persons for recovery
#' based on sex-specific rates, resets infection attrs if recovered
#'
#' @inheritParams vitals
#'
#' @export
mod_recovery_mgen <- function(dat, at) {
  # Ensure this module run after infection and clinical modules
  # (if clinical module is included in the model)
  se_flow <- get_epi(dat, "se_flow", at)

  clinical <- "clinical.FUN" %in% names(get_control_list(dat))

  if (clinical) {
    n_naat_tests <- get_epi(dat, "n_naat_tests", at)
    if (is.na(se_flow) || is.na(n_naat_tests)) {
      stop(
        "Recovery module must be run after both infection and clinical modules in the same time step.",
        call. = FALSE
      )
    }
  } else {
    if (is.na(se_flow)) {
      stop(
        "Recovery module must be run after the infection module in the same time step.",
        call. = FALSE
      )
    }
  }

  # Parameters
  rec_state <- get_param(dat, "rec_state") # recovery state (e.g., "s" for susceptible)

  # Attributes
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  female <- get_attr(dat, "female")
  ## Attrs that will need to be updated for recovered nodes
  sympt <- get_attr(dat, "sympt")
  inf_time <- get_attr(dat, "inf_time")
  ei_time <- get_attr(dat, "ei_time")
  rec_time <- get_attr(dat, "rec_time")
  amr_q <- get_attr(dat, "amr_q")
  amr_m <- get_attr(dat, "amr_m")
  curr_tx <- get_attr(dat, "curr_tx")
  tx_end_day <- get_attr(dat, "tx_end_day")
  tx_success <- get_attr(dat, "tx_success")

  # Setup initial vectors
  n_recov <- n_recov_m <- n_recov_f <- 0

  # Identify eligible infected nodes for recovery
  ids_recov <- which(active == 1 & status == "i" & rec_time == at)

  if (clinical) {
    # If clinical module is included, add nodes with tx success
    ids_recov_tx <- which(active == 1 & status == "i" & tx_success == 1)
    ids_recov <- unique(c(ids_recov, ids_recov_tx))
  }

  if (length(ids_recov) > 0) {
    n_recov <- length(ids_recov)
    n_recov_m <- sum(female[ids_recov] == 0, na.rm = TRUE)
    n_recov_f <- sum(female[ids_recov] == 1, na.rm = TRUE)

    # Update attributes for recovered nodes
    status[ids_recov] <- rec_state
    sympt[ids_recov] <- NA
    inf_time[ids_recov] <- NA
    ei_time[ids_recov] <- NA
    rec_time[ids_recov] <- NA
    amr_q[ids_recov] <- NA
    amr_m[ids_recov] <- NA
    curr_tx[ids_recov] <- NA
    tx_end_day[ids_recov] <- NA
    tx_success[ids_recov] <- NA
  }

  # Update attrs
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "sympt", sympt)
  dat <- set_attr(dat, "rec_time", rec_time)
  dat <- set_attr(dat, "inf_time", inf_time)
  dat <- set_attr(dat, "ei_time", ei_time)
  dat <- set_attr(dat, "amr_q", amr_q)
  dat <- set_attr(dat, "amr_m", amr_m)
  dat <- set_attr(dat, "curr_tx", curr_tx)
  dat <- set_attr(dat, "tx_end_day", tx_end_day)
  dat <- set_attr(dat, "tx_success", tx_success)

  # Update epi
  flow_name <- paste0("i", rec_state, "_flow")
  dat <- set_epi(dat, flow_name, at, n_recov)
  dat <- set_epi(dat, paste0(flow_name, "_m"), at, n_recov_m)
  dat <- set_epi(dat, paste0(flow_name, "_f"), at, n_recov_f)

  # Return dat
  dat
}
