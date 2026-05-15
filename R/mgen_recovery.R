#' @title Recover from M. genitalium infection
#'
#' @description Evaluates eligible infected persons for recovery
#' based on sex-specific rates, resets infection attrs if recovered
#'
#' @inheritParams vitals
#' @importFrom EpiModel get_attr set_attr set_epi get_param
#'
#' @export
mod_mgen_recovery <- function(dat, at) {
  # Parameters
  rec_rate_m <- 1 / get_param(dat, "inf_dur_m")
  rec_rate_f <- 1 / get_param(dat, "inf_dur_f")
  rec_state <- get_param(dat, "rec_state")

  # Attributes
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  inf_time <- get_attr(dat, "inf_time")
  rec_time <- get_attr(dat, "rec_time")
  female <- get_attr(dat, "female")
  sympt <- get_attr(dat, "sympt")

  # Setup initial vectors
  n_recov <- n_recov_m <- n_recov_f <- 0
  ids_elig <- which(active == 1 & status == "i")
  n_elig <- length(ids_elig)

  # For now, assume exponential dist for recovery
  # leaving inf_dur calc here in case that changes
  inf_dur <- at - inf_time[active == 1 & status == "i"]
  inf_dur[inf_dur == 0] <- 1

  if (n_elig > 0) {
    female_vec <- female[ids_elig]
    rate_vec <- ifelse(female_vec == 1, rec_rate_f, rec_rate_m)
    recov_vec <- which(rbinom(n_elig, 1, rate_vec) == 1)
    if (length(recov_vec) > 0) {
      ids_recov <- ids_elig[recov_vec]
      n_recov <- length(ids_recov)
      n_recov_f <- sum(female[ids_recov] == 1)
      n_recov_m <- sum(female[ids_recov] == 0)
      status[ids_recov] <- rec_state
      rec_time[ids_recov] <- at
      inf_time[ids_recov] <- NA
      sympt[ids_recov] <- NA
    }
  }

  # Update attrs
  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "sympt", sympt)
  dat <- set_attr(dat, "rec_time", rec_time)
  dat <- set_attr(dat, "inf_time", inf_time)

  # Update epi
  flow_name <- paste0("i", rec_state, "_flow")
  dat <- set_epi(dat, flow_name, at, n_recov)
  dat <- set_epi(dat, paste0(flow_name, "_m"), at, n_recov_m)
  dat <- set_epi(dat, paste0(flow_name, "_f"), at, n_recov_f)

  # Return dat
  dat
}
