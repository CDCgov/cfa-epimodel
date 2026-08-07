#' @title Transition from Exposed to Infected
#'
#' @description Evaluates eligible exposed persons for infection
#' based on sex-specific rates, updates infection status accordingly
#'
#' @inheritParams vitals
#'
#' @export
mod_ei_mgen <- function(dat, at) {
  # Get attributes
  active <- get_attr(dat, "active")
  female <- get_attr(dat, "female")
  status <- get_attr(dat, "status")
  ei_time <- get_attr(dat, "ei_time")

  # Initialize epi trackers
  tot_ei <- n_ei_m <- n_ei_f <- 0

  # Identify nodes transitioning from E to I
  ids_new_i <- which(active == 1 & status == "e" & ei_time == at)
  l_ei <- length(ids_new_i)

  # Update status for nodes transitioning from E to I
  if (l_ei > 0) {
    tot_ei <- l_ei
    n_ei_m <- sum(female[ids_new_i] == 0, na.rm = TRUE)
    n_ei_f <- sum(female[ids_new_i] == 1, na.rm = TRUE)
    status[ids_new_i] <- "i"
    dat <- set_attr(dat, "status", status)
  }

  # Epi Trackers
  dat <- set_epi(dat, "ei_flow", at, tot_ei)
  dat <- set_epi(dat, "ei_flow_m", at, n_ei_m)
  dat <- set_epi(dat, "ei_flow_f", at, n_ei_f)
}
