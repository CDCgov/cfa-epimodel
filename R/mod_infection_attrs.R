#' @title Assigns Infection-Related Attributes for Newly Exposed Individuals (M. genitalium)
#'
#' @description Assigns symptomatic status, calculates the time of transition from exposed
#' to infected based on the incubation period, calculates time of natural recovery.
#' This module should be run after the infection module in the same time step.
#' Note: This module assumes a normal distribution for the incubation period and infection duration,
#' with the mean specified in the parameters and a standard deviation equal to half of the mean.
#'
#' @inheritParams vitals
#'
#' @export

mod_infection_attrs_mgen <- function(dat, at) {
  # Ensure this module run after infection module
  se_flow <- get_epi(dat, "se_flow", at)
  if (is.na(se_flow)) {
    stop(
      "Symptoms module must be run after infection module in the same time step.",
      call. = FALSE
    )
  }

  # Proceed
  # Get attributes
  active <- get_attr(dat, "active")
  female <- get_attr(dat, "female")
  status <- get_attr(dat, "status")
  sympt <- get_attr(dat, "sympt")
  inf_time <- get_attr(dat, "inf_time")
  ei_time <- get_attr(dat, "ei_time")
  rec_time <- get_attr(dat, "rec_time")

  # Get parameters
  sympt_prob_m <- get_param(dat, "sympt_prob_m")
  sympt_prob_f <- get_param(dat, "sympt_prob_f")
  mean_incubation_period <- get_param(dat, "mean_incubation_period")
  mean_inf_dur_m <- get_param(dat, "mean_infection_duration_m")
  mean_inf_dur_f <- get_param(dat, "mean_infection_duration_f")

  # Assign E to I time for newly exposed nodes
  ## Get IDs of newly exposed nodes
  ids_new_e <- which(active == 1 & status == "e" & inf_time == at)
  n_new_e <- length(ids_new_e)

  if (n_new_e > 0) {
    ## Assign symptomatic status for newly exposed nodes based on sex-specific probabilities
    sympt_prob_vec <- ifelse(
      female[ids_new_e] == 1,
      sympt_prob_f,
      sympt_prob_m
    )
    sympt_ids <- get_successful_ids_binom(ids_new_e, sympt_prob_vec)
    asympt_ids <- setdiff(ids_new_e, sympt_ids)

    ## Calculate incubation period for each newly exposed node
    ## (not sex-specific, but could be modified to be if desired)
    incubation_period <- ceiling(rnorm(
      n_new_e,
      mean = mean_incubation_period,
      sd = mean_incubation_period / 2
    ))

    ## Calculate infection duration for newly infected nodes
    ## based on sex-specific means
    inf_dur_vec <- ifelse(
      female[ids_new_e] == 1,
      ceiling(rnorm(n_new_e, mean = mean_inf_dur_f, sd = mean_inf_dur_f / 2)),
      ceiling(rnorm(n_new_e, mean = mean_inf_dur_m, sd = mean_inf_dur_m / 2))
    )

    ## Update sympt, ei_time, and rec_time attributes
    sympt[sympt_ids] <- 1
    sympt[asympt_ids] <- 0
    dat <- set_attr(dat, "sympt", sympt)
    ei_time[ids_new_e] <- at + incubation_period
    dat <- set_attr(dat, "ei_time", ei_time)
    rec_time[ids_new_e] <- at + inf_dur_vec
    dat <- set_attr(dat, "rec_time", rec_time)
  }

  # Return dat
  dat
}
