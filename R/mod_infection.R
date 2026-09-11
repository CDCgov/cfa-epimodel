#' @title Model Infection between Partners
#'
#' @description Disease transmission module accounting for
#' age-specific activity rates and increased transmission prob
#' for symptomatic persons
#'
#' @inheritParams vitals
#' @importFrom EpiModel get_attr set_attr set_epi discord_edgelist set_transmat
#'
#' @export

mod_infection_mgen <- function(dat, at) {
  # Notes
  ## NEED TESTS FOR SYMPTOMATIC MODIFIERS
  ## needs AMR tracker implementation (to add when AMR module is done)

  # Variables ---------------------------------------------------------------

  active <- get_attr(dat, "active")
  inf_time <- get_attr(dat, "inf_time")
  status <- get_attr(dat, "status")
  sympt <- get_attr(dat, "sympt")
  female <- get_attr(dat, "female")
  age_group <- get_attr(dat, "age_group")

  inf_prob_mtf <- get_param(dat, "inf_prob_mtf")
  inf_prob_ftm <- get_param(dat, "inf_prob_ftm")
  sympt_inf_modifier <- get_param(dat, "sympt_inf_modifier")
  act_rate_vec <- get_param(dat, "act_rate_vec")
  cond_prob_vec <- get_param(dat, "cond_prob_vec")
  cond_eff <- get_param(dat, "cond_eff")
  sympt_prob_m <- get_param(dat, "sympt_prob_m")
  sympt_prob_f <- get_param(dat, "sympt_prob_f")

  # Parameter Checks -------------------------------------------------------
  inf_probs <- c(
    inf_prob_mtf,
    inf_prob_ftm,
    cond_prob_vec,
    cond_eff,
    sympt_prob_m,
    sympt_prob_f
  )

  if (any(inf_probs < 0) || any(inf_probs > 1)) {
    stop(
      "All infection-related probabilities must be >=0 and <=1 (or null).",
      call. = FALSE
    )
  }

  if (sympt_inf_modifier < 1) {
    stop("sympt_inf_modifier parameter must be >= 1.", call. = FALSE)
  }

  if (
    (inf_prob_mtf * sympt_inf_modifier > 1 ||
      inf_prob_ftm * sympt_inf_modifier > 1)
  ) {
    stop(
      "Infection probabilities during symptomatic infection will exceed 1.
    Adjust inf_prob_mtf, inf_prob_ftm, or sympt_inf_modifier parameters.",
      call. = FALSE
    )
  }

  # Check that act_rate_vec length is valid
  n_age_groups <- length(unique(age_group[active == 1]))
  if (!(length(act_rate_vec) == 1 || length(act_rate_vec) == n_age_groups)) {
    stop(
      "act_rate_vec parameter length must be either 1 or
      equal to the number of age groups in the population.",
      call. = FALSE
    )
  }
  # Check that cond_prob_vec length is valid
  if (!(length(cond_prob_vec) == 1 || length(cond_prob_vec) == dat$num.nw)) {
    stop(
      "cond_prob_vec parameter length must be either 1 or
      equal to the number of networks in the simulation.",
      call. = FALSE
    )
  }
  # If a single condom-use probability is provided, apply it to all networks
  if (length(cond_prob_vec) == 1 && dat$num.nw > 1) {
    cond_prob_vec <- rep(cond_prob_vec, dat$num.nw)
  }

  # Process -----------------------------------------------------------------
  # Vector of infected and susceptible IDs
  ids_inf <- which(active == 1 & status == "i")
  n_active <- sum(active == 1)
  n_elig <- length(ids_inf)

  # Initialize vectors
  # G2 = female postscript (female attr == 1)
  n_inf <- n_inf_g2 <- tot_inf <- 0

  # If some infected AND some susceptible, then proceed
  if (n_elig > 0 && n_elig < n_active) {
    # Get discordant edgelist
    del_list <- lapply(
      seq_len(dat$num.nw),
      discord_edgelist,
      dat = dat,
      at = at,
      include.network = TRUE
    )
    del <- dplyr::bind_rows(del_list)

    # If some discordant edges, then proceed
    if (NROW(del) > 0) {
      #browser()
      # Infection duration to at
      del$inf_dur <- at - inf_time[del$inf]
      del$inf_dur[del$inf_dur == 0] <- 1

      # Calculate transmission rates based on sex & symptom status
      del$sympt_modifier <-
        ifelse(
          sympt[del$inf] == 1,
          sympt_inf_modifier,
          1
        )
      del$transProb <-
        ifelse(
          female[del$sus] == 1,
          inf_prob_mtf,
          inf_prob_ftm
        ) *
        del$sympt_modifier

      # Add age-group specific act rates
      # Using age of youngest partner in partnerships
      l_act_rate <- length(act_rate_vec)
      if (l_act_rate == 1) {
        del$act_rate <- act_rate_vec
      } else {
        age_inf <- age_group[del$inf]
        age_sus <- age_group[del$sus]
        age_youngest <- ifelse(age_inf < age_sus, age_inf, age_sus)
        del$act_rate <- act_rate_vec[age_youngest]
      }

      # STILL TO BE IMPLEMENTED
      # AMR tracker
      # AMR initial condition - add to attrs

      # Add probability of effective condom use based on rel type (network)
      cond_use_vec <- cond_prob_vec[del$network]
      ## Do they use condoms this time period?
      del$condUse <- stats::rbinom(nrow(del), 1, cond_use_vec)
      ## If using condoms
      ## apply effectiveness to reduce transmission probability
      ## condFinal = 0, no reduction in transmission probability
      ## condFinal = 1, 100% reduction
      del$condFinal <- del$condUse * cond_eff

      # Calculate final transmission probability per timestep
      del$finalProb <- 1 -
        (1 - (del$transProb * (1 - del$condFinal)))^del$act_rate

      # Randomize transmissions and subset df
      transmit <- stats::rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]

      # Set new infections attribute vectors
      ## status, infection time
      ids_new_inf <- unique(del$sus)

      if (length(ids_new_inf) > 0) {
        dat <- set_infection_attrs_mgen(dat, at, ids_new_inf)

        # Count new infections
        female_attrs <- unique(female)
        if (!all(female_attrs %in% c(0, 1))) {
          stop(
            "Female attribute must be coded as 0/1 only,
            which is required for sex-stratified epi slots.",
            call. = FALSE
          )
        }
        ## Calculate new infections among each sex explicitly by code:
        n_inf <- sum(female[ids_new_inf] == 0, na.rm = TRUE)
        n_inf_g2 <- sum(female[ids_new_inf] == 1, na.rm = TRUE)
        tot_inf <- n_inf + n_inf_g2
      }
    } # end some discordant edges condition
  } # end some active discordant nodes condition

  # Output ------------------------------------------------------------------

  # Save transmission matrix
  if (tot_inf > 0) {
    dat <- set_transmat(dat, del, at)
  }

  ## Save incidence vector
  dat <- set_epi(dat, "se_flow", at, tot_inf)
  dat <- set_epi(dat, "se_flow_m", at, n_inf)
  dat <- set_epi(dat, "se_flow_f", at, n_inf_g2)

  # Return
  dat
}

#' @title Assigns Infection-Related Attributes for Newly Exposed Individuals (M. genitalium)
#'
#' @description Assigns symptomatic status, calculates the time of transition from exposed
#' to infected based on the incubation period, calculates time of natural recovery.
#' This module should be run after the infection module in the same time step.
#' Note: This module assumes a normal distribution for the incubation period and infection duration,
#' with the mean specified in the parameters and a standard deviation equal to half of the mean.
#'
#' @inheritParams vitals
#' @param ids_new_e A vector of IDs corresponding to newly exposed individuals.
#'
#' @importFrom stats rnorm
#' @importFrom EpiModel get_attr set_attr get_param
#' @export
set_infection_attrs_mgen <- function(dat, at, ids_new_e) {
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

  n_new_e <- length(ids_new_e)

  # Assign infection status and time to newly exposed nodes
  status[ids_new_e] <- "e"
  inf_time[ids_new_e] <- at

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "inf_time", inf_time)

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
  ei_time[ids_new_e] <- at + incubation_period
  rec_time[ids_new_e] <- at + inf_dur_vec

  dat <- set_attr(dat, "ei_time", ei_time)
  dat <- set_attr(dat, "sympt", sympt)
  dat <- set_attr(dat, "rec_time", rec_time)

  # Return dat
  dat
}
