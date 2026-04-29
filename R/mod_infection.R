#' @title Model Infection between Partners
#'
#' @description Disease transmission module accounting for
#' age-specific activity rates and increased transmission prob
#' for symptomatic persons
#'
#' @inheritParams vitals
#' @importFrom EpiModel get_attr set_attr set_epi get_param
#' set_transmat discord_edgelist
#'
#' @export

mod_infection <- function(dat, at) {
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
      status[ids_new_inf] <- "i"
      dat <- set_attr(dat, "status", status)
      inf_time[ids_new_inf] <- at
      dat <- set_attr(dat, "inf_time", inf_time)
      ## symptomatic status
      sympt_prob_vec <- ifelse(
        female[ids_new_inf] == 1,
        sympt_prob_f,
        sympt_prob_m
      )
      sympt_vec <- which(rbinom(length(ids_new_inf), 1, sympt_prob_vec) == 1)
      ids_sympt <- ids_new_inf[sympt_vec]
      ids_asympt <- setdiff(ids_new_inf, ids_sympt)
      sympt[ids_sympt] <- 1
      sympt[ids_asympt] <- 0
      dat <- set_attr(dat, "sympt", sympt)

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
      ## si.flow.female0 corresponds to new infs among female == 0
      ## si.flow.female1 corresponds to new infs among female == 1
      n_inf <- sum(female[ids_new_inf] == 0, na.rm = TRUE)
      n_inf_g2 <- sum(female[ids_new_inf] == 1, na.rm = TRUE)
      tot_inf <- n_inf + n_inf_g2
    } # end some discordant edges condition
  } # end some active discordant nodes condition

  # Output ------------------------------------------------------------------

  # Save transmission matrix
  if (tot_inf > 0) {
    dat <- set_transmat(dat, del, at)
  }

  ## Save incidence vector
  dat <- set_epi(dat, "si.flow", at, tot_inf)
  dat <- set_epi(dat, "si.flow.female0", at, n_inf)
  dat <- set_epi(dat, "si.flow.female1", at, n_inf_g2)

  # Return
  dat
}
