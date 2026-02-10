#' @title Model Infection between Partners
#'
#' @description Disease transmission module accounting for age-specific activity rates
#'
#' @inheritParams vitals
#' @importFrom EpiModel get_attr set_attr set_epi get_param set_transmat discord_edgelist
#'
#' @export

mod_infection <- function(dat, at) {
  # Notes
  ## NEED TESTS FOR INFECTION STAGE MODIFIERS
  ## NEEDS TESTS FOR CONDOM USE PROBABILITIES & EFFECTIVENESS
  ## needs AMR tracker implementation (to add when AMR module is done)

  # Variables ---------------------------------------------------------------

  active <- get_attr(dat, "active")
  infTime <- get_attr(dat, "infTime")
  status <- get_attr(dat, "status")
  female <- get_attr(dat, "female")
  age_group <- get_attr(dat, "age_group")

  inf_prob_mtf <- get_param(dat, "inf_prob_mtf")
  inf_prob_ftm <- get_param(dat, "inf_prob_ftm")
  acute_inf_modifier <- get_param(dat, "acute_inf_modifier")
  acute_duration <- get_param(dat, "acute_duration")
  act_rate_vec <- get_param(dat, "act_rate_vec")
  cond_prob_vec <- get_param(dat, "cond_prob_vec")
  cond_eff <- get_param(dat, "cond_eff")

  # Parameter Checks -------------------------------------------------------
  inf_probs <- c(inf_prob_mtf, inf_prob_ftm, cond_prob_vec, cond_eff)
  if (any(inf_probs < 0) || any(inf_probs > 1)) {
    stop("All infection-related probabilities must be >=0 and <=1 (or null).", call. = FALSE)
  }

  if (acute_inf_modifier < 1) {
    stop("acute_inf_modifier parameter must be >= 1.", call. = FALSE)
  }

  # Generate inf prob vectors if acute modifier > 1 & acute_duration defined
  if (acute_inf_modifier > 1 && !is.null(acute_duration)) {
    inf_prob_mtf <- c(rep(inf_prob_mtf * acute_inf_modifier, acute_duration), inf_prob_mtf)
    inf_prob_ftm <- c(rep(inf_prob_ftm * acute_inf_modifier, acute_duration), inf_prob_ftm)
    if (any(inf_prob_mtf > 1) || any(inf_prob_ftm > 1)) {
      stop("Infection probabilities during acute infection stage exceed 1.
      Adjust inf_prob_mtf, inf_prob_ftm, or acute_inf_modifier parameters.", call. = FALSE)
    }
  }

  # Check that act_rate_vec length is valid
  n_age_groups <- length(unique(age_group[active == 1]))
  if (!(length(act_rate_vec) == 1 || length(act_rate_vec) == n_age_groups)) {
    stop("act_rate_vec parameter length must be either 1 or equal to the number of age groups in the population.",
      call. = FALSE
    )
  }
  # Check that cond_prob_vec length is valid
  if (!(length(cond_prob_vec) == 1 || length(cond_prob_vec) == dat$num.nw)) {
    stop("cond_prob_vec parameter length must be either 1 or equal to the number of networks in the simulation.",
      call. = FALSE
    )
  }

  # Process -----------------------------------------------------------------
  # Vector of infected and susceptible IDs
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)

  # Initialize vectors
  # G2 = female postscript (female attr == 1)
  nInf <- nInfG2 <- totInf <- 0

  # If some infected AND some susceptible, then proceed
  if (nElig > 0 && nElig < nActive) {
    # Get discordant edgelist
    del_list <- lapply(seq_len(dat$num.nw), discord_edgelist, dat = dat, at = at, include.network = TRUE)
    del <- dplyr::bind_rows(del_list)

    # If some discordant edges, then proceed
    if (NROW(del) > 0) {
      # Infection duration to at
      del$infDur <- at - infTime[del$inf]
      del$infDur[del$infDur == 0] <- 1

      # Calculate infection-stage transmission rates
      linf.prob <- length(inf_prob_mtf)
      if (is.null(inf_prob_ftm)) {
        del$transProb <- ifelse(del$infDur <= linf.prob,
          inf_prob_mtf[del$infDur],
          inf_prob_mtf[linf.prob]
        )
      } else {
        del$transProb <- ifelse(female[del$sus] == 1,
          ifelse(del$infDur <= linf.prob,
            inf_prob_mtf[del$infDur],
            inf_prob_mtf[linf.prob]
          ),
          ifelse(del$infDur <= linf.prob,
            inf_prob_ftm[del$infDur],
            inf_prob_ftm[linf.prob]
          )
        )
      }

      # Add age-group specific act rates
      # Using age of youngest partner in partnerships
      lact.rate <- length(act_rate_vec)
      if (length(act_rate_vec) == 1) {
        del$actRate <- act_rate_vec
      } else {
        age_inf <- age_group[del$inf]
        age_sus <- age_group[del$sus]
        age_youngest <- ifelse(age_inf < age_sus, age_inf, age_sus)
        del$actRate <- act_rate_vec[age_youngest]
      }

      # STILL TO BE IMPLEMENTED
      # AMR tracker
      # AMR initial condition - add to attrs

      # Add probability of effective condom use based on rel type (network)
      cond_use_vec <- cond_prob_vec[del$network]
      ## Do they use condoms this time period?
      del$condUse <- stats::rbinom(nrow(del), 1, cond_use_vec)
      ## If using condoms, apply effectiveness to reduce transmission probability
      ## condFinal = 0, no reduction in transmission probability
      ## condFinal = 1, 100% reduction
      del$condFinal <- del$condUse * cond_eff

      # Calculate final transmission probability per timestep
      del$finalProb <- 1 - (1 - (del$transProb * (1 - del$condFinal)))^del$actRate

      # Randomize transmissions and subset df
      transmit <- stats::rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]

      # Set new infections attribute vectors
      idsNewInf <- unique(del$sus)
      status[idsNewInf] <- "i"
      dat <- set_attr(dat, "status", status)
      infTime[idsNewInf] <- at
      dat <- set_attr(dat, "infTime", infTime)

      # Count new infections
      female_attrs <- unique(female)
      if (length(female_attrs) > 2) {
        stop("Female attribute has more than 2 unique values, which is not currently supported.", call. = FALSE)
      }
      ## If there are both females and males in the population,
      ## calculate new infections among each sex
      if (length(female_attrs) == 2) {
        nInf <- sum(female[idsNewInf] == female_attrs[1])
        nInfG2 <- sum(female[idsNewInf] == female_attrs[2])
      } else {
        nInf <- sum(female[idsNewInf] == female_attrs[1])
        nInfG2 <- 0
      }
      totInf <- nInf + nInfG2
    } # end some discordant edges condition
  } # end some active discordant nodes condition


  # Output ------------------------------------------------------------------

  # Save transmission matrix
  if (totInf > 0) {
    dat <- set_transmat(dat, del, at)
  }

  ## Save incidence vector
  dat <- set_epi(dat, "si.flow", at, totInf)
  dat <- set_epi(dat, "si.flow.female0", at, nInf)
  dat <- set_epi(dat, "si.flow.female1", at, nInfG2)

  # Return
  dat
}
