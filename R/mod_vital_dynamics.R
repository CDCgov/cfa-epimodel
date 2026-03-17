#' @title Modules for vital dynamics
#'
#' @description Handles node aging, departure, and arrivals
#'
#' @param dat Main \code{netsim_dat} object containing a \code{networkDynamic}
#'        object and other initialization information passed from
#'        \code{\link{netsim}}.
#' @param at Current time step.
#'
#' @importFrom EpiModel get_attr set_attr get_epi set_epi get_param append_core_attr append_attr
#' get_edgelist apportion_lr
#'
#' @name vitals
NULL

#' @rdname vitals
#' @export
mod_aging <- function(dat, at) {
  # Calc Updated Age Attributes
  age <- get_attr(dat, "age")
  age_group <- get_attr(dat, "age_group")
  units <- get_param(dat, "units_per_year")

  # Update age and age_group vectors
  age <- age + (1 / units)
  age_group <- dplyr::case_when(
    age < 50 & age >= 45 ~ 7,
    age < 45 & age >= 40 ~ 6,
    age < 40 & age >= 35 ~ 5,
    age < 35 & age >= 30 ~ 4,
    age < 30 & age >= 25 ~ 3,
    age < 25 & age >= 19 ~ 2,
    age < 19 ~ 1
  )

  # Update Attributes
  dat <- set_attr(dat, "age", age)
  dat <- set_attr(dat, "age_group", age_group)

  ## Summary statistics ##
  dat <- set_epi(dat, "meanAge", at, mean(age, na.rm = TRUE))

  # Return
  dat
}

# Departures Module ----------------------------------------------------
#' @rdname vitals
#' @export
mod_departures <- function(dat, at) {
  ## Attributes
  active <- get_attr(dat, "active")
  exitTime <- get_attr(dat, "exitTime")
  age <- get_attr(dat, "age")

  ## Parameters
  exit_age <- get_param(dat, "exit_age")

  ## if we had ASMR we would add that here

  ## Query alive but past simulation age range
  ## this setup a little odd make it easier to include ASMR later
  idsElig <- which(active == 1 & ceiling(age) >= exit_age)
  nElig <- length(idsElig)
  nDepts <- 0

  if (nElig > 0) {
    idsDept <- idsElig
    nDepts <- length(idsDept)
    ## Update nodal attributes
    active[idsDept] <- 0
    exitTime[idsDept] <- at
  }

  ## Reset attr
  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  ## Summary statistics
  dat <- set_epi(dat, "d.flow", at, nDepts)

  # Return
  dat
}


# Arrivals Module ----------------------------------------------------
#' @rdname vitals
#' @export
mod_arrivals <- function(dat, at) {
  ## Parameters
  n <- sum(get_attr(dat, "active") == 1)
  aType <- get_param(dat, "arrivalType")
  female_prob <- get_param(dat, "entry_female_prob")
  race_probs <- get_param(dat, "entry_race_probs")
  race_names <- get_param(dat, "entry_race_names")
  entry_age <- get_param(dat, "entry_age")

  nArrivals <- 0

  if (!aType %in% c("rate", "departures")) {
    stop("Arrival Type must be either 'rate' or 'departures'")
  }

  if (aType == "rate") {
    a_rate <- get_param(dat, "arrival.rate")

    ## Process
    nArrivalsExp <- n * a_rate
    nArrivals <- stats::rpois(1, nArrivalsExp)
  }

  if (aType == "departures") {
    nArrivals <- get_epi(dat, "d.flow", at)
  }

  if (nArrivals > 0) {
    ## Determine sex, race
    if (nArrivals <= 5) {
      # for small nArrivals, sample individually
      # 5 is arbitrary cutoff but seems to work well in testing
      arrival_sex <- sample(
        c(0, 1),
        nArrivals,
        prob = c(1 - female_prob, female_prob),
        replace = TRUE
      )
      arrival_race <- sample(
        race_names,
        nArrivals,
        prob = race_probs,
        replace = TRUE
      )
    } else {
      # use base EpiModel apportion_lr function if nArrivals > 5
      arrival_sex <- apportion_lr(
        nArrivals,
        c(0, 1),
        c(1 - female_prob, female_prob)
      )
      arrival_race <- apportion_lr(nArrivals, race_names, race_probs)
    }

    ## Update attributes for new arrivals
    dat <- append_core_attr(dat, at, nArrivals)
    dat <- append_attr(dat, "status", "s", nArrivals)
    dat <- append_attr(dat, "infTime", NA, nArrivals)
    dat <- append_attr(dat, "age", entry_age, nArrivals)
    dat <- append_attr(dat, "age_group", 1, nArrivals)
    dat <- append_attr(dat, "race", arrival_race, nArrivals)
    dat <- append_attr(dat, "female", arrival_sex, nArrivals)
  }

  ## Summary statistics
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  # Return
  dat
}
