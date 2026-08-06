#' @title Modules for vital dynamics
#'
#' @description Handles node aging, departure, and arrivals
#'
#' @param dat Main \code{netsim_dat} object containing a \code{networkDynamic}
#'        object and other initialization information passed from
#'        \code{\link{EpiModel::netsim}}.
#' @param at Current time step.
#'
#' @name vitals
#' @rdname vitals
#' @export
mod_aging_mgen <- function(dat, at) {
  # Calc Updated Age Attributes
  age <- get_attr(dat, "age")
  age_group <- get_attr(dat, "age_group")
  units <- get_param(dat, "units_per_year")

  # Update age and age_group vectors
  age <- age + (1 / units)
  splits <- get_param(dat, "age_group_splits")
  splits_epi <- get_param(dat, "age_group_splits_epi")

  age_group <- cut(
    age,
    breaks = c(-Inf, splits, Inf),
    labels = FALSE,
    right = FALSE
  )
  age_group_epi <- cut(
    age,
    breaks = c(-Inf, splits_epi, Inf),
    labels = FALSE,
    right = FALSE
  )

  # Update Attributes
  dat <- set_attr(dat, "age", age)
  dat <- set_attr(dat, "age_group", age_group)
  dat <- set_attr(dat, "age_group_epi", age_group_epi)

  ## Summary statistics ##
  dat <- set_epi(dat, "meanAge", at, mean(age, na.rm = TRUE))

  # Return
  dat
}

# Departures Module ----------------------------------------------------
#' @rdname vitals
#' @export
mod_departures_mgen <- function(dat, at) {
  ## Attributes
  active <- get_attr(dat, "active")
  exitTime <- get_attr(dat, "exitTime")
  age <- get_attr(dat, "age")

  ## Parameters
  exit_age <- get_param(dat, "exit_age")

  ## if we had ASMR we would add that here

  ## Query alive but past simulation age range
  ## this setup a little odd make it easier to include ASMR later
  idsElig <- which(active == 1 & age >= exit_age)
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
mod_arrivals_mgen <- function(dat, at) {
  ## Parameters
  n <- sum(get_attr(dat, "active") == 1)
  aType <- get_param(dat, "arrivalType")

  ## Demographic attributes for new arrivals
  female_values <- get_param(dat, "entry_female_values")
  female_probs <- get_param(dat, "entry_female_probs")
  race_values <- get_param(dat, "entry_race_values")
  race_probs <- get_param(dat, "entry_race_probs")
  entry_age <- get_param(dat, "entry_age")
  entry_age_group <- get_param(dat, "age_group_splits")[[1]]
  entry_age_group_epi <- get_param(dat, "age_group_splits_epi")[[1]]

  ## Set up for new arrivals
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
      ## for small nArrivals, sample individually
      ## 5 is arbitrary cutoff but seems to work well in testing
      arrival_sex <- sample(
        female_values,
        nArrivals,
        prob = female_probs,
        replace = TRUE
      )
      arrival_race <- sample(
        race_values,
        nArrivals,
        prob = race_probs,
        replace = TRUE
      )
    } else {
      ## use base EpiModel apportion_lr function if nArrivals > 5
      arrival_sex <- apportion_lr(
        nArrivals,
        female_values,
        female_probs
      )
      arrival_race <- apportion_lr(nArrivals, race_values, race_probs)
    }

    ## Record length of attr vectors before new arrivals
    l <- length(get_attr_list(dat)[[1]])

    ## Update attributes for new arrivals
    ## EpiModel default core attrs: active, entryTime, exitTime, unique_id
    dat <- append_core_attr(dat, at, nArrivals)
    ## Custom attrs
    dat <- append_attr(dat, "status", "s", nArrivals)

    ### Required attrs for network formation: age, age_group, female, race
    dat <- append_attr(dat, "age", entry_age, nArrivals)
    dat <- append_attr(dat, "age_group", entry_age_group, nArrivals)
    dat <- append_attr(dat, "age_group_epi", entry_age_group_epi, nArrivals)
    dat <- append_attr(dat, "race", arrival_race, nArrivals)
    dat <- append_attr(dat, "female", arrival_sex, nArrivals)

    # Assign all other attrs NA (e.g. inf_time, rec_time, sympt, etc)
    # Attrs that need assignment have length equal to l,
    # the length of attr vectors before new arrivals
    attr_list <- get_attr_list(dat)
    attr_names <- names(attr_list)
    for (attr_name in attr_names) {
      if (length(attr_list[[attr_name]]) == l) {
        dat <- append_attr(dat, attr_name, NA, nArrivals)
      }
    }

    # Check that all attr vectors are now same length by pulling attr_list
    # again and checking lengths
    attr_lengths <- lengths(get_attr_list(dat))
    if (length(unique(attr_lengths)) != 1) {
      stop(paste0(
        "Not all attr vectors are same length after new arrivals at time ",
        at
      ))
    }
  }

  ## Summary statistics
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  ## Return
  dat
}
