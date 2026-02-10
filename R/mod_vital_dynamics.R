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
  age_group_width <- get_param(dat, "age_group_width")
  entry_age <- get_param(dat, "entry_age")
  exit_age <- get_param(dat, "exit_age")
  units <- get_param(dat, "units_per_year")

  # Update age and age_group vectors
  age <- age + (1 / units)
  ngrps <- ceiling((exit_age - entry_age) / age_group_width)
  for (i in seq_len(ngrps)) {
    nodes_in_group <- which(
      age >= (entry_age + (age_group_width * (i - 1))) &
        age < (entry_age + (age_group_width * i))
    )
    age_group[nodes_in_group] <- i
  }

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
  status <- get_attr(dat, "status")

  ## Parameters
  exit_age <- get_param(dat, "exit_age")

  ## Double check that entry age and attr.rules age are consistent
  ## Do this here because we currently use default epimodel arrival module
  entry_age <- get_param(dat, "entry_age")
  attr_rules_age <- get_control(dat, "attr.rules")$age
  if (!is.null(attr_rules_age) && entry_age != attr_rules_age) {
    stop("entry_age parameter and attr.rules age value must be the same.")
  }

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

    ## If departing nodes are in main relationship,
    ## Set partner's olderpartner attr to 1
    el <- get_edgelist(dat = dat, network = 1) # get edgelist
    rels1 <- which(el[, 1] %in% idsDept) # any departing id in row 1
    rels2 <- which(el[, 2] %in% idsDept) # any departing id in row 2
    allrels <- c(rels1, rels2) # combine

    if (length(allrels) > 0) {
      allParts <- el[allrels, ] # get all IDs (departing and partner)
      idsParts <- setdiff(allParts, idsDept) # extract partner IDs
      # Update Partner Attribute
      olderpartner <- get_attr(dat, "olderpartner")
      olderpartner[idsParts] <- 1
      dat <- set_attr(dat, "olderpartner", olderpartner)
    }
  }

  ## Reset attr
  dat <- set_attr(dat, "active", active)
  dat <- set_attr(dat, "exitTime", exitTime)

  ## Summary statistics
  dat <- set_epi(dat, "d.flow", at, nDepts)
  ### track number of edges in each network
  for (i in seq_along(dat$run$el)) {
    dat <- set_epi(dat, paste0("edges_net", i), at, nrow(dat$run$el[[i]]))
  }

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
    if (nArrivals <= 5) { # for small nArrivals, sample individually
      # 5 is arbitrary cutoff but seems to work well in testing
      arrivalSex <- sample(c(0, 1), nArrivals, prob = c(1 - female_prob, female_prob), replace = TRUE)
      arrivalRace <- sample(race_names, nArrivals, prob = race_probs, replace = TRUE)
    } else { # use base EpiModel apportion_lr function if nArrivals > 5
      arrivalSex <- apportion_lr(nArrivals, c(0, 1), c(1 - female_prob, female_prob))
      arrivalRace <- apportion_lr(nArrivals, race_names, race_probs)
    }


    ## Update attributes for new arrivals
    dat <- append_core_attr(dat, at, nArrivals)
    dat <- append_attr(dat, "status", "s", nArrivals)
    dat <- append_attr(dat, "infTime", NA, nArrivals)
    dat <- append_attr(dat, "age", entry_age, nArrivals)
    dat <- append_attr(dat, "age_group", 1, nArrivals)
    dat <- append_attr(dat, "race", arrivalRace, nArrivals)
    dat <- append_attr(dat, "female", arrivalSex, nArrivals)
    dat <- append_attr(dat, "olderpartner", 0, nArrivals)
  }

  ## Summary statistics
  dat <- set_epi(dat, "a.flow", at, nArrivals)

  # Return
  dat
}
