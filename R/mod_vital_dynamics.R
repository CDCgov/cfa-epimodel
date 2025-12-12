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
  entryAge <- get_param(dat, "entryAge")
  exitAge <- get_param(dat, "exitAge")
  units <- get_param(dat, "units_per_year")

  # Update age and age_group vectors
  age <- age + (1 / units)
  ngrps <- ceiling((exitAge - entryAge) / age_group_width)
  # age groups: 0-19,20-24,25-29,30-34,35-39,40-44,45+, hardcoded
  for (i in seq_len(ngrps)) {
    nodes_in_group <- which(
      age >= (entryAge + (age_group_width * (i - 1))) &
        age < (entryAge + (age_group_width * i))
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
  exitAge <- get_param(dat, "exitAge")
  ## if we had ASMR we would add that here

  ## Query alive but past simulation age range
  ## this setup a little odd make it easier to include ASMR later
  idsElig <- which(active == 1 & ceiling(age) >= exitAge)
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
  entryAge <- get_param(dat, "entryAge")
  femaleProb <- get_param(dat, "entryFemaleProb")
  raceNames <- get_param(dat, "entryRaceNames")
  raceProbs <- get_param(dat, "entryRaceProbs")

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
    if (nArrivals == 1) {
      arrivalSex <- sample(c(0, 1), 1, prob = c(1 - femaleProb, femaleProb))
      arrivalRace <- sample(raceNames, 1, prob = raceProbs)
    } else { # use base EpiModel apportion_lr function if nArrivals > 1
      arrivalSex <- apportion_lr(nArrivals, c(0, 1), c(1 - femaleProb, femaleProb))
      arrivalRace <- apportion_lr(nArrivals, raceNames, raceProbs)
    }


    ## Update attributes for new arrivals
    dat <- append_core_attr(dat, at, nArrivals)
    dat <- append_attr(dat, "status", "s", nArrivals)
    dat <- append_attr(dat, "infTime", NA, nArrivals)
    dat <- append_attr(dat, "age", entryAge, nArrivals)
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
