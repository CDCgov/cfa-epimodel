#' @title Updates to Network Attributes During Mpox Resimulation
#'
#' @description This function updates the network-related degree
#' attributes on the three-layer MSM sexual network that account for the
#' resimulated network structure.
#'
#' @inheritParams vitals
#' @param network Integer for network number (values 1, 2, or 3 for main,
#' casual, and one-time networks).
#'
#' @details
#' This function is called between network resimulations in
#' [`EpiModel::resim_nets`], passed into [`EpiModel::control.net`]
#' through the `dat.updates` argument. This implementation updates degree
#' attributes calculated from the current network snapshot for use as ERGM
#' terms in the other network layers (e.g., degree in the casual network
#' is a function of the degree in the main network). See the
#' general documentation for `dat.updates` at [`EpiModel::control.net`].
#'
#'
#' @export
#'
resimnet_updates_sti <- function(dat, at, network) {
  if (network == 0L) {
    dat <- EpiModel::set_attr(
      dat,
      "deg_casual",
      EpiModel::get_degree(dat$run$el[[2]])
    )
  } else if (network == 1L) {
    dat <- EpiModel::set_attr(
      dat,
      "deg_main",
      EpiModel::get_degree(dat$run$el[[1]])
    )
  } else if (network == 2L) {
    dat <- EpiModel::set_attr(
      dat,
      "deg_casual",
      EpiModel::get_degree(dat$run$el[[2]])
    )
    dat <- EpiModel::set_attr(
      dat,
      "deg_tot",
      EpiModel::get_attr(dat, "deg_main") +
        EpiModel::get_degree(dat$run$el[[2]])
    )
  }
  dat
}

#' @title Helper Function for Binomial Sampling of Eligible IDs
#' @description This function takes a vector of eligible IDs and a corresponding
#' vector of probabilities, and returns a vector of IDs that are "successful"
#' based on a binomial draw. If there are no eligible IDs, it returns an empty integer vector.
#' @param elig_ids A vector of eligible IDs (e.g., indices of symptomatic infected individuals).
#' @param eval_prob A vector of probabilities corresponding to each eligible ID, or a single
#' probability to be applied to all eligible IDs.
#' @return A vector of IDs or an empty integer vector.
#' @export
get_successful_ids_binom <- function(elig_ids, eval_prob) {
  l_elig_ids <- length(elig_ids)
  l_eval_prob <- length(eval_prob)
  if (length(elig_ids) > 0) {
    if (l_eval_prob %in% c(1, l_elig_ids)) {
      elig_ids[rbinom(length(elig_ids), 1, eval_prob) == 1]
    } else {
      stop(
        "Length of eval_prob must be either 1 or equal to the number of eligible IDs."
      )
    }
  } else {
    integer(0) # return empty integer vector if no eligible IDs
  }
}
