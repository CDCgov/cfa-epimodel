#' @title Track Number of Relationships in Each Network
#'
#' @description This module tracks the number of edges in each network at each
#' time step, storing the counts in the epidemiological summary statistics list.
#'
#' @inheritParams vitals
#'
#' @export

track_edges <- function(dat, at) {
  # Loop through edgeslists for each network and set summary stat
  for (i in seq_along(dat$run$el)) {
    dat <- set_epi(dat, paste0("edges_net_", i), at, nrow(dat$run$el[[i]]))
  }

  # Return dat object
  dat
}
