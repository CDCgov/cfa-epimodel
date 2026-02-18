#' @title Initialization Module
#'
#' @description Initialize networks, starting nodal attributes, infection status,
#' and other necessary components for STI transmission models. Most of this module
#' is adapted from \code{EpiModel}'s \code{initialize.net} function, with bespoke
#' checks for required attributes and infection assignment (in future)
#'
#'
#' @param x If \code{control$start == 1}, either a fitted network model object
#'        of class \code{netest} or a list of such objects. If
#'        \code{control$start > 1}, an object of class \code{netsim}. When
#'        multiple networks are used, the node sets (including network size
#'        and nodal attributes) are assumed to be the same for all networks.
#' @param param An \code{EpiModel} object of class \code{\link{param.net}}.
#' @param init An \code{EpiModel} object of class \code{\link{init.net}}.
#' @param control An \code{EpiModel} object of class \code{\link{control.net}}.
#' @param s Simulation number, used for restarting dependent simulations.
#' @details When re-initializing a simulation, the \code{netsim} object passed
#'          to \code{initialize.net} must contain the elements \code{param},
#'          \code{nwparam}, \code{epi}, \code{coef.form}, and \code{num.nw}.
#'
#' @return A \code{netsim_dat} class main data object.
#'
#' @importFrom EpiModel get_attr set_attr get_epi set_epi get_param append_core_attr append_attr
#' get_edgelist init_status.net create_dat_object init_nets sim_nets_t1 summary_nets
#' padded_vector get_control get_attr_prop set_param
#' @export

mod_sti_initialize <- function(x, param, init, control, s) {
  if (control$start == 1) {
    # Main Data List --------------------------------------------------------
    dat <- create_dat_object(param, init, control)

    # network and stats initialization
    dat <- init_nets(dat, x)

    ## Store current proportions of attr
    if (!is.null(dat$run$nwterms)) {
      dat$run$t1.tab <- get_attr_prop(dat, dat$run$nwterms)
    }

    # simulate first time step
    dat <- sim_nets_t1(dat)
    dat <- summary_nets(dat, at = 1L)

    ## Check for attrs required on initialization for use in vital dynamics modules
    required_attrs_init <- c("age", "age_group", "olderpartner", "female", "race")
    current_attrs <- dat$run$nwterms
    missing_init_attrs <- setdiff(required_attrs_init, current_attrs)
    if (length(missing_init_attrs) > 0) {
      stop(
        "The following required attributes are missing from the network at initialization: ",
        paste.and(missing_init_attrs),
        call. = FALSE
      )
    }

    ## Infection Status and Time
    dat <- init_status.net(dat)

    # Summary Stats -----------------------------------------------------------
    dat <- do.call(control[["prevalence.FUN"]], list(dat, at = 1))

    # Restart/Reinit Simulations ----------------------------------------------
  } else if (control$start > 1) {
    ## check that required names are present
    required_names <- c("param", "nwparam", "epi", "run", "coef.form", "num.nw")
    missing_names <- setdiff(required_names, names(x))
    if (length(missing_names) > 0) {
      stop(
        "x is missing the following elements required for re-initialization: ",
        paste.and(missing_names),
        call. = FALSE
      )
    }

    # recycle sims in the restart object
    # e.g. 5 sim out of a size 3 restart object we will give: 1, 2, 3, 1, 2
    s <- (s - 1) %% length(x$run) + 1

    dat <- create_dat_object(
      param = param,
      control = control,
      run = x$run[[s]]
    )

    missing_params <- setdiff(names(x$param), names(param))
    for (mp in missing_params) {
      dat <- set_param(dat, mp, x$param[[mp]])
    }

    dat$num.nw <- x$num.nw

    dat$nwparam <- x$nwparam
    for (network in seq_len(dat$num.nw)) {
      dat$nwparam[[network]]$coef.form <- x$coef.form[[s]][[network]]
    }
    dat$epi <- sapply(x$epi, function(var) var[s])
    names(dat$epi) <- names(x$epi)

    dat$stats <- lapply(x$stats, function(var) var[[s]])
    if (get_control(dat, "save.nwstats")) {
      nsteps <- get_control(dat, "nsteps")
      start <- get_control(dat, "start")
      dat$stats$nwstats <- lapply(
        dat$stats$nwstats,
        function(oldstats) padded_vector(list(oldstats), nsteps - start + 2L)
      )
    }
    if (is.data.frame(dat$stats$transmat)) {
      nsteps <- get_control(dat, "nsteps")
      dat$stats$transmat <- padded_vector(list(dat$stats$transmat), nsteps)
    }
  }

  # Return
  dat
}
