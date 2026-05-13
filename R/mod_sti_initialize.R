#' @title Initialization Module
#'
#' @description Initialize networks, starting nodal attributes,
#' infection status,and other necessary components for STI
#' transmission models. Most of this module is adapted from
#' \code{EpiModel}'s \code{initialize.net} function, with bespoke
#' checks for required attributes and infection assignment (in future)
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
#' @importFrom EpiModel get_attr set_attr get_epi set_epi
#' get_param append_core_attr append_attr
#' get_edgelist init_status.net create_dat_object
#' init_nets sim_nets_t1 summary_nets
#' padded_vector get_control get_attr_prop set_param
#' @export

mod_sti_initialize <- function(x, param, init, control, s) {
  if (control$start == 1) {
    dat <- create_dat_object(param, init, control)
    dat <- init_nets(dat, x)
    if (!is.null(dat$run$nwterms)) {
      dat$run$t1.tab <- get_attr_prop(dat, dat$run$nwterms)
    }
    dat <- sim_nets_t1(dat)
    dat <- summary_nets(dat, at = 1L)

    # THIS IS THE ONLY CUSTOM CODE IN THIS FUNCTION, ALL ELSE IS FROM EpiModel::initialize.net
    dat <- init_mgen_status(dat) # initialize infection status and related attributes for M. genitalium
    # END OF CUSTOM CODE

    dat <- do.call(control[["prevalence.FUN"]], list(dat, at = 1))
  } else if (control$start > 1) {
    required_names <- c("param", "nwparam", "epi", "run", "coef.form", "num.nw")
    missing_names <- setdiff(required_names, names(x))
    if (length(missing_names) > 0) {
      stop(
        "x is missing the following elements required for re-initialization: ",
        paste.and(missing_names),
        call. = FALSE
      )
    }
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

  # Return dat
  dat
}

#' @title Initialize Infection Status
#' @description Initialize infection status and related attributes for STI
#' transmission models. This function is called within \code{mod_sti_initialize}
#' to set up the initial infection status of the population based on the number
#' of initial infections specified in the \code{init} object.
#' @inheritParams mod_sti_initialize
#' @return A modified \code{dat} object with initialized infections
#' @rdname mod_sti_initialize
#' @importFrom EpiModel get_init
#' @export
init_mgen_status <- function(dat) {
  num <- sum(get_attr(dat, "active") == 1)
  female <- get_attr(dat, "female")

  i_num <- get_init(dat, "i.num", override.null.error = TRUE)

  # Defaults
  status <- rep("s", num)
  inf_time <- sympt <- rec_time <- rep(NA, num)

  # Get initial infs
  ids_inf <- sample(num, i_num)
  status[ids_inf] <- "i"
  inf_time[ids_inf] <- 1

  ## symptomatic status
  sympt_prob_f <- get_param(dat, "sympt_prob_f")
  sympt_prob_m <- get_param(dat, "sympt_prob_m")

  sympt_prob_vec <- ifelse(
    female[ids_inf] == 1,
    sympt_prob_f,
    sympt_prob_m
  )
  sympt_vec <- which(rbinom(length(ids_inf), 1, sympt_prob_vec) == 1)
  ids_sympt <- ids_inf[sympt_vec]
  ids_asympt <- setdiff(ids_inf, ids_sympt)
  sympt[ids_sympt] <- 1
  sympt[ids_asympt] <- 0

  # Set attrs

  dat <- set_attr(dat, "status", status)
  dat <- set_attr(dat, "inf_time", inf_time)
  dat <- set_attr(dat, "sympt", sympt)
  dat <- set_attr(dat, "rec_time", rec_time)
  dat <- set_attr(dat, "curr_tx", rep(NA, num)) # initialize treatment status attribute
  dat <- set_attr(dat, "tx_end_day", rep(NA, num)) # initialize treatment evaluation day
  # THESE WILL NEED TO CHANGE ONCE AMR STATUS INCLUDED IN INITIALIZATION
  dat <- set_attr(dat, "amr_q", rep(0, num)) # initialize quinolone resistance status attribute (0 = susceptible, 1 = resistant)
  dat <- set_attr(dat, "amr_m", rep(0, num)) # initialize macrolide resistance status attribute (0 = susceptible, 1 = resistant)

  # Optional, save dat object for testing
  saveout <- get_control(dat, "save_dat")
  if (saveout) {
    folder_loc <- get_control(dat, 'save_dat_folder')
    saveRDS(
      dat,
      file.path(
        folder_loc,
        paste0("dat.rds")
      )
    )
  }

  # Return dat
  dat
}
