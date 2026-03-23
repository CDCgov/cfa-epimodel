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
    dat <- init_mgen_status(dat)
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

#' @export
init_mgen_status <- function(dat) {
  type <- get_control(dat, "type", override.null.error = TRUE)
  type <- if (is.null(type)) {
    "None"
  } else {
    type
  }
  nsteps <- get_control(dat, "nsteps")
  tergmLite <- get_control(dat, "tergmLite")
  vital <- get_param(dat, "vital")
  groups <- get_param(dat, "groups")
  status.vector <- get_init(dat, "status.vector", override.null.error = TRUE)
  if (type %in% c("SIS", "SIR")) {
    rec.rate <- get_param(dat, "rec.rate")
  }
  if (vital) {
    di.rate <- get_param(dat, "di.rate")
  }
  i.num <- get_init(dat, "i.num", override.null.error = TRUE)
  if (type == "SIR" && is.null(status.vector)) {
    r.num <- get_init(dat, "r.num")
  }
  num <- sum(get_attr(dat, "active") == 1)
  if (groups == 2) {
    group <- get_attr(dat, "group")
    if (!all(group %in% c(1, 2))) {
      stop(
        "When using the `group` attribute, the only authorized values",
        " are 1 and 2.\n",
        "The values found were: ",
        paste0(unique(group), collapse = ", ")
      )
    }
    i.num.g2 <- get_init(dat, "i.num.g2")
    if (type == "SIR" && is.null(status.vector)) {
      r.num.g2 <- get_init(dat, "r.num.g2", override.null.error = TRUE)
    }
  } else {
    group <- rep(1, num)
  }
  statOnNw <- "status" %in% dat$run$nwterms
  if (!statOnNw) {
    if (!is.null(status.vector)) {
      status <- status.vector
    } else {
      status <- rep("s", num)
      status[sample(which(group == 1), size = i.num)] <- "i"
      if (groups == 2) {
        status[sample(which(group == 2), size = i.num.g2)] <- "i"
      }
      if (type == "SIR") {
        status[sample(which(group == 1 & status == "s"), size = r.num)] <- "r"
        if (groups == 2) {
          status[sample(
            which(group == 2 & status == "s"),
            size = r.num.g2
          )] <- "r"
        }
      }
    }
    dat <- set_attr(dat, "status", status)
  } else {
    status <- get_attr(dat, "status")
  }
  if (!tergmLite) {
    if (!statOnNw) {
      for (network in seq_len(dat$num.nw)) {
        dat$run$nw[[network]] <- set_vertex_attribute(
          dat$run$nw[[network]],
          "status",
          status
        )
      }
    }
    for (network in seq_len(dat$num.nw)) {
      dat$run$nw[[network]] <- activate.vertex.attribute(
        dat$run$nw[[network]],
        prefix = "testatus",
        value = status,
        onset = 1,
        terminus = Inf
      )
    }
  }
  if (type == "None") {
    infTime <- rep(NA, num)
    idsInf <- idsInf <- which(status == "i")
    infTime[idsInf] <- 1
    dat <- set_attr(dat, "infTime", infTime)
  } else {
    idsInf <- which(status == "i")
    infTime <- rep(NA, length(status))
    infTime.vector <- get_init(
      dat,
      "infTime.vector",
      override.null.error = TRUE
    )
    if (!is.null(infTime.vector)) {
      infTime <- infTime.vector
    } else {
      if (vital && di.rate > 0) {
        if (type == "SI") {
          infTime[idsInf] <- -rgeom(n = length(idsInf), prob = di.rate) + 2
        } else {
          infTime[idsInf] <- -rgeom(
            n = length(idsInf),
            prob = di.rate + (1 - di.rate) * mean(rec.rate)
          ) +
            2
        }
      } else {
        if (type == "SI" || mean(rec.rate) == 0) {
          infTime[idsInf] <- ssample(
            1:(-nsteps + 2),
            length(idsInf),
            replace = TRUE
          )
        } else {
          infTime[idsInf] <- -rgeom(n = length(idsInf), prob = mean(rec.rate)) +
            2
        }
      }
    }
    dat <- set_attr(dat, "infTime", infTime)
  }
  # For now, sympt status = 1 for all currenly infected
  # NA otherwise
  sympt_vec <- rep(NA, num)
  if (length(idsInf) > 0) {
    sympt_vec[idsInf] <- 1
  }
  dat <- set_attr(dat, "sympt", sympt_vec)

  # Return dat
  dat
}
