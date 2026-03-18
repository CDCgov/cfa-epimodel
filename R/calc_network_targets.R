#' @title Generate Starting Network for ERGM Estimation
#'
#' @description Generates empty network of specified size
#' and adds nodal attributes.
#' Currently these attributes are specific to this project:
#' a sex variable ("female"), race, age group, age
#' @param params a list of parameters with a sublist called "pop"
#' with the size of the network
#' and the named attributes and their distribution in the population
#' @param seed optional, maintains distribution of nodal attributes
#' @param assign_deg_casual default = FALSE, if TRUE, uses additional
#' parameters information in "params" to set
#' a plausible casual degree for each node
#' (can help fit main partnership network when cross-network terms present)
#' @param assign_deg_main default = FALSE, if TRUE, uses additional
#' parameters information in "params" to set
#' a plausible main degree for each node
#' (can help fit casual partnership network when cross-network terms present)
#' when casual network fit first in workflow)
#' @importFrom rlang .data
#' @importFrom stats rbinom rpois
#' @return an object of class network
#'
#' @export
#'

generate_init_network <- function(
  params,
  seed = NULL,
  assign_deg_casual = FALSE,
  assign_deg_main = FALSE,
  olderpartner = FALSE
) {
  if (is.null(seed)) {
    warning("No seed specified")
  } else {
    set.seed(seed)
  }

  # desired population size
  if (!is.numeric(params$pop$size)) {
    stop("Population size must be numeric")
  }
  num <- as.integer(params$pop$size)

  # Generate empty network of desired population size
  nw <- network::network.initialize(num, directed = FALSE)

  # Set up list for recording attributes
  attr_values <- list()

  # Create initial attribute vectors
  ### sex variable
  # ensure that dist sums to 1
  if (sum(params$pop$female$dist) != 1) {
    stop("Distribution of sex attribute must sum to 1")
  }
  female <- rep(params$pop$female$levels, params$pop$female$dist * num)
  attr_values$female <- female

  ### race variable
  race <- sample(EpiModel::apportion_lr(
    num,
    params$pop$race$levels,
    params$pop$race$dist
  ))

  attr_values$race <- race

  ### age variables
  # assign age
  age <- rep(NA, num)

  ages <- seq(
    from = params$pop$age$min,
    to = params$pop$age$max + 0.99,
    by = 1 / 365
  )

  # assign age
  age <- sample(ages, num, replace = TRUE)

  attr_values$age <- age
  attr_values$agesq <- age^2
  attr_values$age_group <- dplyr::case_when(
    age < 50 & age >= 45 ~ 7,
    age < 45 & age >= 40 ~ 6,
    age < 40 & age >= 35 ~ 5,
    age < 35 & age >= 30 ~ 4,
    age < 30 & age >= 25 ~ 3,
    age < 25 & age >= 19 ~ 2,
    age < 19 ~ 1
  )

  # Assign additional degree flags (optional)

  # re-create prediction pop table
  pop <- tibble(
    age = params$prediction_pop$age,
    age_group = params$prediction_pop$age_group,
    race = params$prediction_pop$race
  )

  if (isTRUE(assign_deg_casual) && !is.null(params$casual)) {

    pop$casual <- params$casual$nodefactor

    cas_probs <- attr_values |>
      bind_cols() |>
      mutate(age = floor(age)) |>
      left_join(pop, by = c("age", "race")) |>
      pull(.data$casual)

    deg_casual <- rbinom(num, 1, cas_probs)

    attr_values$deg_casual <- deg_casual
  }
  if (isTRUE(assign_deg_main) && !is.null(params$main)) {

    pop$main <- params$main$nodefactor

    main_probs <- attr_values |>
      bind_cols() |>
      mutate(age = floor(age)) |>
      left_join(pop, by = c("age", "race")) |>
      pull(main)

    deg_main <- rbinom(num, 1, main_probs)
    attr_values$deg_main <- deg_main
  }
  if (isTRUE(assign_deg_casual) && is.null(params$casual)) {
    warning(
      "assign_deg_casual = TRUE, but there are no casual parameters in yaml.
      Not setting deg_casual in network attributes."
    )
  }
  if (isTRUE(assign_deg_main) && is.null(params$main)) {
    warning(
      "assign_deg_main = TRUE, but there are no main parameters in yaml.
      Not setting deg_main in network attributes."
    )
  }

  # Set attributes on network
  nw <- network::set.vertex.attribute(nw, names(attr_values), attr_values)
}

#' @title Calculate Network Target Statistics
#'
#' @description Calculates network targets given, at minimum,
#' a starting network, a list of parameters calculated from data,
#' plus desired ERGM term and attribute.
#' Function currently supports calculations for edges, nodefactor,
#' nodematch, absdiff, and concurrent ERGM terms.
#' @param nw the network object outputted from "generate_init_network"
#' (usually already in environment during workflow)
#' @param params the parameter list object with parameters information for
#' all networks estimated from data
#' @param rel string, which relationship/network targets are calculated
#' @param count_type string, which ERGM term to generate targets
#' Currently can only be in the form:
#' "edges", "nodefactor", "nodematch", "absdiff_sqrt_age", or "concurrent".
#' @param attr_name string, used for nodefactor and nodematch targets.
#' @param diff default = FALSE, if TRUE, calc nodematch targets among each group
#' @param inst_correct default = FALSE, if TRUE, adjust year-long
#' cumulative reporting of one-time partnerships to daily or weekly counts
#' @param time_unit default = "weeks",
#' the desired time unit for inst reporting conversion
#' @param level additional statification for nodedov target calculation function
#' @param attr_squared for nodecov target calculation,
#' should squared version of attribute be calculated? (usually, age)
#' @importFrom rlang .data
#' @export
#'

calc_targets <- function(
  nw,
  params,
  rel,
  count_type,
  attr_name = NULL,
  diff = FALSE,
  inst_correct = FALSE,
  level = NULL,
  attr_squared = FALSE,
  time_unit = "weeks"
) {
  # Test that conditions are met to calculate specified target -----------------
  check_conditions(nw, params, rel, count_type, attr_name)

  # Calculate Targets -------------------------------
  ## (currently assumes we use joint distribution to make edges calculations)
  ## get pop size, attribute vectors
  num <- network::network.size(nw)
  attrs <- get_nw_attr_vecs(nw)

  # calculate expected nodefactor targets
  nf_counts <- calc_nodefactor(
    params,
    attrs,
    rel,
    grouping_vars = attr_name
  )

  # calculate expected total edges based on sum of nodefactor
  edges <- calc_edges(
    params,
    attrs,
    rel,
    grouping_vars = NULL
  )

  if (count_type == "edges") {
    final_targets <- edges
  }

  if (count_type == "nodecov") {
    final_targets <- calc_nodecov_age(
      params,
      rel,
      attr_name,
      edges,
      level,
      nf_counts,
      attr_squared
    )
  }

  # if true, calc target for absdiff sqrt age
  if (count_type == "absdiff_sqrt_age") {
    final_targets <- calc_absdiff(params, rel, count_type, edges)
  }

  if (count_type == "concurrent") {
    # calc number of people in rels
    # number of people with > 1 partners
    final_targets <- calc_concurrent(params, rel, num)
  }

  if (count_type == "cross_network") {
    # calc number of people in rels
    # number of people with > 1 partners
    final_targets <- calc_cross_network(params, rel)
  }

  # if true, calc nodefactor and then if true nodematch
  if (count_type %in% c("nodefactor", "nodematch")) {
    # targets for full joint distribution
    attr_targets <- as.numeric(nf_counts)

    # if nodefactor, leave targets as-is
    if (count_type == "nodefactor") {
      final_targets <- attr_targets
    }
    # if nodematch, use nodefactor targets using nodefactor info
    if (count_type == "nodematch") {
      final_targets <- calc_nodematch(
        params,
        attr_name,
        attr_targets,
        rel,
        diff
      )
    }
  }

  # Check that these are reasonable targets
  check_targets(edges, final_targets, count_type)

  # Instantaneous Rel Correction
  ## correct for survey data reflecting number of one-times in last year
  if (rel == "inst" && isTRUE(inst_correct)) {
    final_targets <- inst_correction(final_targets, time_unit)
  }

  # Output targets
  round(final_targets)
}

#' @title Extract Nodal Attribute Vectors from Network
#'
#' @description Outputs nodal attribute vectors as list from network object
#'
#' @param nw a network object,
#' usually the network generated from "generate_init_network"
#'
#' @importFrom network list.vertex.attributes %v%
#' @return A list of nodal attribute vectors
#' @export
#'

get_nw_attr_vecs <- function(nw) {
  if (!"network" %in% class(nw)) {
    stop("input must be a network object")
  }

  n <- list.vertex.attributes(nw)

  attrs <- list()

  for (i in n) {
    if (i == "age") {
      attrs[[i]] <- floor(nw %v% i)
    } else {
      attrs[[i]] <- nw %v% i
    }
  }

  return(attrs) # nolint
}

#' @title Make Empirical Mixing Matrix Symmetrical
#'
#' @description Convert empirical mixing matrix to symmetrical
#' matrix based on mean of upper/lower sections
#'
#' @param mat empirical mixing matrix
#'
#' @return matrix
#' @export

matrix_symmetrical <- function(mat) {
  if (dim(mat)[1] != dim(mat)[2]) {
    stop("Matrix must be square.")
  }

  ncats <- dim(mat)[1]

  newmat <- matrix(NA, ncats, ncats)
  for (i in 1:ncats) {
    for (j in 1:ncats) {
      newmat[i, j] <- mean(c(mat[i, j], mat[j, i]))
    }
  }

  return(newmat) # nolint
}

#' @title Target correction for instantaneous network
#'
#' @description In NSFG (and other surveys), information about
#' the frequency of one-time (instantaneous)
#' partnerships are reported as the cumulative number over the
#' last 12 months. This function converts the
#' empirical yearly target to per-week or per-day counts.
#'
#' @param targets a numeric vector of the ergm target statistics
#' estimated from cumulative year data
#' @param time_unit either "days" or "weeks", the unit of time
#' represented by each discrete step in simulation model
#'
#' @return A numeric vector
#' @export

inst_correction <- function(targets, time_unit = NULL) {
  # if available time unit not specifed, return targets unmodified
  # assumes default target time_unit is a year
  if (is.null(time_unit) || !time_unit %in% c("weeks", "days")) {
    warning("Specified time_unit not available, returning unmodified targets.")
    unit_correction <- 1
  } else {
    if (time_unit == "weeks") {
      unit_correction <- 52
    } else {
      if (time_unit == "days") {
        unit_correction <- 365
      }
    }
  }

  return(targets / unit_correction) # nolint
}

#' @title Target Stats Calculation Helpers
#'
#' @description Small helper functions used as part of calc_targets()
#'
#' @inheritParams calc_targets
#' @param num the number of nodes in the network (population size)
#'
#'
#'
#' @name targets
NULL

#' @rdname targets
#' @param attrs list of network attributes
#' @param grouping_vars attributes to group results by, NULL for all edges
#' @importFrom tibble tibble
#' @importFrom dplyr group_by summarize bind_cols left_join n pick
#' @importFrom rlang .data
#' @export
calc_nodefactor <- function(params, nw_attrs, rel, grouping_vars = NULL) {
  # Re-create prediction pop table
  pop <- tibble(
    age = params$prediction_pop$age,
    agesq = params$prediction_pop$agesq,
    age_group = params$prediction_pop$age_group,
    race = params$prediction_pop$race,
    deg = params[[rel]][["nodefactor"]]
  )

  # Convert nw attrs list to df, add nf probs, calculate by group
  nf <- nw_attrs |>
    bind_cols() |>
    group_by(.data$age, .data$race) |>
    summarize(pop_counts = n(), .groups = "drop") |>
    left_join(pop, by = c("age", "race")) |>
    mutate(rel_counts = .data$pop_counts * .data$deg) |>
    group_by(pick(!!grouping_vars)) |>
    summarize(rels = sum(.data$rel_counts)) |>
    pull(.data$rels) |>
    as.numeric()

  # Return
  nf
}

#' @rdname targets
#' @param nf_targets calculated from calc_single_attr_nodefactor()
#' @export
calc_nodematch <- function(params, attr_name, nf_targets, rel, diff) {
  if (!attr_name %in% names(params[[rel]][["nodematch"]])) {
    stop("Attr name not available for nodematch statistic.")
  }
  attr_probs_nodematch <- params[[rel]][["nodematch"]][[attr_name]]
  if (diff) {
    final_targets <- attr_probs_nodematch$each * nf_targets / 2
  } else {
    # if diff = FALSE, then use the same attr_probs_nodematch for all edges
    final_targets <- attr_probs_nodematch$global * sum(nf_targets) / 2
  }
  final_targets
}

#' @rdname targets
#' @export
calc_edges <- function(params, nw_attrs, rel, grouping_vars = NULL) {
  calc_nodefactor(params, nw_attrs, rel, grouping_vars) / 2
}

#' @rdname targets
#' @param edges output from calc_edges()
#' @export
calc_absdiff <- function(params, rel, count_type, edges) {
  avg <- params[[rel]][[count_type]]
  avg * edges
}

#' @rdname targets
#' @export
calc_concurrent <- function(params, rel, num) {
  num * params[[rel]][["concurrent"]]
}

#' @rdname targets
#' @export
calc_cross_network <- function(params, rel) {
  params$pop$size * params[[rel]][["cross_network"]]
}

#' @rdname targets
#' @param nf_counts output from calc_joint_nodefactor()
#' @param edges output from calc_edges()
#' @export
calc_nodecov_age <- function(
  params,
  rel,
  attr_name,
  edges,
  level = NULL,
  nf_counts,
  attr_squared
) {
  # first check if cutoff exists
  cutoff <- params[[rel]][["nodecov"]][["cutoff"]]

  if (is.null(cutoff)) {
    # get mean of value for attr(i) + attr(j) from data
    if (attr_squared) {
      attr_name <- paste0(attr_name, "sq")
    }
    attr_mean <- params[[rel]][["nodecov"]][[attr_name]]
  }

  if (!is.null(cutoff)) {
    # Need to re-calculate edges for that cutoff group
    nf_counts <- calc_nodefactor(
      params,
      grouping_vars = "age",
      nf_counts
    )
    if (level == "low") {
      age_range <- params$pop$age$min:(cutoff - 1)
    }
    if (level == "high") {
      age_range <- cutoff:params$pop$age$max
    }
    counts <- nf_counts[which(names(nf_counts) %in% age_range)]
    edges <- sum(counts) / 2

    # if cutoff exists, need level to clarify nodecov target
    # get mean of value for attr(i) + attr(j) from data
    if (attr_squared) {
      attr_name <- paste0(attr_name, "sq")
    }
    attr_mean <- params[[rel]][["nodecov"]][[level]][[attr_name]]
  }

  # nodecov target is the attr_mean multiplied by the number of edges
  attr_mean * edges
}

#' @rdname targets
#' @param edges output from calc_edges()
#' @param final_targets vector, final ergm term targets
#' to be checked before output
#' @param threshold default = 0.01, proportion of expected
#' activity based on edges that
#' calulated target is allowed within (+/- threshold)
#' @export
check_targets <- function(edges, final_targets, count_type, threshold = 0.01) {
  # check that all targets are positive
  if (sum(final_targets < 0) > 0) {
    stop("All targets must be >= 0")
  }

  expected_activity <- edges * 2

  if (count_type == "nodefactor") {
    high_threshold <- expected_activity * (1 + threshold)
    low_threshold <- expected_activity * (1 - threshold)
    if (
      (sum(final_targets) > high_threshold) ||
        (sum(final_targets) < low_threshold)
    ) {
      stop(
        "Sum of nodefactor targets do not match expected activity,
        check distribiton of attribute and activity levels of attribute"
      )
    }
  }

  if (count_type == "nodematch") {
    if (sum(final_targets) > expected_activity) {
      stop("Sum of nodematch targets cannot exceed expected activity level")
    }
  }

  if (
    count_type %in%
      c("concurrent", "absdiff_sqrt_age", "nodecov", "cross_network")
  ) {
    if (length(final_targets) > 1) {
      stop(
        "target must be of length 1 for nodecov, concurrent,
        cross_network and absdiff_sqrt_age, targets"
      )
    }
  }
}

#' @rdname targets
#' @export
# nolint start
check_conditions <- function(
  nw,
  params,
  rel,
  count_type,
  attr_name
) {
  if (!"network" %in% class(nw)) {
    stop("inputted initial network must be a network object")
  }
  if (!"list" %in% class(params)) {
    stop("Parameters must be in the form of a list")
  }
  if (!rel %in% names(params)) {
    stop("Specified relationship type does not appear in parameters list")
  }
  if (
    !count_type %in%
      c(
        "edges",
        "nodefactor",
        "nodematch",
        "absdiff_sqrt_age",
        "concurrent",
        "nodecov",
        "cross_network"
      )
  ) {
    stop(
      "This function is only designed to estimate targets
                for edges, nodefactor, nodematch, absdiff by sqrt age, nodecov, cross_network or concurrent ergm terms"
    )
  }
  if (count_type != "edges" && !count_type %in% names(params[[rel]])) {
    stop(
      "Specified count type does not appear in parameter list for this relationship type"
    )
  }
}
# nolint end
