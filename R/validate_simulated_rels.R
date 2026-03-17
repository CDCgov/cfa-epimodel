#' @title Get Target Mean Degrees for Main and Casual Networks
#' @description Extracts target degrees for main and casual networks
#' from a YAML file. Currently assumes the YAML file
#' includes sexual activity terms under "nodefactor" for each network
#' stratified by a joint set of attributes,
#' kept as vectors under "prediction_pop"
#' @param yaml_params_loc Path to the YAML file containing network targets.
#' @param nets Character vector specifying which networks to extract targets
#' for.
#' Default is c("main", "casual").
#' @return A data frame where each cell represents the mean degree that
#' set of attrs and network type.
#' @importFrom yaml read_yaml
#' @importFrom dplyr mutate group_by summarize
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @export
get_target_degrees <- function(
  yaml_params_loc,
  nets = c("main", "casual")
) {
  # load network targets from yaml file
  x <- read_yaml(yaml_params_loc)

  # Validate inputs, currently only supports main and casual networks
  if (sum(nets == c("main", "casual")) != 2) {
    stop(
      "Currently only 'main' and 'casual' networks are supported"
    )
  }

  # Make sure the specified networks and attributes exist in the YAML file
  if (!all(nets %in% names(x))) {
    stop("Specified networks not found in the YAML file.")
  }

  if (is.null(x$pop)) {
    stop("Population attributes not found in the YAML file.")
  }

  # Extract the attribute values
  dat <- tibble(
    age = x$prediction_pop$age,
    race = x$prediction_pop$race,
    main = x$main$nodefactor,
    casual = x$casual$nodefactor
  )


  # Rreshape to long format
  dat |>
    pivot_longer(
      cols = c("main", "casual"),
      names_to = "network",
      values_to = "degree"
    ) |>
    mutate(data = "target")
}

#' @title Get Edges History from Simulation
#' @description Extracts edges history as the number of edges per time step from
#' a simulation object, calculating the difference between the simulated edges
#' and the target edges for all networks.
#' Assumes edge history is tracked in the simulation object via the "epi" list,
#' not via setting the control parameter 'nwstats' TRUE.
#' @param sim A simulation object of class `EpiModel::netsim`.
#' @return A df with time, simulation id, edges for all networks,
#' the difference from target edges, and the percentage diff from target edges.
#' @importFrom rlang .data := !!
#' @importFrom dplyr select mutate filter group_by ungroup all_of bind_rows
#'             rename
#' @importFrom tidyr pivot_longer
#' @export
get_edges_history <- function(sim, edge_prefix = "edges_net_",
                              net_names = c("main", "casual", "inst")) {
  # Get number of networks tracked in the simulation object
  n_nets <- sim$num.nw

  edges <- paste0(edge_prefix, seq_len(n_nets))

  if (!all(edges %in% names(sim$epi))) {
    stop("Edges history for all networks not found in the simulation object.")
  }

  edges_list <- list()

  for (i in seq_len(n_nets)) {
    # Extract edges target for given network
    # always stored as first element in target.stats for each network
    target_edges <- sim$nwparam[[i]]$target.stats[1]
    edge <- edges[i]
    name <- net_names[i]

    # Extract edges history from simulation object
    edges_list[[i]] <- sim |>
      as.data.frame() |>
      select(.data$time, .data$sim, !!edge) |>
      rename(!!name := !!edge) |>
      pivot_longer(cols = !!name, names_to = "net",
                   values_to = "edges") |>
      mutate(
        target = !!target_edges,
        absolute = .data$edges - .data$target,
        percent = (.data$absolute / .data$target) * 100
      ) |>
      pivot_longer(
        cols = c("absolute", "percent", "edges"),
        names_to = "diff_type",
        values_to = "diff"
      ) |>
      group_by(.data$time, .data$net, .data$diff_type) |>
      mutate(mean = mean(.data$diff, na.rm = TRUE)) |>
      ungroup() |>
      mutate(target = ifelse(.data$diff_type == "edges", .data$target, 0))
  }

  edges_df <- bind_rows(edges_list)

  # Return df with edges history and diffs from target
  edges_df
}

#' @title Plot Edges History
#' @description Plots the edges history for a specified network and type
#' of difference (absolute, percent, or edges).
#' @param x Either a data frame containing the edges history
#' (the output of `get_edges_history()`),
#' or a simulation object of class `EpiModel::netsim`.
#' If a simulation object is provided, edges history will be extracted
#' using `get_edges_history()`.
#' @param network A character string specifying the network type,
#' either "main" or "casual".
#' @param type A character string specifying the type of difference to plot,
#' either "percent", "absolute", or "edges".
#' @return A ggplot object showing the edges history over time for the
#' specified network and type.
#' @importFrom ggplot2 ggplot aes geom_line geom_hline labs theme
#' @importFrom dplyr filter pull
#' @importFrom rlang .data
#' @importFrom viridis scale_color_viridis
#' @export
plot_edges_history <- function(x, network, type) {
  if (!class(x) %in% c("netsim", "data.frame")) {
    stop("x must be a netsim object or a data frame.")
  }
  if (!type %in% c("percent", "absolute", "edges")) {
    stop("type must be one of 'percent', 'absolute', or 'edges'.")
  }
  if (!network %in% c("main", "casual", "inst")) {
    stop("network must be either 'main', 'casual', or 'inst'.")
  }

  if (inherits(x, "netsim")) {
    edges_df <- get_edges_history(x)
  } else {
    edges_df <- x
  }

  if (
    !all(
      c("time", "sim", "net", "target", "diff_type", "diff", "mean") %in%
        names(edges_df)
    )
  ) {
    stop(
      "edges_df must contain the columns: time, sim, net, target, diff_type, diff, and mean."
    )
  }

  target_val <- edges_df |>
    filter(.data$net == !!network, .data$diff_type == !!type) |>
    pull(.data$target) |>
    unique()

  edges_df |>
    filter(.data$net == !!network, .data$diff_type == !!type,
           .data$time >= 2) |> # time 1 = NA, so start at time 2 for plotting
    mutate(sim = as.factor(.data$sim)) |>
    ggplot(aes(x = .data$time, y = .data$diff, color = .data$sim)) +
    geom_line(alpha = 0.5) +
    geom_line(aes(y = .data$mean), color = "black", linewidth = 3) +
    geom_hline(aes(yintercept = target_val)) +
    labs(
      title = paste(
        "Edges history for ",
        network,
        " network (",
        type,
        ")",
        sep = ""
      ),
      y = paste(type, "difference"),
      x = "time"
    ) +
    scale_color_viridis(discrete = TRUE) +
    theme(legend.position = "none")
}

#' @title Summarize Final Degrees from Simulation
#' @description Summarizes the final degrees of individuals in the
#' main and casual networks
#' at the end of the simulation and calculates the mean degree for
#' each age and race combination.
#' For `EpiModel::netest` objects, a network is simulated using the
#' "newnetwork" slot as its basis.
#' For `EpiModel::netdx` objects, data is extracted from the
#' "tedgelist" (timed edgelist) slot.
#' For `EpiModel::netsim` objects, data is extracted from the "network" slot.
#' @param input An object of class `EpiModel::netest`, `EpiModel::netsim`
#' or `EpiModel::netdx`.
#' @param network A character string specifying the network type, either
#' "main" or "casual".
#' @return A data frame summarizing the mean degree, interquartile range
#' (IQR), and data source
#' for each age and race combination
#' @importFrom dplyr group_by summarize mutate
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom EpiModel get_degree
#' @importFrom stats quantile
#' @importFrom network %v% network.size
#' @importFrom ergm control.simulate.formula
#' @importFrom stats simulate
#' @export

# frequency of rels by age in networks at end of simulation
summarize_final_degrees <- function(input, network) {
  if (
    !inherits(input, "netsim") &&
      !inherits(input, "netdx") &&
      !inherits(input, "netest")
  ) {
    stop("input sim must be a netest, netdx, or netsim object.")
  }

  if (inherits(input, "netdx") && is.null(input$tedgelist)) {
    stop(
      "netdx object must contain simulated networks
      using option keep.tedgelist = TRUE."
    )
  }

  if (!network %in% c("main", "casual")) {
    stop("network must be either 'main' or 'casual'.")
  }

  simdat <- NULL

  # Extract final degrees from netsim object
  if (inherits(input, "netsim")) {
    nsims <- input$control$nsims
    data_type <- "simulated"
    net_position <- ifelse(network == "main", 1, 2)
    for (i in seq_len(nsims)) {
      this_sim <- paste0("sim", i)
      d <- data.frame(
        age = floor(input[["network"]][[this_sim]][[net_position]] %v% "age"),
        race = input[["network"]][[this_sim]][[net_position]] %v% "race",
        deg = get_degree(input[["network"]][[this_sim]][[net_position]]),
        sim = this_sim
      )
      simdat <- rbind(simdat, d)
    }
  }
  # Extract final degrees from netdx object
  if (inherits(input, "netdx")) {
    nsims <- input$nsims
    data_type <- "netdx"
    for (i in seq_len(nsims)) {
      # convert timed edgelist to edgelist readable by EpiModel::get_degree()
      el <- input$tedgelist[[i]]
      active_el <- as.matrix(
        el[el$terminus.censored, c("tail", "head")],
        rownames.force = FALSE
      )
      attr(active_el, "n") <- network.size(input$nw)

      d <- data.frame(
        age = floor(input$nw %v% "age"),
        race = input$nw %v% "race",
        deg = get_degree(active_el),
        sim = i
      )
      simdat <- rbind(simdat, d)
    }
  }
  # Extract final degrees from netest object
  if (inherits(input, "netest")) {
    data_type <- "fitted"
    nw <- simulate(
      input$formula,
      coef = input$coef.form.crude,
      basis = input$newnetwork,
      constraints = input$constraints,
      control = control.simulate.formula(
        MCMC.prop = ~sparse,
        MCMC.burnin = 2e+05
      ),
      dynamic = FALSE
    )
    simdat <- data.frame(
      age = floor(nw %v% "age"),
      race = nw %v% "race",
      deg = get_degree(nw),
      sim = 1
    )
  }

  simdat |>
    # calc mean degree by sim, age, race
    group_by(.data$sim, .data$age, .data$race) |>
    summarize(
      mean_deg = mean(.data$deg, na.rm = TRUE),
      .groups = "drop"
    ) |>
    # mean and IQR of mean degree across sims
    group_by(.data$age, .data$race) |>
    summarize(
      degree = mean(.data$mean_deg),
      IQR1 = quantile(.data$mean_deg, 1 / 4),
      IQR3 = quantile(.data$mean_deg, 3 / 4),
      .groups = "drop"
    ) |>
    mutate(data = !!data_type, network = !!network)
}

#' @title Plot Final Degrees for Main and Casual Networks
#' @description Plots the final degrees of individuals in the main
#' and casual networks summarized across simulations
#' and compares them to target degrees extracted from a YAML file.
#' @param input A simulation object of class `EpiModel::netsim` or
#' \`EpiModel::netdx`.
#' @param network A character string specifying the network type, either
#' "main" or "casual".
#' @param yaml_params_loc Path to the YAML file containing target degrees.
#' @return A ggplot object showing the final degrees for the specified
#' network type,
#' comparing simulated degrees to target degrees.
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar facet_wrap
#' @importFrom dplyr filter mutate
#' @importFrom rlang .data
#' @export
plot_final_degrees <- function(input, network, yaml_params_loc) {
  if (!network %in% c("main", "casual")) {
    stop("network must be either 'main' or 'casual'.")
  }

  s <- summarize_final_degrees(input, network)
  t <- get_target_degrees(yaml_params_loc) |>
    # targets do not have IQRs
    dplyr::mutate(IQR1 = .data$degree, IQR3 = .data$degree) |>
    dplyr::filter(.data$network == !!network)

  y <- rbind(s, t)

  y |>
    ggplot(aes(x = .data$age, y = .data$degree, color = .data$data)) +
    geom_point() +
    geom_errorbar(aes(ymin = .data$IQR1, ymax = .data$IQR3), width = 0.2) +
    facet_wrap(~ .data$race)
}


#' @title Get Mean Durations of Relationships at End of Simulation
#' @description Calculates the mean durations of relationships
#' in the main and casual
#' networks at the end
#' of the simulation, comparing them to target durations specified
#' in a YAML file.
#' @param sim A simulation object of class `EpiModel::netsim`.
#' @param nets A character vector specifying the networks to calculate
#' durations for, default is c("main", "casual").
#' @param yaml_params_loc Path to the YAML file containing target durations.
#' @return A data frame summarizing the target and simulated mean durations
#' for each network,
#' along with the standard deviation of the simulated durations.
#' @importFrom yaml read_yaml
#' @importFrom stats sd
#' @importFrom network %n%
#' @export
get_mean_durations <- function(
  sim,
  nets = c("main", "casual"),
  yaml_params_loc
) {
  x <- read_yaml(yaml_params_loc)
  main_durs <- NULL
  casual_durs <- NULL
  nsims <- sim$control$nsims
  nsteps <- sim$control$nsteps

  for (i in seq_len(nsims)) {
    this_sim <- paste0("sim", i)
    m <- sim[["network"]][[this_sim]][[1]] %n% "lasttoggle"
    mdurs <- nsteps - m[, 3]
    main_durs <- c(main_durs, mean(mdurs, na.rm = TRUE))

    c <- sim[["network"]][[this_sim]][[2]] %n% "lasttoggle"
    cdurs <- nsteps - c[, 3]
    casual_durs <- c(casual_durs, mean(cdurs, na.rm = TRUE))
  }

  data.frame(
    network = c("main", "casual"),
    target = c(
      x[[nets[1]]]$duration$overall,
      x[[nets[2]]]$duration$overall
    ),
    sim_mean = c(mean(main_durs), mean(casual_durs)),
    sim_sd = c(sd(main_durs), sd(casual_durs))
  )
}
