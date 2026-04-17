test_that("Calculate target stats", {
  # load data & generate init network
  params <- yaml::read_yaml(test_path("input", "full_nw_params_for_test.yaml"))
  nw <- epimodelcfa::generate_init_network(params, seed = 123)

  # Test that we get no errors in this simple case
  expect_no_error(calc_targets(
    nw = nw,
    params = params,
    rel = "main",
    count_type = "edges"
  ))

  # Tests for names that don't exist in parameters / wrong objects
  ## net should be a network object
  expect_error(calc_targets(nw = data.frame()))
  ## x should be a list
  expect_error(calc_targets(nw = nw, params = data.frame()))
  ## needs count_type
  expect_error(calc_targets(nw = nw, params = params, rel = "main"))
  ## unavailable count_type
  expect_error(calc_targets(
    nw = nw,
    params = params,
    rel = "main",
    count_type = "test"
  ))
  ## assumes no concurrent rels in main net
  expect_error(calc_targets(
    nw = nw,
    params = params,
    rel = "main",
    count_type = "concurrent"
  ))
  ## unavailable attr names
  expect_error(calc_targets(
    nw = nw,
    params = params,
    rel = "main",
    count_type = "nodefactor",
    attr_name = "test"
  ))
  ## nm only avail for race or age_group
  expect_error(calc_targets(
    nw = nw,
    params = params,
    rel = "main",
    count_type = "nodematch",
    attr_name = "age",
    diff = TRUE
  ))

  # expect correct length of returned results
  rels <- names(params)[-c(1:2)] # ignore pop/prediction_pop
  for (rel in rels) {
    # edges
    edge_tar <- calc_targets(
      nw = nw,
      params = params,
      rel = rel,
      count_type = "edges"
    )
    expect_equal(length(edge_tar), 1)

    if (!rel %in% "inst") {
      nm_tar <- calc_targets(
        nw = nw,
        params = params,
        rel = rel,
        count_type = "nodematch",
        attr_name = "race",
        diff = TRUE
      )
      nm_data <- params[[rel]][["nodematch"]][["race"]][["each"]]
      expect_equal(length(nm_tar), length(nm_data))
    }
  }

  # concurrent (only used for casual)
  expect_equal(
    length(calc_targets(
      nw = nw,
      params = params,
      rel = "casual",
      count_type = "concurrent"
    )),
    1
  )
})
