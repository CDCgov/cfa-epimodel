# generate nw and fit object for testing
# current required attributes: age, olderpartner, female, race
# number of edges and duration are arbitrary
size <- 100
nw <- network::network.initialize(size, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "age", seq_len(size) + 18)
nw <- network::set.vertex.attribute(nw, "olderpartner", rbinom(size, 1, 0.2))
nw <- network::set.vertex.attribute(nw, "female", rbinom(size, 1, 0.5))
nw <- network::set.vertex.attribute(nw, "race", apportion_lr(size, c("A", "B"), c(0.6, 0.4)))
fit <- EpiModel::netest(nw,
  formation = ~edges,
  target.stats = 15,
  coef.diss = dissolution_coefs(~ offset(edges), 10),
  edapprox = TRUE
)

# Prep input functions
## 1 time step = 1 year, to speed aging processes
## Parameters
params <- EpiModel::param.net(
  units_per_year = 1,
  exitAge = 50, entryAge = 15, arrivalType = "departures",
  entryFemaleProb = 0.5, entryRaceNames = c("A", "B"),
  entryRaceProbs = c(0.6, 0.4)
)

## Initial Conditions
inits <- EpiModel::init.net(i.num = 5)

## Control Settings & Modules
controls_aging <- EpiModel::control.net(
  nsims = 1, nsteps = 10,
  aging.FUN = mod_aging,
  save.other = c("attr")
)
controls_all_vitals <- EpiModel::control.net(
  nsims = 1, nsteps = 30, # longer time to allow for aging out and new arrivals
  arrivals.FUN = mod_arrivals,
  departures.FUN = mod_departures,
  aging.FUN = mod_aging,
  save.other = c("attr")
)


test_that("mod_aging updates age and age_group correctly, arrivalType = departures returns static pop size", {
  # Run simulation with only aging module
  sim <- EpiModel::netsim(fit, params, inits, controls_aging)

  # Test age and age_group attributes at end of simulation
  sim_age <- sim$attr$sim1$age
  sim_age_group <- sim$attr$sim1$age_group

  ## compare to expected values
  ## age should increase at 1 / units_per_year per time step
  ## but first time step is at time 0, so only nsteps - 1 increments
  age_compare <- nw %v% "age" + ((1 / sim$param$units_per_year) * (sim$control$nsteps - 1))
  age_group_compare <- dplyr::case_when(
    sim_age < 20 ~ 1,
    sim_age >= 20 & sim_age < 25 ~ 2,
    sim_age >= 25 & sim_age < 30 ~ 3,
    sim_age >= 30 & sim_age < 35 ~ 4,
    sim_age >= 35 & sim_age < 40 ~ 5,
    sim_age >= 40 & sim_age < 45 ~ 6,
    sim_age >= 45 ~ 7
  )

  expect_equal(sim_age_group, age_group_compare)
  expect_equal(sim_age, age_compare)
})

test_that("vital dynamics and arrival attr assignment working", {
  # Run simulation with all vital dynamics modules
  sim <- EpiModel::netsim(fit, params, inits, controls_all_vitals)

  ## First check that there are some arrivals and departures
  total_arrivals <- sum(sim$epi$a.flow$sim1, na.rm = TRUE)
  total_departures <- sum(sim$epi$d.flow$sim1, na.rm = TRUE)
  expect_gt(total_arrivals, 0)
  expect_gt(total_departures, 0)

  # Test that total arrivals equal total departures when arrivalType = "departures"
  # and pop size remains the same over time
  expect_equal(total_arrivals, total_departures)
  ## "num" is the total active population size at each time step
  ## regardless of whether or not departed nodes are removed from tracking & networks
  pop_sizes <- sim$epi$num$sim1
  expect_true(all(pop_sizes == pop_sizes[1]))

  # get indices of arrivals
  arrivals_indices <- which(sim$attr$sim1$entrTime > 1 & sim$attr$sim1$active == 1)

  ## Test that ages of arrivals are within expected range
  arrivals_duration <- sim$control$nsteps - (which(sim$epi$a.flow$sim1 > 0)[1])
  expected_age_range <- sim$param$entryAge:(sim$param$entryAge + (arrivals_duration / sim$param$units_per_year))

  arrivals_ages <- sim$attr$sim1$age[arrivals_indices]

  expect_true(all(arrivals_ages %in% expected_age_range))

  ## test that race apportionment for new arrivals matches specified probabilities
  n_arrivals <- length(arrivals_indices)
  arrivals_races <- sim$attr$sim1$race[arrivals_indices]
  arrivals_races_props <- table(arrivals_races) / n_arrivals

  # expect all races to be represented within some tolerance
  ## Note: with small sample sizes and many categories, this test may fail by chance.
  ## May need to increase size of network and/or nsteps if that occurs frequently.
  tolerance <- 0.1
  for (i in seq_along(sim$param$entryRaceNames)) {
    race <- sim$param$entryRaceNames[i]
    expected_prop <- sim$param$entryRaceProbs[i]
    expect_true(race %in% names(arrivals_races_props))
    expect_equal(arrivals_races_props[[i]], expected_prop, tolerance = tolerance)
  }

  ## Test that sex assignment for new arrivals matches specified probability
  ## Assumes female attr is binary 0/1
  arrivals_sex <- sim$attr$sim1$female[arrivals_indices]
  arrivals_sex_props <- table(arrivals_sex) / n_arrivals
  expected_sex_prop <- c(1 - sim$param$entryFemaleProb, sim$param$entryFemaleProb)
  for (i in seq_along(expected_sex_prop)) {
    expect_equal(arrivals_sex_props[[i]], expected_sex_prop[i], tolerance = tolerance)
  }
})
