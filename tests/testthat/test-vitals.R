# generate nw and fit object for testing
# current required attributes: age, olderpartner, female, race
# number of edges and duration are arbitrary
size <- 10
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

test_that("mod_aging updates age and age_group correctly, arrivalType = departures returns static pop size", {
  # 1 time step = 1 year, to speed aging processes
  params <- EpiModel::param.net(
    units_per_year = 1,
    exitAge = 50, entryAge = 15, arrivalType = "departures",
    entryFemaleProb = 0.5, entryRaceNames = c("A", "B"),
    entryRaceProbs = c(0.6, 0.4)
  )
  inits <- EpiModel::init.net(i.num = 5)
  controls <- EpiModel::control.net(
    nsims = 1, nsteps = 10,
    aging.FUN = mod_aging,
    arrivals.FUN = mod_arrivals,
    departures.FUN = mod_departures,
    save.other = c("attr")
  )

  sim <- EpiModel::netsim(fit, params, inits, controls)

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

  # Test that population size is static when arrivalType = "departures"
  ## when arrivalType = "departures", n new arrivals should equal n departures at each time step
  pop_sizes <- sim$epi$num$sim1
  expect_true(all(pop_sizes == pop_sizes[1]))
})
