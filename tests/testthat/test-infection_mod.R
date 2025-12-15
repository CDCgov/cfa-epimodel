# Load test network fit object
fit <- readRDS(test_path("input", "test_nw.RDS"))

ag <- fit$newnetwork %v% "age_group"
ngrps <- length(unique(ag))

# Prep input functions
## 1 time step = 1 year, to speed aging processes
## Parameters
params_single_rate <- EpiModel::param.net(
  infProbMTF = 1,
  infProbFTM = 1,
  act.rate = 2 # single value, should work regardless of n age groups
)
params_rate_per_group <- EpiModel::param.net(
  infProbMTF = 1,
  infProbFTM = 1,
  act.rate = rep(2, ngrps) # correct length
)
params_rate_length_long <- EpiModel::param.net(
  infProbMTF = 1,
  infProbFTM = 1,
  act.rate = rep(2, ngrps + 3) # too long (should be length 1 or length of age groups)
)

params_rate_length_short <- EpiModel::param.net(
  infProbMTF = 1,
  infProbFTM = 1,
  act.rate = rep(2, ngrps - 1) # too short (should be length 1 or length of age groups)
)

## Initial Conditions
inits <- EpiModel::init.net(i.num = 5)

## Control Settings & Modules
controls <- EpiModel::control.net(
  nsims = 1, nsteps = 10,
  infection.FUN = mod_infection,
  save.other = c("attr"),
  verbose = FALSE
)

test_that("mod_infection works with single act.rate value", {
  expect_no_error(
    sim <- EpiModel::netsim(fit, params_single_rate, inits, controls)
  )
  # And double check that infections occurred
  init_inum <- sim$epi$i.num[[1]][1]
  expect_gt(sum(sim$epi$i.num[[1]]), init_inum)
})

test_that("mod_infection works with act.rate value per age group", {
  expect_no_error(
    sim <- EpiModel::netsim(fit, params_rate_per_group, inits, controls)
  )
  # And double check that infections occurred
  init_inum <- sim$epi$i.num[[1]][1]
  expect_gt(sum(sim$epi$i.num[[1]]), init_inum)
})
test_that("mod_infection errors with wrong length act.rate", {
  expect_error(
    suppressMessages(EpiModel::netsim(fit, params_rate_length_long, inits, controls)),
    regexp = "act.rate parameter length must be either 1 or equal to the number of age groups in the population"
  )

  expect_error(
    suppressMessages(EpiModel::netsim(fit, params_rate_length_short, inits, controls)),
    regexp = "act.rate parameter length must be either 1 or equal to the number of age groups in the population"
  )
})
