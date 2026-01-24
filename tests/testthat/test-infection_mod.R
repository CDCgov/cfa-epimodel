# Load test network fit object
fit <- readRDS(test_path("input", "test_nw.RDS"))

ag <- fit$newnetwork %v% "age_group"
ngrps <- length(unique(ag))

# Prep input functions
## 1 time step = 1 year, to speed aging processes
## Parameters
params_single_rate <- EpiModel::param.net(
  inf_prob_mtf = 0.5,
  inf_prob_ftm = 0.5,
  acute_inf_modifier = 2,
  acute_duration = 5,
  act_rate_vec = 2, # single value, should work regardless of n age groups
  cond_prob_vec = 0,
  cond_eff = 0
)
params_rate_per_group <- EpiModel::param.net(
  inf_prob_mtf = 0.5,
  inf_prob_ftm = 0.5,
  acute_inf_modifier = 2,
  acute_duration = 5,
  act_rate_vec = rep(2, ngrps), # correct length
  cond_prob_vec = 0,
  cond_eff = 0
)
params_rate_length_long <- EpiModel::param.net(
  inf_prob_mtf = 0.5,
  inf_prob_ftm = 0.5,
  acute_inf_modifier = 2,
  acute_duration = 5,
  act_rate_vec = rep(2, ngrps + 3), # too long (should be length 1 or length of age groups)
  cond_prob_vec = 0,
  cond_eff = 0
)
params_rate_length_short <- EpiModel::param.net(
  inf_prob_mtf = 0.5,
  inf_prob_ftm = 0.5,
  acute_inf_modifier = 2,
  acute_duration = 5,
  act_rate_vec = rep(2, ngrps - 1), # too short (should be length 1 or length of age groups)
  cond_prob_vec = 0,
  cond_eff = 0
)

## Initial Conditions
inits <- EpiModel::init.net(i.num = 5)

## Control Settings & Modules
controls <- EpiModel::control.net(
  nsims = 1, nsteps = 10,
  infection.FUN = mod_infection,
  epi.by = "female",
  save.other = c("attr"),
  verbose = FALSE
)

test_that("mod_infection works with single act_rate_vec value", {
  expect_no_error(
    sim <- EpiModel::netsim(fit, params_single_rate, inits, controls)
  )
  # And double check that infections occurred
  init_inum <- sim$epi$i.num[[1]][1]
  expect_gt(sum(sim$epi$i.num[[1]]), init_inum)
})

test_that("mod_infection works with act_rate_vec value per age group", {
  expect_no_error(
    sim <- EpiModel::netsim(fit, params_rate_per_group, inits, controls)
  )
  # And double check that infections occurred
  init_inum <- sim$epi$i.num[[1]][1]
  expect_gt(sum(sim$epi$i.num[[1]]), init_inum)
})

test_that("mod_infection errors with wrong length act_rate_vec", {
  expect_error(
    suppressMessages(EpiModel::netsim(fit, params_rate_length_long, inits, controls)),
    regexp = "act_rate_vec parameter length must be either 1 or equal to the number of age groups in the population"
  )

  expect_error(
    suppressMessages(EpiModel::netsim(fit, params_rate_length_short, inits, controls)),
    regexp = "act_rate_vec parameter length must be either 1 or equal to the number of age groups in the population"
  )
})
