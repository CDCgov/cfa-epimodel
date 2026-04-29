# Load test network fit object
fit <- readRDS(test_path("input", "test_nw.RDS"))

# --------------------------------------------------------------
# TESTING AGE-SPECIFIC ACT RATES IN mod_infection --------------
#--------------------------------------------------------------
ag <- fit$newnetwork %v% "age_group"
ngrps <- length(unique(ag))

## 1 time step = 1 year, to speed aging processes
## Set up reused parameters
single_rate_params <- rate_per_group_params <-
  rate_length_long_params <- rate_length_short_params <- list(
    inf_prob_mtf = 1,
    inf_prob_ftm = 1,
    acute_duration = 5,
    cond_prob_vec = 0,
    cond_eff = 0,
    sympt_prob_m = 0.4,
    sympt_prob_f = 0.6,
    sympt_inf_modifier = 1
  )

## Initial Conditions
inits <- EpiModel::init.net(i.num = 50)

## Control Settings & Modules
controls <- EpiModel::control.net(
  nsims = 1,
  nsteps = 10,
  initialize.FUN = mod_sti_initialize,
  infection.FUN = mod_infection,
  epi.by = "female",
  save.other = c("attr"),
  verbose = FALSE
)

test_that("mod_infection works with single act_rate_vec value", {
  single_rate_params$act_rate_vec <- 2
  params_single_rate <- do.call(
    EpiModel::param.net,
    single_rate_params
  )
  sim <- EpiModel::netsim(fit, params_single_rate, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  # And double check that infections occurred
  df <- as.data.frame(sim)
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_gt(sum_inc_vec, 0)
})

test_that("mod_infection works with act_rate_vec value per age group", {
  rate_per_group_params$act_rate_vec <- rep(2, ngrps)
  params_rate_per_group <- do.call(
    EpiModel::param.net,
    rate_per_group_params
  )
  sim <- EpiModel::netsim(fit, params_rate_per_group, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  # And double check that infections occurred
  df <- as.data.frame(sim)
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_gt(sum_inc_vec, 0)
})

test_that("mod_infection errors with wrong length act_rate_vec", {
  rate_length_long_params$act_rate_vec <- rep(2, ngrps + 3)
  params_rate_length_long <- do.call(
    EpiModel::param.net,
    rate_length_long_params
  )
  expect_error(
    suppressMessages(EpiModel::netsim(
      fit,
      params_rate_length_long,
      inits,
      controls
    )),
    regexp = "act_rate_vec parameter length must be either 1 or
      equal to the number of age groups in the population."
  )

  rate_length_short_params$act_rate_vec <- rep(2, ngrps - 1)
  params_rate_length_short <- do.call(
    EpiModel::param.net,
    rate_length_short_params
  )
  expect_error(
    suppressMessages(EpiModel::netsim(
      fit,
      params_rate_length_short,
      inits,
      controls
    )),
    regexp = "act_rate_vec parameter length must be either 1 or
      equal to the number of age groups in the population."
  )
})

# --------------------------------------------------------------
# TESTING DIRECTIONAL ACT RATES IN mod_infection --------------
#--------------------------------------------------------------

## Parameters
only_mtf_params <- only_ftm_params <- list(
  cond_prob_vec = 0,
  cond_eff = 0,
  act_rate_vec = 2,
  sympt_prob_m = 0.4,
  sympt_prob_f = 0.6,
  sympt_inf_modifier = 1
)

only_mtf_params$inf_prob_mtf <- 1
only_mtf_params$inf_prob_ftm <- 0
params_only_mtf <- do.call(EpiModel::param.net, only_mtf_params)

only_ftm_params$inf_prob_mtf <- 0
only_ftm_params$inf_prob_ftm <- 1
params_only_ftm <- do.call(EpiModel::param.net, only_ftm_params)

## Initial Conditions
inits <- EpiModel::init.net(i.num = 50)

## Control Settings & Modules
controls <- EpiModel::control.net(
  nsims = 1,
  nsteps = 10,
  initialize.FUN = mod_sti_initialize,
  infection.FUN = mod_infection,
  epi.by = "female",
  save.other = c("attr"),
  verbose = FALSE
)

test_that("mod_infection works with directional infection probabilities", {
  # Testing MTF directionality --------------------------------------
  sim_mtf <- EpiModel::netsim(fit, params_only_mtf, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  # Convert sim to data frame for easier manipulation
  # (only extracts epi list, not attributes)
  df <- as.data.frame(sim_mtf)
  # Sum of vector of new infections over time should be > 0
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_gt(sum_inc_vec, 0)
  # Sum of vector of new infections among females should = sum
  # of vector of all new infections
  sum_female_inc_vec <- sum(df$si.flow.female1, na.rm = TRUE)
  expect_equal(sum_female_inc_vec, sum_inc_vec)

  # Testing FTM directionality --------------------------------------
  sim_ftm <- EpiModel::netsim(fit, params_only_ftm, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  df <- as.data.frame(sim_ftm)
  # Double check that infections occurred
  # Sum of vector of new infections over time should be > 0
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_gt(sum_inc_vec, 0)
  # Sum of vector of new infections among males should = sum
  # of vector of all new infections
  sum_male_inc_vec <- sum(df$si.flow.female0, na.rm = TRUE)
  expect_equal(sum_male_inc_vec, sum_inc_vec)
})

# --------------------------------------------------------------
# TESTING CONDOM USE AND EFFECTIVENESS IN mod_infection --------------
#--------------------------------------------------------------
## Static parameters
condom_no_trans_params <- condom_no_eff_params <-
  condom_no_cond_params <- condom_invalid_params <-
    condom_invalid_length_params <- list(
      inf_prob_mtf = 1,
      inf_prob_ftm = 1,
      act_rate_vec = 2,
      sympt_prob_m = 0.4,
      sympt_prob_f = 0.6,
      sympt_inf_modifier = 1
    )

## Initial Conditions
inits <- EpiModel::init.net(i.num = 50)

## Control Settings & Modules
controls <- EpiModel::control.net(
  nsims = 1,
  nsteps = 10,
  initialize.FUN = mod_sti_initialize,
  infection.FUN = mod_infection,
  epi.by = "female",
  save.other = c("attr"),
  verbose = FALSE
)

test_that("when condom use & effectiveness = 1, we get no transmissions", {
  condom_no_trans_params$cond_prob_vec <- 1
  condom_no_trans_params$cond_eff <- 1
  params_condom_no_trans <- do.call(
    EpiModel::param.net,
    condom_no_trans_params
  )

  sim <- EpiModel::netsim(fit, params_condom_no_trans, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  df <- as.data.frame(sim)
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_equal(sum_inc_vec, 0)
})

test_that("with condom use = 1 but effectiveness = 0, we get transmissions", {
  condom_no_eff_params$cond_prob_vec <- 1
  condom_no_eff_params$cond_eff <- 0
  params_condom_no_eff <- do.call(
    EpiModel::param.net,
    condom_no_eff_params
  )

  sim <- EpiModel::netsim(fit, params_condom_no_eff, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  df <- as.data.frame(sim)
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_gt(sum_inc_vec, 0)
})

test_that("with condom use = 0 but effectiveness = 1, we get transmissions", {
  condom_no_cond_params$cond_prob_vec <- 0
  condom_no_cond_params$cond_eff <- 1
  params_condom_no_cond <- do.call(
    EpiModel::param.net,
    condom_no_cond_params
  )

  sim <- EpiModel::netsim(fit, params_condom_no_cond, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  df <- as.data.frame(sim)
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_gt(sum_inc_vec, 0)
})

test_that("invalid parameter values throw errors", {
  condom_invalid_params$cond_prob_vec <- 1.5
  condom_invalid_params$cond_eff <- -0.2
  params_condom_invalid <- do.call(
    EpiModel::param.net,
    condom_invalid_params
  )

  expect_error(
    suppressMessages(EpiModel::netsim(
      fit,
      params_condom_invalid,
      inits,
      controls
    )),
    regexp = "All infection-related probabilities must be >=0 and <=1"
  )

  condom_invalid_length_params$cond_prob_vec <- c(1, 1, 1)
  condom_invalid_length_params$cond_eff <- 0
  params_condom_invalid_length <- do.call(
    EpiModel::param.net,
    condom_invalid_length_params
  )

  expect_error(
    suppressMessages(EpiModel::netsim(
      fit,
      params_condom_invalid_length,
      inits,
      controls
    )),
    regexp = "cond_prob_vec parameter length must be either 1 or
      equal to the number of networks in the simulation."
  )
})

#--------------------------------------------------------------
# TESTING SYMPTOMATIC INFECTION MODIFIER IN mod_infection --------------
#--------------------------------------------------------------

## Static parameters
sympt_mod_invalid_params <- sympt_mod_invalid_params2 <- list(
  inf_prob_mtf = 1,
  inf_prob_ftm = 1,
  act_rate_vec = 2,
  cond_prob_vec = 0,
  cond_eff = 0,
  sympt_prob_m = 0.4,
  sympt_prob_f = 0.6
)

## Initial Conditions
inits <- EpiModel::init.net(i.num = 50)

## Control Settings & Modules
controls <- EpiModel::control.net(
  nsims = 1,
  nsteps = 10,
  initialize.FUN = mod_sti_initialize,
  infection.FUN = mod_infection,
  epi.by = "female",
  save.other = c("attr"),
  verbose = FALSE
)

test_that("symptomatic infection modifier throws error when too large", {
  sympt_mod_invalid_params$sympt_inf_modifier <- 5
  params_sympt_mod_invalid <- do.call(
    EpiModel::param.net,
    sympt_mod_invalid_params
  )

  expect_error(
    sim <- suppressMessages(EpiModel::netsim(
      fit,
      params_sympt_mod_invalid,
      inits,
      controls
    )),
    regexp = "Infection probabilities during symptomatic infection will exceed 1.
    Adjust inf_prob_mtf, inf_prob_ftm, or sympt_inf_modifier parameters."
  )
})

test_that("symptomatic infection modifier throws error when < 1", {
  sympt_mod_invalid_params2$sympt_inf_modifier <- 0.5
  params_sympt_mod_invalid2 <- do.call(
    EpiModel::param.net,
    sympt_mod_invalid_params2
  )

  expect_error(
    sim <- suppressMessages(EpiModel::netsim(
      fit,
      params_sympt_mod_invalid2,
      inits,
      controls
    )),
    regexp = "sympt_inf_modifier parameter must be >= 1."
  )
})
