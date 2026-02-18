# Load test network fit object
fit <- readRDS(test_path("input", "test_nw.RDS"))

# --------------------------------------------------------------
# TESTING AGE-SPECIFIC ACT RATES IN mod_infection --------------
#--------------------------------------------------------------
ag <- fit$newnetwork %v% "age_group"
ngrps <- length(unique(ag))

## 1 time step = 1 year, to speed aging processes
## Parameters
static_params <- list(
  inf_prob_mtf = 0.5,
  inf_prob_ftm = 0.5,
  acute_inf_modifier = 2,
  acute_duration = 5,
  cond_prob_vec = 0,
  cond_eff = 0
)
params_single_rate <- with(
  static_params,
  EpiModel::param.net(
    inf_prob_mtf = inf_prob_mtf,
    inf_prob_ftm = inf_prob_ftm,
    acute_inf_modifier = acute_inf_modifier,
    acute_duration = acute_duration,
    cond_prob_vec = cond_prob_vec,
    cond_eff = cond_eff,
    # single value, should work regardless of n age groups
    act_rate_vec = 2
  )
)

params_rate_per_group <- with(
  static_params,
  EpiModel::param.net(
    inf_prob_mtf = inf_prob_mtf,
    inf_prob_ftm = inf_prob_ftm,
    acute_inf_modifier = acute_inf_modifier,
    acute_duration = acute_duration,
    cond_prob_vec = 0,
    cond_eff = 0,
    # correct length
    act_rate_vec = rep(2, ngrps)
  )
)

params_rate_length_long <- with(
  static_params,
  EpiModel::param.net(
    inf_prob_mtf = inf_prob_mtf,
    inf_prob_ftm = inf_prob_ftm,
    acute_inf_modifier = acute_inf_modifier,
    acute_duration = acute_duration,
    cond_prob_vec = 0,
    cond_eff = 0,
    # too long (should be length 1 or length of age groups)
    act_rate_vec = rep(2, ngrps + 3)
  )
)
params_rate_length_short <- with(
  static_params,
  EpiModel::param.net(
    inf_prob_mtf = inf_prob_mtf,
    inf_prob_ftm = inf_prob_ftm,
    acute_inf_modifier = acute_inf_modifier,
    acute_duration = acute_duration,
    cond_prob_vec = cond_prob_vec,
    cond_eff = cond_eff,
    # too short (should be length 1 or length of age groups)
    act_rate_vec = rep(2, ngrps - 1)
  )
)

## Initial Conditions
inits <- EpiModel::init.net(i.num = 5)

## Control Settings & Modules
controls <- EpiModel::control.net(
  nsims = 1,
  nsteps = 10,
  infection.FUN = mod_infection,
  epi.by = "female",
  save.other = c("attr"),
  verbose = FALSE
)

test_that("mod_infection works with single act_rate_vec value", {
  sim <- EpiModel::netsim(fit, params_single_rate, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  # And double check that infections occurred
  df <- as.data.frame(sim)
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_gt(sum_inc_vec, 0)
})

test_that("mod_infection works with act_rate_vec value per age group", {
  sim <- EpiModel::netsim(fit, params_rate_per_group, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  # And double check that infections occurred
  df <- as.data.frame(sim)
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_gt(sum_inc_vec, 0)
})

test_that("mod_infection errors with wrong length act_rate_vec", {
  expect_error(
    suppressMessages(EpiModel::netsim(
      fit,
      params_rate_length_long,
      inits,
      controls
    )),
    regexp = "act_rate_vec parameter length must be either 1 or equal to the number of age groups in the population"
  )

  expect_error(
    suppressMessages(EpiModel::netsim(
      fit,
      params_rate_length_short,
      inits,
      controls
    )),
    regexp = "act_rate_vec parameter length must be either 1 or equal to the number of age groups in the population"
  )
})

# --------------------------------------------------------------
# TESTING DIRECTIONAL ACT RATES IN mod_infection --------------
#--------------------------------------------------------------

## Parameters
static_directional_params <- list(
  acute_inf_modifier = 1,
  acute_duration = 5,
  cond_prob_vec = 0,
  cond_eff = 0,
  act_rate_vec = 2
)

param_only_mtf <- with(
  static_directional_params,
  EpiModel::param.net(
    inf_prob_mtf = 1,
    inf_prob_ftm = 0,
    acute_inf_modifier = acute_inf_modifier,
    acute_duration = acute_duration,
    cond_prob_vec = cond_prob_vec,
    cond_eff = cond_eff,
    act_rate_vec = act_rate_vec
  )
)

param_only_ftm <- with(
  static_directional_params,
  EpiModel::param.net(
    inf_prob_mtf = 0,
    inf_prob_ftm = 1,
    acute_inf_modifier = acute_inf_modifier,
    acute_duration = acute_duration,
    cond_prob_vec = cond_prob_vec,
    cond_eff = cond_eff,
    act_rate_vec = act_rate_vec
  )
)

## Initial Conditions
inits <- EpiModel::init.net(i.num = 50)

## Control Settings & Modules
controls <- EpiModel::control.net(
  nsims = 1,
  nsteps = 10,
  infection.FUN = mod_infection,
  epi.by = "female",
  save.other = c("attr"),
  verbose = FALSE
)

test_that("mod_infection works with directional infection probabilities", {
  # Testing MTF directionality --------------------------------------
  sim_mtf <- EpiModel::netsim(fit, param_only_mtf, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  # Convert sim to data frame for easier manipulation (only extracts epi list, not attributes)
  df <- as.data.frame(sim_mtf)
  # Sum of vector of new infections over time should be > 0
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_gt(sum_inc_vec, 0)
  # Sum of vector of new infections among females should = sum of vector of all new infections
  sum_female_inc_vec <- sum(df$si.flow.female1, na.rm = TRUE)
  expect_equal(sum_female_inc_vec, sum_inc_vec)

  # Testing FTM directionality --------------------------------------
  sim_mtf <- EpiModel::netsim(fit, param_only_ftm, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  df <- as.data.frame(sim_ftm)
  # Double check that infections occurred
  # Sum of vector of new infections over time should be > 0
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_gt(sum_inc_vec, 0)
  # Sum of vector of new infections among males should = sum of vector of all new infections
  sum_male_inc_vec <- sum(df$si.flow.female0, na.rm = TRUE)
  expect_equal(sum_male_inc_vec, sum_inc_vec)
})

# --------------------------------------------------------------
# TESTING CONDOM USE AND EFFECTIVENESS IN mod_infection --------------
#--------------------------------------------------------------
## Static parameters
static_params_cond <- list(
  inf_prob_mtf = 1,
  inf_prob_ftm = 1,
  act_rate_vec = 2,
  acute_inf_modifier = 1,
  acute_duration = 5
)

## Initial Conditions
inits <- EpiModel::init.net(i.num = 50)

## Control Settings & Modules
controls <- EpiModel::control.net(
  nsims = 1,
  nsteps = 10,
  infection.FUN = mod_infection,
  epi.by = "female",
  save.other = c("attr"),
  verbose = FALSE
)

test_that("when condom use & effectiveness = 1, we get no transmissions", {
  params_condom_no_trans <- with(
    static_params_cond,
    EpiModel::param.net(
      cond_prob_vec = 1,
      cond_eff = 1,
      inf_prob_mtf = inf_prob_mtf,
      inf_prob_ftm = inf_prob_ftm,
      acute_inf_modifier = acute_inf_modifier,
      acute_duration = acute_duration,
      act_rate_vec = act_rate_vec
    )
  )

  sim <- EpiModel::netsim(fit, params_condom_no_trans, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  df <- as.data.frame(sim)
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_equal(sum_inc_vec, 0)
})

test_that("with condom use = 1 but effectiveness = 0, we get transmissions", {
  params_condom_no_eff <- with(
    static_params_cond,
    EpiModel::param.net(
      cond_prob_vec = 1,
      cond_eff = 0,
      inf_prob_mtf = inf_prob_mtf,
      inf_prob_ftm = inf_prob_ftm,
      acute_inf_modifier = acute_inf_modifier,
      acute_duration = acute_duration,
      act_rate_vec = act_rate_vec
    )
  )

  sim <- EpiModel::netsim(fit, params_condom_no_eff, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  df <- as.data.frame(sim)
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_gt(sum_inc_vec, 0)
})

test_that("with condom use = 0 but effectiveness = 1, we get transmissions", {
  params_condom_no_cond <- with(
    static_params_cond,
    EpiModel::param.net(
      cond_prob_vec = 0,
      cond_eff = 1,
      inf_prob_mtf = inf_prob_mtf,
      inf_prob_ftm = inf_prob_ftm,
      acute_inf_modifier = acute_inf_modifier,
      acute_duration = acute_duration,
      act_rate_vec = act_rate_vec
    )
  )

  sim <- EpiModel::netsim(fit, params_condom_no_cond, inits, controls) |>
    suppressMessages() |>
    expect_no_error()

  df <- as.data.frame(sim)
  sum_inc_vec <- sum(df$si.flow, na.rm = TRUE)
  expect_gt(sum_inc_vec, 0)
})

test_that("invalid parameter values throw errors", {
  params_condom_invalid <- with(
    static_params_cond,
    EpiModel::param.net(
      cond_prob_vec = 1.5,
      cond_eff = -0.2,
      inf_prob_mtf = inf_prob_mtf,
      inf_prob_ftm = inf_prob_ftm,
      acute_inf_modifier = acute_inf_modifier,
      acute_duration = acute_duration,
      act_rate_vec = act_rate_vec
    )
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

  params_condom_invalid_length <- with(
    static_params_cond,
    EpiModel::param.net(
      cond_prob_vec = c(1, 1, 1),
      cond_eff = 0,
      inf_prob_mtf = inf_prob_mtf,
      inf_prob_ftm = inf_prob_ftm,
      acute_inf_modifier = acute_inf_modifier,
      acute_duration = acute_duration,
      act_rate_vec = act_rate_vec
    )
  )

  expect_error(
    suppressMessages(EpiModel::netsim(
      fit,
      params_condom_invalid_length,
      inits,
      controls
    )),
    regexp = "cond_prob_vec parameter length must be either 1 or equal to the number of networks in the simulation."
  )
})

#--------------------------------------------------------------
# TESTING ACUTE DURATION INF MODIFER IN mod_infection --------------
#--------------------------------------------------------------

## Static parameters
static_params_acute <- list(
  inf_prob_mtf = 1,
  inf_prob_ftm = 1,
  act_rate_vec = 2,
  cond_prob_vec = 0,
  cond_eff = 0
)

## Initial Conditions
inits <- EpiModel::init.net(i.num = 50)

## Control Settings & Modules
controls <- EpiModel::control.net(
  nsims = 1,
  nsteps = 10,
  infection.FUN = mod_infection,
  epi.by = "female",
  save.other = c("attr"),
  verbose = FALSE
)

test_that("acute infection modifier throws error when too large", {
  params_acute_mod_invalid <- with(
    static_params_acute,
    EpiModel::param.net(
      acute_inf_modifier = 5,
      acute_duration = 5,
      cond_prob_vec = cond_prob_vec,
      cond_eff = cond_eff,
      inf_prob_mtf = inf_prob_mtf,
      inf_prob_ftm = inf_prob_ftm,
      act_rate_vec = act_rate_vec
    )
  )

  expect_error(
    sim <- suppressMessages(EpiModel::netsim(
      fit,
      params_acute_mod_invalid,
      inits,
      controls
    )),
    regexp = "Infection probabilities during acute infection stage exceed 1.
      Adjust inf_prob_mtf, inf_prob_ftm, or acute_inf_modifier parameters."
  )
})

test_that("acute infection modifier throws error when < 1", {
  params_acute_mod_invalid2 <- with(
    static_params_acute,
    EpiModel::param.net(
      acute_inf_modifier = 0.5,
      acute_duration = 5,
      cond_prob_vec = cond_prob_vec,
      cond_eff = cond_eff,
      inf_prob_mtf = inf_prob_mtf,
      inf_prob_ftm = inf_prob_ftm,
      act_rate_vec = act_rate_vec
    )
  )

  expect_error(
    sim <- suppressMessages(EpiModel::netsim(
      fit,
      params_acute_mod_invalid2,
      inits,
      controls
    )),
    regexp = "acute_inf_modifier parameter must be >= 1."
  )
})
