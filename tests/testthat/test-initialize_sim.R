# Load test network fit object
fit <- readRDS(test_path("input", "test_nw.RDS"))

# Create a version of the fit object with a missing attribute for testing
nw_missing_attr <- fit$newnetwork
network::delete.vertex.attribute(nw_missing_attr, "age")
fit_with_missing_attr <- fit
fit_with_missing_attr$newnetwork <- nw_missing_attr

# Prep input functions
params <- EpiModel::param.net()

## Initial Conditions
inits <- EpiModel::init.net(i.num = 5)

## Control Settings & Modules
controls <- EpiModel::control.net(
  nsims = 1, nsteps = 1,
  initialize.FUN = mod_sti_initialize
)

test_that("mod_sti_initialize throws error if required nodal attributes are missing", {
  # Test normal initialization
  expect_no_error(
    EpiModel::netsim(fit, params, inits, controls)
  )

  # Test initialization with missing required attribute
  expect_error(
    EpiModel::netsim(fit_with_missing_attr, params, inits, controls),
    regexp = "The following required attributes are missing from the network at initialization: age"
  )
})
