# Load test network fit object
fit <- readRDS(test_path("input", "test_nw.RDS"))

# Create a version of the fit object with a missing attribute for testing (that should fail)
nw_missing_attr <- fit$newnetwork
network::delete.vertex.attribute(nw_missing_attr, "age")
fit_with_missing_attr <- fit
fit_with_missing_attr$newnetwork <- nw_missing_attr

# Created a version with an additional attribute (that should pass)
nw_additional_attr <- fit$newnetwork
network::set.vertex.attribute(nw_additional_attr, "new_attr", rep(1, network::network.size(nw_additional_attr)))
fit_with_additional_attr <- fit
fit_with_additional_attr$newnetwork <- nw_additional_attr

# Prep input functions
params <- EpiModel::param.net()

## Initial Conditions
inits <- EpiModel::init.net(i.num = 5)

## Control Settings & Modules
controls <- EpiModel::control.net(
  nsims = 1, nsteps = 1,
  initialize.FUN = mod_sti_initialize,
  verbose = FALSE
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

  # Test initialization with additional attribute (should pass)
  expect_no_error(
    EpiModel::netsim(fit_with_additional_attr, params, inits, controls)
  )
})
