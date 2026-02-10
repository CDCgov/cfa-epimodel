# generate nw and fit object for testing
# current required attributes: age, age_group, olderpartner, female, race
# number of edges and duration are arbitrary
size <- 100
age_min <- 15
age_max <- 50
age_vec <- seq(age_min, age_max, length.out = size)
age_group_vec <- dplyr::case_when(
  age_vec < 20 ~ 1,
  age_vec >= 20 & age_vec < 25 ~ 2,
  age_vec >= 25 & age_vec < 30 ~ 3,
  age_vec >= 30 & age_vec < 35 ~ 4,
  age_vec >= 35 & age_vec < 40 ~ 5,
  age_vec >= 40 & age_vec < 45 ~ 6,
  age_vec >= 45 ~ 7
)
nw <- network::network.initialize(size, directed = FALSE)
nw <- network::set.vertex.attribute(nw, "age", age_vec)
nw <- network::set.vertex.attribute(nw, "age_group", age_group_vec)
nw <- network::set.vertex.attribute(nw, "olderpartner", rbinom(size, 1, 0.2))
nw <- network::set.vertex.attribute(nw, "female", rbinom(size, 1, 0.5))
nw <- network::set.vertex.attribute(nw, "race", EpiModel::apportion_lr(size, c("A", "B"), c(0.6, 0.4)))
fit <- EpiModel::netest(nw,
  formation = ~edges,
  target.stats = 15,
  coef.diss = EpiModel::dissolution_coefs(~ offset(edges), 10),
  constraints = ~ blocks(attr = ~female, levels2 = diag(TRUE, 2)),
  edapprox = TRUE
)

saveRDS(fit, file = "tests/testthat/input/test_nw.RDS")
