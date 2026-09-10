# ---------------------------------------------------------------------------
# Requirement 1:
# With no missing species, n-sampling = rho-sampling with rho = 1.
#
# Requirement 2:
# With one trait state, n-sampling = original DAISIE n-sampling.
# ---------------------------------------------------------------------------


parameter <- list(
  c(2.546591, 1.2, 1, 0.2),
  c(2.678781, 2, 1.9, 3),
  c(0.009326754, 0.003, 0.002, 0.2),
  c(1.008583, 1, 2, 1.5),
  matrix(
    c(0, .1, .05, 0,
      .33, 0, 0, .0086,
      .005, 0, 0, .005,
      0, .5, .35, 0),
    nrow = 4,
    byrow = TRUE
  ),
  1
)


test_that("complete n-sampling equals full rho-sampling", {

  brts <- c(23, 13)
  trait <- 0

  # One sampled singleton in state 0; none in state 1
  S <- c(1, 0)

  log_f <- make_log_f_ES(
    brts = brts,
    trait = trait,
    status = 2,
    parameter = parameter,
    num_observed_states = 2,
    num_hidden_states = 2,
    S = S,
    trait_mainland_ancestor = c(1, 0)
  )

  n_result <- DAISIE_DE_trait_n(
    log_f = log_f,
    missnumspec = c(0, 0),
    S = S
  )$loglik

  rho_result <- DAISIE_DE_trait_logpES(
    brts = brts,
    sampling_fraction = c(1, 1),
    trait = trait,
    status = 2,
    parameter = parameter,
    num_observed_states = 2,
    num_hidden_states = 2,
    trait_mainland_ancestor = c(1, 0)
  )$loglik

  expect_equal(n_result, rho_result)
})


test_that("one-state n-sampling equals original DAISIE", {

  data("Galapagos_datalist", package = "DAISIE")

  trait_parameter <- list(
    2.546591,
    2.678781,
    0.009326754,
    1.008583,
    matrix(0, nrow = 1),
    0
  )

  daisie_pars <- c(
    2.546591,
    2.678781,
    Inf,
    0.009326754,
    1.008583
  )

  brts <- Galapagos_datalist[[4]]$branching_times
  traits <- rep(0, length(brts) - 1)
  phy <- DDD::brts2phylo(brts[-c(1, 2)])

  # Every sampled species has the single possible trait state
  S <- c(length(traits))

  log_f <- make_log_f_EC(
    brts = brts,
    traits = traits,
    status = 2,
    parameter = trait_parameter,
    num_observed_states = 1,
    num_hidden_states = 1,
    S = S,
    trait_mainland_ancestor = 1,
    phy = phy,
    num_threads = 1
  )

  for (missing in 0:10) {

    trait_result <- DAISIE_DE_trait_n(
      log_f = log_f,
      missnumspec = missing,
      S = S
    )$loglik

    daisie_result <- DAISIE:::DAISIE_loglik_CS_choice(
      pars1 = daisie_pars,
      pars2 = c(100, 11, 0, 0),
      brts = brts,
      stac = 2,
      missnumspec = missing,
      datalist = Galapagos_datalist
    )

    expect_equal(
      trait_result,
      daisie_result,
      tolerance = 1e-5,
      info = paste("Missing species:", missing)
    )
  }
})
