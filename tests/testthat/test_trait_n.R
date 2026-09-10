# ---------------------------------------------------------------------------
# Requirement 1:
# With no missing species, n-sampling = rho-sampling with rho = 1.
#
# Requirement 2:
# With one trait state, n-sampling = original DAISIE n-sampling.
# ---------------------------------------------------------------------------


# ---- Common four-state TRAISIE parameters ---------------------------------
parameter <- list(
  c(2.546591, 1.2, 1, 0.2),        # lambda_c
  c(2.678781, 2, 1.9, 3),          # mu
  c(0.009326754, 0.003, 0.002, .2),# gamma
  c(1.008583, 1, 2, 1.5),          # lambda_a
  matrix(c(0, .1, .05, 0,
           .33, 0, 0, .0086,
           .005, 0, 0, .005,
           0, .5, .35, 0),
         nrow = 4, byrow = TRUE),  # transition matrix
  1
)


test_that("complete n-sampling equals full rho-sampling", {

  # Endemic singleton: one species, with trait state 0
  brts <- c(23, 13)
  trait <- 0
  d <- 2
  h <- 2

  # n-sampling: no species are missing
  n_result <- DAISIE_DE_trait_n(
    DAISIE_DE_trait_function = DAISIE_DE_trait_logpES,
    brts = brts,
    missnumspec = c(0, 0),
    traits = trait,
    status = 2,
    parameter = parameter,
    num_observed_states = d,
    num_hidden_states = h,
    trait_mainland_ancestor = c(1, 0)
  )$loglik

  # rho-sampling: all species are sampled
  rho_result <- DAISIE_DE_trait_logpES(
    brts = brts,
    sampling_fraction = c(1, 1),
    trait = trait,
    status = 2,
    parameter = parameter,
    num_observed_states = d,
    num_hidden_states = h,
    trait_mainland_ancestor = c(1, 0)
  )$loglik

  expect_equal(n_result, rho_result)
})


test_that("one-state n-sampling equals original DAISIE", {

  data("Galapagos_datalist", package = "DAISIE")

  # One state means TRAISIE reduces to DAISIE
  trait_parameter <- list(
    2.546591, 2.678781, 0.009326754, 1.008583,
    matrix(0, nrow = 1),
    0
  )

  # Original DAISIE parameters:
  # lambda_c, mu, K = Inf, gamma, lambda_a
  daisie_pars <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  # Use one endemic clade from the Galápagos data
  brts <- Galapagos_datalist[[4]]$branching_times
  traits <- rep(0, length(brts) - 1)
  phy <- DDD::brts2phylo(brts[-c(1, 2)])

  for (missing in 0:10) {

    trait_result <- DAISIE_DE_trait_n(
      DAISIE_DE_trait_function = DAISIE_DE_trait_logpEC,
      brts = brts,
      missnumspec = missing,
      traits = traits,
      status = 2,
      parameter = trait_parameter,
      num_observed_states = 1,
      num_hidden_states = 1,
      datalist = Galapagos_datalist,
      trait_mainland_ancestor = 1,
      phy = phy,
      num_threads = 1
    )$loglik

    daisie_result <- DAISIE:::DAISIE_loglik_CS_choice(
      pars1 = daisie_pars,
      pars2 = c(100, 11, 0, 0),
      brts = brts,
      stac = 2,
      missnumspec = missing,
      datalist = Galapagos_datalist
    )

    expect_equal(trait_result, daisie_result, tolerance = 1e-5)
  }
})
