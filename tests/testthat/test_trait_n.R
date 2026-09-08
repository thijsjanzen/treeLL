# ---------------------------------------------------------------------------
# A single trait state: one observed and one hidden state, so the rates are
# scalars and there are no transitions. That is the DAISIE model itself.
# ---------------------------------------------------------------------------
single_state_parameter <- list(2.546591, 2.678781, 0.009326754, 1.008583,
                               matrix(0, nrow = 1), 0)
single_state_pars1 <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

# n-sampling likelihood of an endemic singleton
single_state_n <- function(brts, missnumspec, datalist) {
  DAISIE_DE_trait_n(
    DAISIE_DE_trait_function = DAISIE_DE_trait_logpES,
    brts                     = brts,
    missnumspec              = missnumspec,
    traits                   = 0,
    status                   = 2,
    parameter                = single_state_parameter,
    num_observed_states      = 1,
    num_hidden_states        = 1,
    datalist                 = datalist,
    trait_mainland_ancestor  = 1,
    atol                     = 1e-15,
    rtol                     = 1e-15,
    methode                  = "ode45",
    use_Rcpp                 = 2)$loglik
}

# rho-sampling likelihood of the same clade
single_state_rho <- function(brts, sampling_fraction, datalist) {
  DAISIE_DE_trait_logpES(
    datalist                = datalist,
    brts                    = brts,
    trait                   = 0,
    status                  = 2,
    sampling_fraction       = sampling_fraction,
    parameter               = single_state_parameter,
    trait_mainland_ancestor = 1,
    num_observed_states     = 1,
    num_hidden_states       = 1,
    atol                    = 1e-15,
    rtol                    = 1e-15,
    methode                 = "ode45",
    use_Rcpp                = 2)$loglik
}


test_that("trait_n equals the rho-sampling likelihood for complete data", {

  # Complete data means no missing species and full sampling, and then the two
  # are the same likelihood.
  if (requireNamespace("DAISIE")) {
    data("Galapagos_datalist", package = "DAISIE")
    datalist <- Galapagos_datalist
    brts <- datalist[[9]]$branching_times

    testthat::expect_equal(single_state_n(brts, missnumspec = 0, datalist),
                           single_state_rho(brts, sampling_fraction = 1,
                                            datalist))
  }
})


test_that("trait_n equals DAISIE for a single trait state", {

  # With one trait state both packages compute the same likelihood, for any
  # number of missing species.
  if (requireNamespace("DAISIE")) {
    data("Galapagos_datalist", package = "DAISIE")
    datalist <- Galapagos_datalist
    brts <- datalist[[9]]$branching_times

    for (missnumspec in 0:3) {
      expected <- DAISIE:::DAISIE_loglik_CS_choice(
        pars1       = single_state_pars1,
        pars2       = c(100, 11, 0, 0),
        brts        = brts,
        stac        = 2,
        missnumspec = missnumspec,
        datalist    = datalist)

      # the finite differences get less accurate as missnumspec grows
      testthat::expect_equal(single_state_n(brts, missnumspec, datalist),
                             expected, tolerance = 1e-5)
    }
  }
})
