test_that("logpEC", {

  if (requireNamespace("DAISIE")) {
    data("Galapagos_datalist", package = "DAISIE")
    datalist <- Galapagos_datalist

    i <- 4
    phy <- DDD::brts2phylo(datalist[[i]]$branching_times[-c(1, 2)])
    brts <- datalist[[i]]$branching_times
    traits <- sample(c(0, 0), length(brts), replace = TRUE)

    sampling_fraction <- sample(c(1, 1), length(1), replace = TRUE)

    parameter <- list(2.546591, 2.678781, 0.009326754, 1.008583,
                      matrix(c(0), nrow = 1), 0, 1)

    res1 <-  DAISIE_DE_trait_logpEC(
      brts                    = brts,
      phy                     = phy,
      traits                  = traits,
      status                  = 2,
      sampling_fraction       = sampling_fraction,
      parameter               = parameter,
      trait_mainland_ancestor = NA,
      num_observed_states     = 1,
      num_hidden_states       = 1,
      atol                    = 1e-10,
      rtol                    = 1e-10,
      methode                 = "ode45")

    res2 <-  DAISIE:::DAISIE_loglik_CS_choice(
                      pars1 = c(2.546591, 2.678781, Inf, 0.009326754, 1.008583),
                      pars2 = c(100, 11, 0, 2),
                      brts = brts,
                      stac = 2,
                      missnumspec = 0,
                      datalist = datalist)
    testthat::expect_equal(res1, res2)

    res3 <-  DAISIE_DE_trait_logpEC(
      brts                    = brts,
      phy                     = phy,
      traits                  = traits,
      status                  = 2,
      sampling_fraction       = sampling_fraction,
      parameter               = parameter,
      num_observed_states     = 1,
      num_hidden_states       = 1,
      atol                    = 1e-10,
      rtol                    = 1e-10,
      methode                 = "ode45",
      use_Rcpp                = 0)

    testthat::expect_equal(res1, res3)

    res3 <-  DAISIE_DE_trait_logpEC(
      brts                    = brts,
      phy                     = phy,
      traits                  = traits,
      status                  = 2,
      sampling_fraction       = sampling_fraction,
      parameter               = parameter,
      num_observed_states     = 1,
      num_hidden_states       = 1,
      atol                    = 1e-10,
      rtol                    = 1e-10,
      methode                 = "ode45",
      use_Rcpp                = 2)

    testthat::expect_equal(res1, res3)
  }
})
