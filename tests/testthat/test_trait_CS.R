test_that("CS", {

  if (requireNamespace("DAISIE")) {

    if (packageVersion("DAISIE") < "4.6.0") {

      data("Galapagos_datalist", package = "DAISIE")
      datalist <- Galapagos_datalist

      parameter <- list(2.546591, 2.678781, 0.009326754, 1.008583,
                        matrix(c(0), nrow = 1), 0, 1)

      res1 <-  DAISIE_DE_trait_loglik_CS(
        datalist            = datalist,
        parameter           = parameter,
        num_observed_states     = 2,
        num_hidden_states       = 1,
        atol                = 1e-10,
        rtol                = 1e-10,
        methode             = "lsodes",
        use_Rcpp            = 0,
        verbose = TRUE
      )


      res3 <-  DAISIE_DE_trait_loglik_CS(
        datalist            = datalist,
        parameter           = parameter,
        num_observed_states     = 2,
        num_hidden_states       = 1,
        atol                = 1e-10,
        rtol                = 1e-10,
        methode             = "lsodes",
        use_Rcpp = 2
      )
      testthat::expect_equal(res1, res3, tol = 0.005)
    }
  }
})
