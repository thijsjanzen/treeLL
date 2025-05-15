test_that("identical R and Rcpp", {
  phy <- ape::read.tree(text = "((t3:0.459,t2:0.459):2.963,(t4:2.884,t1:2.884):0.538);")

  traits <- rep(0, length(phy$tip.label))

  parameter <- list( c(2.546591, 0), c (2.678781, 0), c(0.009326754,0),
                     c(1.008583, 0), c(0.000,0.000), 0 )

  brts <- c(4.000, 3.958, 3.422, 2.884, 0.459)

  loglik_DAISIE_trait <- treeLL::DAISIE_DE_logpEC_trait1(brts = brts,
                                                        missnumspec = 0,
                                                        parameter = parameter,
                                                        phy = phy,
                                                        traits = traits,
                                                        cond = "proper_cond",
                                                        root_state_weight = "proper_weights",
                                                        see_ancestral_states = TRUE,
                                                        atol = 1e-10,
                                                        rtol = 1e-10,
                                                        methode = "lsodes",
                                                        rhs_func = loglik_rhs)

  pars2 <- c(100, 11,0,2)
  pars <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

  loglik_DAISIE <- DAISIE:::DAISIE_loglik_CS_choice( pars1 = pars,
                                                     pars2 = pars2,
                                                     datalist = NULL,
                                                     brts = brts,
                                                     stac = 2,
                                                     missnumspec = 0,
                                                     methode = "lsodes",
                                                     CS_version = 1,
                                                     abstolint = 1E-16,
                                                     reltolint = 1E-10,
                                                     verbose = FALSE)

  testthat::expect_equal(loglik_DAISIE_trait,
                         loglik_DAISIE)


  # and now we test with 1 hidden trait (which should be the same)
  qs <- c(0.000, 0.000)
  parameters2 <- parameters
  parameters2[[5]] <- matrix(0, 2, 2)
  parameters2[[5]][1, 2] <- qs[1]
  parameters2[[5]][2, 1] <- qs[2]

  loglik_hidden <- DAISIE_DE_logpEC_trait1_hidden(brts = brts,
                                                  missnumspec = 0,
                                                  parameter = parameters2,
                                                  phy = phy,
                                                  traits = traits,
                                                  num_hidden_traits = 1)

  testthat::expect_equal(loglik_DAISIE_trait,
                         loglik_hidden)
})
