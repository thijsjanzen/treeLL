test_that("identical R and Rcpp", {

  phy <- ape::read.tree(text = "((A:2.1483291,B:2.1483291):0.209405,C:2.3577341):11.7431696;")

  traits <- c(1, 1, 1)

  lambda_c <- c(1, 1)
  lambda_a <- c(0.0, 0.0)
  mus <- c(0, 0)
  gammas <- c(0.01, 0.01)
  qs <- c(0, 0)
  p <- 0

  parameters <- list()
  parameters[[1]] <- lambda_c
  parameters[[2]] <- lambda_a
  parameters[[3]] <- mus
  parameters[[4]] <- gammas
  parameters[[5]] <- qs
  parameters[[6]] <- p


  res_R <- calc_loglik(parameter = parameters,
                       phy = phy,
                       traits = traits,
                       see_ancestral_states = TRUE,
                       use_R_version = TRUE)


  res_Rcpp <- calc_loglik(parameter = parameters,
                          phy = phy,
                          traits = traits,
                          see_ancestral_states = TRUE,
                          use_R_version = FALSE)

  testthat::expect_true(all.equal(res_R, res_Rcpp))

  parameters[[5]] <- matrix(0, 2, 2)
  parameters[[5]][1, 2] <- qs[1]
  parameters[[5]][2, 1] <- qs[2]

  res_hidden <- loglik_R_hidden(parameter = parameters,
                                phy = phy,
                                traits = traits,
                                num_hidden_traits = 1,
                                see_ancestral_states = TRUE)

  testthat::expect_true(all.equal(res_R, res_hidden))
})
