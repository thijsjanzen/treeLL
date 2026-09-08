#' Maximum Likelihood Estimation for the DAISIE Trait Model
#'
#' @description
#' This function estimates the maximum likelihood (ML) for a given dataset using a trait-dependent model.
#' The model includes parameters for cladogenesis, extinction, colonization, anagenesis, and state transitions.
#' The optimization is done using the `treeLL::calc_ml` function.
#'
#'
#' @param datalist A list containing the data used for likelihood calculation.
#' @param num_observed_states The number of observed trait states.
#' @param num_hidden_states The number of hidden trait states.
#' @param idparslist A list containing the model parameters required for the trait-dependent diversification process.
#' Each element in the list corresponds to a specific set of parameters:
#'
#' \itemize{
#'   \item \code{idparslist[[1]]}: A numeric vector specifying the **cladogenesis rates** (\( \eqn{\lambda_c} \)) for each trait state. Each element in the vector corresponds to the cladogenesis rate for a different state.
#'   \item \code{idparslist[[2]]}: A numeric vector specifying the **extinction rates** (\( \eqn{\mu} \)) for each trait state. Each element in the vector corresponds to the extinction rate for a different state.
#'   \item \code{idparslist[[3]]}: A numeric vector specifying the **immigration rates** (\( \eqn{\gamma} \)) for each trait state. Each element in the vector corresponds to the immigration rate for a different state.
#'   \item \code{idparslist[[4]]}: A numeric vector specifying the **anagenesis rates** (\( \eqn{\lambda_a} \)) for each trait state. Each element in the vector corresponds to the anagenesis rate for a different state.
#'   \item \code{idparslist[[5]]}: A **transition matrix** for the model, describing the transition rates between all trait states. The matrix is square, and the entries represent the transition rates between states:
#'     \itemize{
#'       \item \code{idparslist[[5]][i, j]}: The transition rate from state \(i\) to state \(j\), where \(i\) and \(j\) range from 1 to the total number of states.
#'       \item Diagonal elements \code{idparslist[[5]][i, i]} represent the self-transition rate are are equal to 0.
#'     }
#'   \item \code{idparslist[[6]]}: A fixed value for the parameter \( p \), which specifies the probability that a trait transition results in a new species. If \( p = 1 \), every transition will result in a new species. If \( p = 0 \), the transition does not lead to the creation of a new species.
#' }

#' @param idparsopt A vector of integers specifying which parameters are to be optimized.
#' @param initparsopt A numeric vector providing the initial values for the parameters to be optimized.
#' @param idparsfix A vector of integers specifying which parameters are to be kept fixed.
#' @param parsfix A vector of values corresponding to the fixed parameters specified in `idparsfix`.
#' @param atol  A numeric specifying the absolute tolerance of integration.
#' @param rtol  A numeric specifying the relative tolerance of integration.
#' @param num_threads number of threads to be used. Default is one thread.
#' @param verbose sets verbose output; default is TRUE when optimmethod is "simplex". If optimmethod is set to "simplex", then even if set to FALSE, optimizer output will be shown.
#' @param use_Rcpp Integer. Specifies whether to use C++ for optimization:
#'   \itemize{
#'     \item 0: Use R implementation (default).
#'     \item 1: Use a mix of R and C++.
#'     \item 2: Use C++ implementation.
#'   }
#'
#' @details
#' This function runs the maximum likelihood estimation for a dataset based on the trait-dependent model.
#' The optimization process uses the `treeLL::calc_ml` function to adjust the parameters for the best fit to the data.
#' The function assumes that the dataset is prepared correctly, and the parameter IDs are provided for both optimization and fixing.
#'
#' @returns
#' A list containing the following elements:
#' \itemize{
#'   \item `MLpars`: The optimized parameter values.
#'   \item `ML`: The calculated log-likelihood value for the model.
#'   \item `conv`: The convergence status of the optimization process.
#' }
#'
#' @param cond conditioning on colonisation, passed on to the likelihood.
#' @param tol tolerances of the optimisation, as c(reltolx, reltolf, abstolx).
#' @param maxiter maximum number of iterations of the optimisation.
#' @param optimmethod optimisation method, "simplex" (default) or "subplex".
#' @param methode used integration method, default is ode45.
#' @param rcpp_methode integration method used when integrating with Rcpp.
#' @param num_cycles number of cycles of the optimisation.
#' @param sampling either "rho" (the default), in which each colonist
#' contributes the rho-sampling likelihood evaluated at its
#' \code{sampling_fraction}, or "n", in which it instead contributes the
#' n-sampling likelihood for its \code{missing_species}. Since
#' \code{missing_species} is a single total that does not say which trait
#' states the unsampled species are in, every split of that total over the
#' observed states contributes and their likelihoods are added. The conversion
#' is done by \code{DAISIE_DE_trait_n}, which differentiates the rho-sampling
#' likelihood around rho = 1.
#'
#' @examples
#' \dontrun{
#' ############## Parameter Estimation
#' ### Binary State Model without Hidden States
#'
#' data("Galapagos_datalist", package = "DAISIE")
#' datalist <- Galapagos_datalist
#'
#' # Mainland species that never colonised, split over the observed trait
#' # states, plus those whose trait state is unknown.
#' datalist[[1]]$not_present_by_state <- rep(datalist[[1]]$not_present / 2, 2)
#' datalist[[1]]$not_present_NA <- 0
#'
#' # Every colonist needs a phylogeny, the trait state of each sampled species
#' # and the trait distribution of its mainland ancestor. It also needs a
#' # sampling fraction per observed state (used when sampling = "rho") and the
#' # total number of unsampled species (used when sampling = "n").
#' for (i in 2:length(datalist)) {
#'   datalist[[i]]$phylogeny <-
#'     DDD::brts2phylo(datalist[[i]]$branching_times[-c(1, 2)])
#'   datalist[[i]]$phylogeny$root.edge <- 0
#'   num_tips <- length(datalist[[i]]$phylogeny$tip.label)
#'   datalist[[i]]$traits <- sample(c(0, 1), size = num_tips, replace = TRUE)
#'   datalist[[i]]$root_state <- c(0.5, 0.5)
#'   datalist[[i]]$sampling_fraction <- rep(1, 2)
#'   datalist[[i]]$missing_species <- 0
#' }
#'
#' # Define the parameter list for the model
#' idparslist <- list()
#' idparslist[[1]] <- c(1, 1)  # lambda_c, one per trait state
#' idparslist[[2]] <- c(2, 2)  # mu
#' idparslist[[3]] <- c(3, 3)  # gamma
#' idparslist[[4]] <- c(4, 4)  # lambda_a
#'
#' # Transition matrix for the model (since there are two trait states). The
#' # diagonal keeps id 0, which is fixed at 0 below.
#' idparslist[[5]] <- matrix(0, 2, 2)
#' colnames(idparslist[[5]]) <- c("0", "1")
#' rownames(idparslist[[5]]) <- colnames(idparslist[[5]])
#' idparslist[[5]][1, 2] <- 5  # 0 -> 1
#' idparslist[[5]][2, 1] <- 6  # 1 -> 0
#'
#' # p, the probability that a transition results in a new species. It must
#' # have the same dimensions as the transition matrix.
#' idparslist[[6]] <- matrix(7, 2, 2)
#'
#' # Parameters to optimize, and the ones that stay fixed. Every id used in
#' # idparslist has to appear in exactly one of the two, id 0 included.
#' idparsopt   <- 1:6
#' initparsopt <- c(1.07, 1.0102, 0.0035, 0.174, 0.001, 0.001)
#' idparsfix   <- c(0, 7)
#' parsfix     <- c(0, 0)
#'
#' # Run the ML estimation with the provided dataset. use_Rcpp = 0 is needed
#' # while p is a matrix, because cpp_solve still takes a scalar p.
#' ml_rho <- treeLL::calc_ml(datalist,
#'                           num_observed_states = 2,
#'                           num_hidden_states = 1,
#'                           idparslist = idparslist,
#'                           idparsopt = idparsopt,
#'                           initparsopt = initparsopt,
#'                           idparsfix = idparsfix,
#'                           parsfix = parsfix,
#'                           atol = 1e-12,
#'                           rtol = 1e-12,
#'                           num_threads = 8,
#'                           verbose = TRUE,
#'                           use_Rcpp = 0,
#'                           sampling = "rho")
#'
#' # The same fit, but with the missing species given per observed trait state
#' # instead of a sampling fraction.
#' ml_n <- treeLL::calc_ml(datalist,
#'                         num_observed_states = 2,
#'                         num_hidden_states = 1,
#'                         idparslist = idparslist,
#'                         idparsopt = idparsopt,
#'                         initparsopt = initparsopt,
#'                         idparsfix = idparsfix,
#'                         parsfix = parsfix,
#'                         atol = 1e-12,
#'                         rtol = 1e-12,
#'                         num_threads = 8,
#'                         verbose = TRUE,
#'                         use_Rcpp = 0,
#'                         sampling = "n")
#'
#' # View the results
#' print(ml_rho)
#' print(ml_n)
#' }
#'
#'
#' @export


calc_ml <- function(datalist,
                    num_observed_states,
                    num_hidden_states,
                    idparslist,
                    idparsopt,
                    initparsopt,
                    idparsfix,
                    parsfix,
                    cond = "proper_cond",
                    tol = c(1e-04, 1e-05, 1e-07),
                    maxiter = 2000 * round((1.25) ^ length(idparsopt)),
                    optimmethod = "simplex",
                    methode = "ode45",
                    rcpp_methode = "odeint::runge_kutta_cash_karp54",
                    num_cycles = 1,
                    verbose = FALSE,
                    num_threads = 1,
                    atol = 1e-15,
                    rtol = 1e-15,
                    use_Rcpp = 0,
                    sampling = "rho"
) {
  if (identical(as.numeric(sort(c(idparsopt, idparsfix))),
                as.numeric(sort(unique(unlist(idparslist))))) == FALSE) {
    stop("All elements in idparslist must be included in either
             idparsopt or idparsfix ")
  }

  see_ancestral_states <- FALSE

  trparsopt <- initparsopt / (1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1

  optimpars <- c(tol, maxiter, verbose)

  initloglik <- loglik_choosepar(trparsopt = trparsopt,
                                 trparsfix = trparsfix,
                                 idparsopt = idparsopt,
                                 idparsfix = idparsfix,
                                 idparslist = idparslist,
                                 datalist = datalist,
                                 num_observed_states = num_observed_states,
                                 num_hidden_states = num_hidden_states,
                                 cond = cond,
                                 atol = atol,
                                 rtol = rtol,
                                 methode = methode,
                                 rcpp_methode = rcpp_methode,
                                 verbose = verbose,
                                 use_Rcpp = use_Rcpp,
                                 num_threads = num_threads,
                                 sampling = sampling)
  # Function here
  if (verbose) print_init_ll(initloglik = initloglik)

  if (initloglik == -Inf) {
    stop("The initial parameter values have a likelihood that is
             equal to 0 or below machine precision.
             Try again with different initial values.")
  } else {
    out <- DDD::optimizer(optimmethod = optimmethod,
                          optimpars = optimpars,
                          fun = loglik_choosepar,
                          trparsopt = trparsopt,
                          num_cycles = num_cycles,
                          idparsopt = idparsopt,
                          trparsfix = trparsfix,
                          idparsfix = idparsfix,
                          idparslist = idparslist,
                          datalist = datalist,
                          num_observed_states = num_observed_states,
                          num_hidden_states = num_hidden_states,
                          cond = cond,
                          atol = atol,
                          rtol = rtol,
                          methode = methode,
                          rcpp_methode = rcpp_methode,
                          verbose = verbose,
                          use_Rcpp = use_Rcpp,
                          num_threads = num_threads,
                          sampling = sampling)
    if (out$conv != 0) {
      stop("Optimization has not converged.
                 Try again with different initial values.")
    } else {
      ml_pars1 <- transform_parameters(as.numeric(unlist(out$par)),
                                       trparsfix,
                                       idparsopt,
                                       idparsfix,
                                       idparslist)
      names(ml_pars1) <- c(
        "ML_lambda_c",
        "ML_mu",
        "ML_gamma",
        "ML_lambda_a",
        "ML_q",
        "p"
      )
      out2 <- list(MLpars = ml_pars1,
                   ML = as.numeric(unlist(out$fvalues)),
                   conv = out$conv)
    }
  }
  return(out2)
}

#' loglik choosepar temp
#' @description
#' temporary export for testing
#' @inheritParams default_params_doc
#' @export
loglik_choosepar <- function(trparsopt,
                             trparsfix,
                             idparsopt,
                             idparsfix,
                             idparslist,
                             datalist,
                             num_observed_states,
                             num_hidden_states,
                             cond = cond,
                             atol,
                             rtol,
                             methode,
                             rcpp_methode,
                             verbose,
                             use_Rcpp,
                             sampling = "rho",
                             num_threads) {
  alltrpars <- c(trparsopt, trparsfix)

  loglik <- NA

  if (max(alltrpars) > 1 || min(alltrpars) < 0) {
    loglik <- -Inf
  } else {
    pars1 <- transform_parameters(trparsopt, trparsfix,
                                  idparsopt, idparsfix,
                                  idparslist)

    loglik <- DAISIE_DE_trait_loglik_CS(parameter = pars1,
                                        datalist = datalist,
                                        methode = methode,
                                        rcpp_methode = rcpp_methode,
                                        atol = atol,
                                        rtol = rtol,
                                        num_observed_states =
                                          num_observed_states,
                                        num_hidden_states = num_hidden_states,
                                        cond = cond,
                                        verbose = verbose,
                                        use_Rcpp = use_Rcpp,
                                        num_threads = num_threads,
                                        sampling = sampling)

    if (is.nan(loglik) || is.na(loglik)) {
      warning("There are parameter values used which cause
                numerical problems.")
      loglik <- -Inf
    }
  }
  return(loglik)
}

