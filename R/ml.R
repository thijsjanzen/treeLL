#' @keywords internal
master_ml <- function(phy,
                      traits,
                      idparslist,
                      idparsopt,
                      initparsopt,
                      idparsfix,
                      parsfix,
                      cond = "proper_cond",
                      root_state_weight = "proper_weights",
                      tol = c(1e-04, 1e-05, 1e-07),
                      maxiter = 1000 * round((1.25) ^ length(idparsopt)),
                      optimmethod = "simplex",
                      num_cycles = 1,
                      verbose = FALSE,
                      num_threads = 1,
                      atol = 1e-8,
                      rtol = 1e-7,
                      method = "odeint::bulirsch_stoer",
                      use_normalization = TRUE) {


  if (identical(as.numeric(sort(c(idparsopt, idparsfix))),
                as.numeric(sort(unique(unlist(idparslist))))) == FALSE) {
    stop("All elements in idparslist must be included in either
             idparsopt or idparsfix ")
  }


  check_ml_conditions(traits,
                      idparslist,
                      initparsopt,
                      idparsopt,
                      idparsfix,
                      parsfix)

  see_ancestral_states <- FALSE

  trparsopt <- initparsopt / (1 + initparsopt)
  trparsopt[which(initparsopt == Inf)] <- 1
  trparsfix <- parsfix / (1 + parsfix)
  trparsfix[which(parsfix == Inf)] <- 1

  optimpars <- c(tol, maxiter, verbose)

  ll_verbose <- ifelse(optimmethod == "subplex",
                       verbose,
                       FALSE)
  initloglik <- loglik_choosepar(trparsopt = trparsopt,
                                        trparsfix = trparsfix,
                                        idparsopt = idparsopt,
                                        idparsfix = idparsfix,
                                        idparslist = idparslist,
                                        phy = phy,
                                        traits = traits,
                                        cond = cond,
                                        root_state_weight = root_state_weight,
                                         see_ancestral_states =
                                          see_ancestral_states,
                                        num_threads = num_threads,
                                        atol = atol,
                                        rtol = rtol,
                                        method = method,
                                        display_warning = FALSE,
                                        verbose = ll_verbose,
                                        use_normalization = use_normalization)
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
                          phy = phy,
                          traits = traits,
                          cond = cond,
                          root_state_weight = root_state_weight,
                          see_ancestral_states = see_ancestral_states,
                          num_threads = num_threads,
                          atol = atol,
                          rtol = rtol,
                          method = method,
                          display_warning = FALSE,
                          verbose = ll_verbose,
                          use_normalization = use_normalization)
    if (out$conv != 0) {
      stop("Optimization has not converged.
                 Try again with different initial values.")
    } else {
      ml_pars1 <- secsse_transform_parameters(as.numeric(unlist(out$par)),
                                              trparsfix,
                                              idparsopt,
                                              idparsfix,
                                              idparslist)
      out2 <- list(MLpars = ml_pars1,
                   ML = as.numeric(unlist(out$fvalues)),
                   conv = out$conv)
    }
  }
  return(out2)
}

#' @keywords internal
loglik_choosepar <- function(trparsopt,
                             trparsfix,
                             idparsopt,
                             idparsfix,
                             idparslist,
                             phy,
                             traits,
                             cond = cond,
                             root_state_weight,
                             see_ancestral_states,
                             num_threads,
                             atol,
                             rtol,
                             method,
                             display_warning,
                             verbose,
                             use_normalization) {
  alltrpars <- c(trparsopt, trparsfix)
  if (max(alltrpars) > 1 || min(alltrpars) < 0) {
    loglik <- -Inf
  } else {
    pars1 <- secsse_transform_parameters(trparsopt, trparsfix,
                                         idparsopt, idparsfix,
                                         idparslist)

    loglik <- master_loglik(parameter = pars1,
                            phy = phy,
                            traits = traits,
                            cond = cond,
                            root_state_weight =
                              root_state_weight,
                             see_ancestral_states =
                              see_ancestral_states,
                            num_threads = num_threads,
                            method = method,
                            atol = atol,
                            rtol = rtol,
                            display_warning = display_warning,
                            use_normalization = use_normalization)

    if (is.nan(loglik) || is.na(loglik)) {
      warning("There are parameter values used which cause
                numerical problems.")
      loglik <- -Inf
    }
  }
  if (verbose) {
    out_print <- c(trparsopt / (1 - trparsopt), loglik)
    message(paste(out_print, collapse = " "))
  }
  return(loglik)
}

#' Maximum likehood estimation for (SecSSE)
#'
#' Maximum likehood estimation under Several examined and concealed
#' States-dependent Speciation and Extinction (SecSSE) with cladogenetic option
#'
#' @inheritParams default_params_doc
#'
#' @return Parameter estimated and maximum likelihood
#' @export
calc_ml <- function(phy,
                    traits,
                    idparslist,
                    idparsopt,
                    initparsopt,
                    idparsfix,
                    parsfix,
                    cond = "proper_cond",
                    root_state_weight = "proper_weights",
                     tol = c(1e-04, 1e-05, 1e-07),
                    maxiter = 1000 * round((1.25)^length(idparsopt)),
                    optimmethod = "simplex",
                    num_cycles = 1,
                    verbose = FALSE,
                    num_threads = 1,
                    atol = 1e-8,
                    rtol = 1e-7,
                    method = "odeint::bulirsch_stoer",
                    use_normalization = TRUE) {
  master_ml(phy = phy,
            traits = traits,
            idparslist = idparslist,
            idparsopt = idparsopt,
            initparsopt = initparsopt,
            idparsfix = idparsfix,
            parsfix = parsfix,
            cond = cond,
            root_state_weight = root_state_weight,
            tol = tol,
            maxiter = maxiter,
            optimmethod = optimmethod,
            num_cycles = num_cycles,
            verbose = verbose,
            num_threads = num_threads,
            atol = atol,
            rtol = rtol,
            method = method,
            use_normalization = use_normalization)
}
