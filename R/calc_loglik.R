calc_init_state <- function(trait) {
  out <- rep(7, 0)
  if (trait == 0) {
    out <- c(1, 0, 0, 0, 0, 0, 1)
  }
  if (trait == 1) {
    out <- c(0, 1, 0, 0, 0, 0, 1)
  }
  return(out)
}


#' @keywords internal
master_loglik <- function(parameter,
                          phy,
                          traits,
                          cond = "proper_cond",
                          root_state_weight = "proper_weights",
                          see_ancestral_states = FALSE,
                          num_threads = 1,
                          atol = 1e-8,
                          rtol = 1e-7,
                          method = "odeint::bulirsch_stoer",
                          display_warning = TRUE,
                          use_normalization = TRUE) {

  lambda_c <- parameter[[1]]
  lambda_a <- parameter[[2]]
  mus <- parameter[[3]]
  gammas <- parameter[[4]]
  qs <- parameter[[5]]
  p_value <- parameter[[6]]

  number_of_lineages <- length(phy$tip.label)

  states <- matrix(nrow = number_of_lineages + phy$Nnode,
                   ncol = 7,
                   data = NA)

  for (i in 1:length(traits)) {
    states[i, ] <- calc_init_state(traits[i])
  }

  phy$node.label <- NULL
  split_times <- sort(event_times(phy), decreasing = FALSE)
  ances <- as.numeric(names(split_times))
  forTime <- cbind(phy$edge, phy$edge.length)

  d <- ncol(states) / 2

  RcppParallel::setThreadOptions(numThreads = num_threads)

  calcul <- calc_ll_cpp(ances = ances,
                        states = states,
                        forTime = forTime,
                        lambda_cs = lambda_c,
                        lambda_as = lambda_a,
                        mus = mus,
                        gammas = gammas,
                        qs = qs,
                        p = p_value,
                        method = method,
                        atol = atol,
                        rtol = rtol,
                        see_states = see_ancestral_states,
                        use_normalization = use_normalization)

  prob_states <- calcul$merge_branch
  prob_states <- matrix(prob_states, nrow = 1,
                        dimnames = list(NULL, c("DE_0", "DE_1", "DM3_0", "DM3_1", "E_0", "E_1", "DA3")))
  return(prob_states)
}

#' @title Likelihood for any model, using Rcpp
#'
#' @inheritParams default_params_doc
#'
#' @return The loglikelihood of the data given the parameters
#' @export
calc_loglik <- function(parameter,
                        phy,
                        traits,
                        cond = "proper_cond",
                        root_state_weight = "proper_weights",
                        see_ancestral_states = FALSE,
                        num_threads = 1,
                        method = "odeint::bulirsch_stoer",
                        atol = 1e-8,
                        rtol = 1e-7,
                        display_warning = TRUE,
                        use_normalization = TRUE,
                        use_R_version = TRUE,
                        methode = "ode45",
                        rhs_func = loglik_rhs) {

  if (use_R_version) {
    return(simple_loglik_R(parameter = parameter,
                           phy = phy,
                           traits = traits,
                           cond = cond,
                           root_state_weight = root_state_weight,
                           see_ancestral_states = see_ancestral_states,
                           atol = atol,
                           rtol = rtol,
                           methode = methode,
                           rhs_func = rhs_func)
    )
  } else {
    return(master_loglik(parameter = parameter,
                phy = phy,
                traits = traits,
                cond = cond,
                root_state_weight = root_state_weight,
                see_ancestral_states = see_ancestral_states,
                num_threads = num_threads,
                atol = atol,
                rtol = rtol,
                method = method,
                display_warning = display_warning,
                use_normalization = use_normalization)
    )
  }
}
