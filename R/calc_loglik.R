#' @keywords internal
master_loglik <- function(parameter,
                          phy,
                          traits,
                          cond = "proper_cond",
                          root_state_weight = "proper_weights",
                          sampling_fraction,
                          setting_calculation = NULL,
                          see_ancestral_states = FALSE,
                          loglik_penalty = 0,
                          is_complete_tree = FALSE,
                          num_threads = 1,
                          atol = 1e-8,
                          rtol = 1e-7,
                          method = "odeint::bulirsch_stoer",
                          display_warning = TRUE,
                          use_normalization = TRUE) {

  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  q_matrix <- parameter[[3]]

  if (is.null(setting_calculation)) {
    check_input(traits,
                phy,
                sampling_fraction,
                root_state_weight,
                is_complete_tree)
    setting_calculation <- build_initStates_time(phy,
                                                 traits,
                                                 num_concealed_states,
                                                 sampling_fraction,
                                                 is_complete_tree,
                                                 mus,
                                                 num_modeled_traits,
                                                 traitStates = traitStates)
  }

  states <- setting_calculation$states
  forTime <- setting_calculation$forTime
  ances <- setting_calculation$ances

  d <- ncol(states) / 2

  RcppParallel::setThreadOptions(numThreads = num_threads)

  calcul <- calc_ll_cpp(ances = ances,
                        states = states,
                        forTime = forTime,
                        lambdas = lambdas,
                        mus = mus,
                        Q = q_matrix,
                        method = method,
                        atol = atol,
                        rtol = rtol,
                        see_states = see_ancestral_states,
                        use_normalization = use_normalization)

  loglik <- calcul$loglik
  nodeM <- calcul$node_M
  mergeBranch <- calcul$merge_branch

  if (length(nodeM) > 2 * d) nodeM <- nodeM[1:(2 * d)]

  ## At the root
  weight_states <- get_weight_states(root_state_weight,
                                     num_concealed_states,
                                     mergeBranch,
                                     lambdas,
                                     nodeM,
                                     d,
                                     is_cla = using_cla)

  wholeLike <- sum((mergeBranch2) * (weight_states))

  LL <- log(wholeLike) +
    loglik -
    penalty(pars = parameter, loglik_penalty = loglik_penalty)

  if (see_ancestral_states == TRUE) {
    states <- calcul$states
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips + 1):(nrow(states)), ]
    ancestral_states <-
      ancestral_states[, -1 * (1:(ncol(ancestral_states) / 2))]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states = ancestral_states, LL = LL, states = states))
  } else {
    return(LL)
  }
}

#' @title Likelihood for any model, using Rcpp
#'
#' @inheritParams default_params_doc
#'
#' @return The loglikelihood of the data given the parameters
#' @examples
#'rm(list=ls(all=TRUE))
#'library(secsse)
#'set.seed(13)
#'phylotree <- ape::rcoal(12, tip.label = 1:12)
#'traits <- sample(c(0,1,2),ape::Ntip(phylotree),replace=TRUE)
#'num_concealed_states <- 3
#'sampling_fraction <- c(1,1,1)
#'phy <- phylotree
#'# the idparlist for a ETD model (dual state inheritance model of evolution)
#'# would be set like this:
#'idparlist <- cla_id_paramPos(traits,num_concealed_states)
#'lambd_and_modeSpe <- idparlist$lambdas
#'lambd_and_modeSpe[1,] <- c(1,1,1,2,2,2,3,3,3)
#'idparlist[[1]] <- lambd_and_modeSpe
#'idparlist[[2]][] <- 0
#'masterBlock <- matrix(4,ncol=3,nrow=3,byrow=TRUE)
#'diag(masterBlock) <- NA
#'idparlist [[3]] <- q_doubletrans(traits,masterBlock,diff.conceal = FALSE)
#'# Now, internally, clasecsse sorts the lambda matrices, so they look like:
#'prepare_full_lambdas(traits,num_concealed_states,idparlist[[1]])
#'# which is a list with 9 matrices, corresponding to the 9 states
#'# (0A,1A,2A,0B,etc)
#'# if we want to calculate a single likelihood:
#'parameter <- idparlist
#'lambda_and_modeSpe <- parameter$lambdas
#'lambda_and_modeSpe[1,] <- c(0.2,0.2,0.2,0.4,0.4,0.4,0.01,0.01,0.01)
#'parameter[[1]] <- prepare_full_lambdas(traits,num_concealed_states,
#'lambda_and_modeSpe)
#'parameter[[2]] <- rep(0,9)
#'masterBlock <- matrix(0.07, ncol=3, nrow=3, byrow=TRUE)
#'diag(masterBlock) <- NA
#'parameter [[3]] <- q_doubletrans(traits,masterBlock,diff.conceal = FALSE)
#'cla_secsse_loglik(parameter, phy, traits, num_concealed_states,
#'                  cond = 'maddison_cond',
#'                  root_state_weight = 'maddison_weights', sampling_fraction,
#'                  setting_calculation = NULL,
#'                  see_ancestral_states = FALSE,
#'                  loglik_penalty = 0)
#'# LL = -42.18407
#' @export
calc_loglik <- function(parameter,
                        phy,
                        traits,
                        cond = "proper_cond",
                        root_state_weight = "proper_weights",
                        sampling_fraction,
                        setting_calculation = NULL,
                        see_ancestral_states = FALSE,
                        loglik_penalty = 0,
                        num_threads = 1,
                        method = "odeint::bulirsch_stoer",
                        atol = 1e-8,
                        rtol = 1e-7,
                        display_warning = TRUE,
                        use_normalization = TRUE) {
  master_loglik(parameter = parameter,
                phy = phy,
                traits = traits,
                cond = cond,
                root_state_weight = root_state_weight,
                sampling_fraction = sampling_fraction,
                setting_calculation = setting_calculation,
                see_ancestral_states = see_ancestral_states,
                loglik_penalty = loglik_penalty,
                num_threads = num_threads,
                atol = atol,
                rtol = rtol,
                method = method,
                display_warning = display_warning,
                use_normalization = use_normalization)
}
