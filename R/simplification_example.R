#' Testing function for comparison with DAISIE
#'
#' @description
#' This function calculates the likelihood of observing a clade with specified species trait states,
#' given known colonization time. It is designed for comparison with DAISIE-based models.
#'
#' @param brts Branching times.
#' @param missnumspec Number of missing species.
#' @param parameter List of model parameters.
#' @param num_observed_states Number of observed trait states.
#' @param num_hidden_traits Number of hidden trait states.
#' @param trait_mainland_ancestor Trait state of the species at the stem (mainland ancestor).
#' @param phy Phylogeny (class 'phylo').
#' @param traits Vector of trait states for the tips.
#' @param cond Conditioning scheme (default = "proper_cond").
#' @param root_state_weight Root state weighting method (default = "proper_weights").
#' @param setting_calculation Argument used in ML optimization routines.
#' @param see_ancestral_states Logical; whether to return ancestral state reconstructions.
#' @param atol Absolute tolerance for numerical integration.
#' @param rtol Relative tolerance for numerical integration.
#' @param methode Numerical integration method (e.g., "ode45").
#'
#' @export
#'
#' @examples
#' library(DAISIE)
#' library(secsse)
#' data("Galapagos_datalist")
#' datalist <- Galapagos_datalist
#'
#' i <- 4
#' phy <- DDD::brts2phylo(datalist[[i]]$branching_times[-c(1, 2)])
#'
#' traits <- sample(c(0, 0), length(phy$tip.label), replace = TRUE)
#' sampling_fraction <- sample(c(1, 1), length(phy$tip.label), replace = TRUE)
#'
#' parameter <- list(
#'   c(2.546591, 0),                      # Cladogenesis rate + hidden
#'   c(2.678781, 0),                      # Extinction rate + hidden
#'   c(0.009326754, 0),                   # Trait-change rate + hidden
#'   c(1.008583, 0),                      # Anagenesis rate + hidden
#'   matrix(rep(0, 4), nrow = 2),         # Transition matrix Q (2x2 zeros)
#'   0                                    # Probability p
#' )
#'
#' DAISIE_DE_trait_logpEC(
#'   brts                    = datalist[[i]]$branching_times,
#'   phy                     = phy,
#'   traits                  = traits,
#'   status                  = 2,
#'   sampling_fraction       = sampling_fraction,
#'   parameter               = parameter,
#'   num_observed_states     = 2,
#'   num_hidden_traits       = 1,
#'   cond                    = "proper_cond",
#'   root_state_weight       = "proper_weights",
#'   see_ancestral_states    = TRUE,
#'   atol                    = 1e-10,
#'   rtol                    = 1e-10,
#'   methode                 = "ode45",
#'   rhs_func                = loglik_hidden_rhs,
#'   get_initial_conditions2 = get_initial_conditions2,
#'   get_initial_conditions4 = get_initial_conditions4,
#'   func_for_solution       = func_for_solution
#' )



DAISIE_DE_trait_logpEC <- function(
                                  brts,
                                  parameter,
                                  phy,
                                  traits,
                                  num_hidden_traits,
                                  num_observed_states,
                                  trait_mainland_ancestor = FALSE,
                                  status,
                                  sampling_fraction,
                                  type = "max_age_hidden", # new argument
                                  cond = "proper_cond",
                                  root_state_weight = "proper_weights",
                                  see_ancestral_states = TRUE,
                                  atol = 1e-10,
                                  rtol = 1e-10,
                                  methode = "ode45",
                                  rhs_func = loglik_hidden_rhs,
                                  get_initial_conditions2,
                                  get_initial_conditions4,
                                  func_for_solution
                              ) {
  # Unpack times from brts
  t0   <- brts[1]
  tmax <- brts[2]
  t2   <- brts[3]
  tp   <- 0

  # Number of states in the system
  num_states <- num_observed_states * num_hidden_traits

  # Solve for interval [tp, t2] (stem phase)
  res <- loglik_R_hidden(
    parameter = parameter,
    phy = phy,
    traits = traits,
    sampling_fraction = sampling_fraction,
    num_hidden_traits = num_hidden_traits,
    see_ancestral_states = see_ancestral_states,
    atol = atol,
    rtol = rtol
  )

  # Select interval functions depending on model status
  if (status == 2 || status == 3) {
    interval2 <- get_func_interval("interval2")
    interval3 <- get_func_interval("interval3")
  }

  # Initial conditions at start of crown phase
  initial_conditions2 <- get_initial_conditions2(status, res)

  # Solve for interval [t2, tmax] (crown to max age)
  time2 <- c(t2, tmax)
  solution2 <- func_for_solution(status, trait_mainland_ancestor, "interval2")

  # Solve for interval [tmax, t0] (after max age to present)
  interval4 <- get_func_interval("interval4")
  initial_conditions4 <- get_initial_conditions4(status, solution2, trait_mainland_ancestor)

  time4 <- c(tmax, t0)
  solution4 <- func_for_solution(status, trait_mainland_ancestor, "interval4")

  # Extract log-likelihood from final solution
  Lk <- solution4[2, length(solution4[2, ])]
  logLkb <- log(Lk)

  return(logLkb)
}
