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
#' @param num_hidden_states Number of hidden trait states.
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
#' i <- 5
#' phy <- DDD::brts2phylo(datalist[[i]]$branching_times[-c(1, 2)])
#' brts <- datalist[[i]]$branching_times
#' traits <- sample(c(0, 0), length(brts), replace = TRUE)
#' sampling_fraction <- sample(c(1, 1), length(brts), replace = TRUE)
#'
#'
#' parameter <- list(
#'   c(2.546591, 1.2, 1, 0.2),
#'   c(2.678781, 2, 1.9, 3),
#'   c(0.009326754, 0.003, 0.002, 0.2),
#'   c(1.008583, 1, 2, 1.5),
#'   matrix(c(
#'     0,    1,    0.5,  0,
#'     0,    0,    0.002,0.005,
#'     rep(0, 8)
#'   ), nrow = 4),
#'   0
#' )
#'
#' parameter <- list(2.546591, 2.678781, 0.009326754, 1.008583, matrix(c(0), nrow = 1), 0 )
#'
#' DAISIE_DE_trait_logpEC(
#'   brts                    = brts,
#'   phy                     = phy,
#'   traits                  = traits,
#'   status                  = 2,
#'   sampling_fraction       = sampling_fraction,
#'   parameter               = parameter,
#'   num_observed_states     = 2,
#'   num_hidden_states       = 2,
#'   atol                    = 1e-10,
#'   rtol                    = 1e-10,
#'   methode                 = "ode45")

DAISIE_DE_trait_logpEC <- function(
    brts,
    parameter,
    phy,
    traits,
    num_observed_states,
    num_hidden_states,
    trait_mainland_ancestor = FALSE,
    status,
    sampling_fraction,
    atol = 1e-10,
    rtol = 1e-10,
    methode = "ode45"
) {

  check_arguments(brts, parameter, phy, traits, num_observed_states,
                  num_hidden_states, status, sampling_fraction)

  if (length(brts) < 3) {
    stop("need at least three branching times")
  }

  # Unpack times from brts
  t0   <- brts[1]
  tmax <- brts[2]
  t1   <- brts[2]
  t2   <- brts[3]
  tp   <- 0

  # Time intervals

  time2 <- c(t2, t1)
  time3 <- c(t2, tmax)
  time4 <- c(tmax, t0)

  # Number of states in the system
  #n <- num_observed_states * num_hidden_states

  # Solve for interval [tp, t2] (stem phase)
  res <- loglik_R_tree(
    parameter = parameter,
    phy = phy,
    traits = traits,
    sampling_fraction = sampling_fraction,
    num_hidden_states = num_hidden_states,
    atol = atol,
    rtol = rtol
  )

  # Run appropriate sequence of intervals
  if ((status == 2 || status == 3) && length(brts) > 2) {

    initial_conditions2 <- get_initial_conditions2(status = status,
                                                   res = res,
                                                   trait = traits,
                                                   num_observed_states = num_observed_states,
                                                   num_hidden_states = num_hidden_states,
                                                   brts = brts,
                                                   sf = sampling_fraction,
                                                   trait_mainland_ancestor = trait_mainland_ancestor)

    solution2 <- solve_branch(interval_func = interval2,
                              initial_conditions = initial_conditions2,
                              time = time2,
                              parameter = parameter,
                              methode = methode,
                              atol = atol,
                              rtol =  rtol)

    initial_conditions4 <- get_initial_conditions4(status = status,
                                                   solution = solution2,
                                                   parameter = parameter,
                                                   trait_mainland_ancestor = trait_mainland_ancestor,
                                                   num_observed_states = num_observed_states,
                                                   num_hidden_states = num_hidden_states)

    solution4 <- solve_branch(interval = interval4,
                              initial_conditions = initial_conditions4,
                              time = time4,
                              parameter = parameter,
                              methode = methode,
                              atol = atol,
                              rtol = rtol)
  }

  if (status == 6) {
    initial_conditions3 <- get_initial_conditions3(status = status,
                                                   res = res,
                                                   num_observed_states = num_observed_states,
                                                   num_hidden_states = num_hidden_states,
                                                   trait = traits)
    solution3 <- solve_branch(interval_func = interval3,
                              initial_conditions = initial_conditions3,
                              time = time3,
                              parameter = parameter,
                              methode = methode,
                              atol = atol,
                              rtol = rtol)


    initial_conditions4 <- get_initial_conditions4(status = status,
                                                   solution = solution3,
                                                   parameter = parameter,
                                                   trait_mainland_ancestor = trait_mainland_ancestor,
                                                   num_observed_states = num_observed_states,
                                                   num_hidden_states = num_hidden_states)
    solution4 <- solve_branch(interval_func = interval4,
                              initial_conditions = initial_conditions4,
                              time = time4,
                              parameter = parameter,
                              methode = methode,
                              atol = atol,
                              rtol = rtol)
  }

  # Extract log-likelihood from final solution
  Lk <- solution4[2, length(solution4[2, ])]
  logLkb <- log(Lk)
  return(logLkb)
}
