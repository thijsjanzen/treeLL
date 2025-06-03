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
#' @param trait trait state of the species at the tip.
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
#' i <- 7
#' brts <- datalist[[i]]$branching_times
#' trait <- 0
#' sf <- 1
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
#'
#' DAISIE_DE_trait_logpES(
#'   brts                    = brts,
#'   trait                  = trait,
#'   status                  = 2,
#'   parameter               = parameter,
#'   num_observed_states     = 2,
#'   num_hidden_states       = 2,
#'   atol                    = 1e-10,
#'   rtol                    = 1e-10,
#'   methode                 = "ode45"
#' )
DAISIE_DE_trait_logpES <- function(brts,
                                   status,
                                   trait,
                                   sf = 1,
                                   trait_mainland_ancestor = FALSE,
                                   num_observed_states,
                                   num_hidden_states,
                                   parameter,
                                   atol  = 1e-10,
                                   rtol  = 1e-10,
                                   methode                 = "ode45",
                                   use_R = TRUE) {

  check_arguments(brts, parameter,
                  phy = 0,
                  trait,
                  num_observed_states,
                  num_hidden_states,
                  status,
                  sampling_fraction = 0)

  # Unpack times from brts
  t0   <- brts[1]
  tmax <- brts[2]
  t1   <- brts[2]
  tp   <- 0

  # Time intervals

  time2 <- c(tp, t1)
  time3 <- c(tp, tmax)
  time4 <- c(tmax, t0)

  # Number of states in the system
  #n <- num_observed_states * num_hidden_states

  # Solve for interval [tp, t2] (stem phase)


  # Run appropriate sequence of intervals
  if ((status == 2 || status == 3) && length(brts) == 2) {
    initial_conditions2 <- get_initial_conditions2(status = status,
                                                   num_observed_states = num_observed_states,
                                                   num_hidden_states = num_hidden_states,
                                                   trait = trait,
                                                   brts = brts,
                                                   sf = sf,
                                                   trait_mainland_ancestor = trait_mainland_ancestor)

    solution2 <- solve_branch(interval_func = interval2,
                              initial_conditions = initial_conditions2,
                              time = time2,
                              parameter = parameter,
                              methode = methode,
                              atol = atol,
                              rtol = rtol,
                              use_R = use_R)


    initial_conditions4 <- get_initial_conditions4(status = status,
                                                   solution = solution2,
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
                              rtol = rtol,
                              use_R = use_R)
  }

  if (status == 5) {
    initial_conditions3 <- get_initial_conditions3(status = status,
                                                   num_observed_states = num_observed_states,
                                                   num_hidden_states = num_hidden_states,
                                                   trait = trait)
    solution3 <- solve_branch(interval_func = interval3,
                              initial_conditions = initial_conditions3,
                              time = time3,
                              parameter = parameter,
                              methode = methode,
                              atol = atol,
                              rtol = rtol,
                              use_R = use_R)

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
                              rtol = rtol,
                              use_R = use_R)
  }

  # Extract log-likelihood from final solution
  Lk <- solution4[2, length(solution4[2, ])]
  logLkb <- log(Lk)
  return(logLkb)
}
