
#' This function calculates the likelihood of observing a non-endemic lineage with specified species trait states,
#'
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
#' i <- 2
#' brts <- datalist[[i]]$branching_times
#' trait <- 0
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
#' DAISIE_DE_trait_logpNE(
#'   brts                    = brts,
#'   trait                  = trait,
#'   status                  = 4,
#'   parameter               = parameter,
#'   num_observed_states     = 1,
#'   num_hidden_states       = 1,
#'   atol                    = 1e-15,
#'   rtol                    = 1e-15,
#'   methode                 = "ode45",
#'   get_initial_conditions2 = get_initial_conditions2,
#'   get_initial_conditions3 = get_initial_conditions3,
#'   get_initial_conditions4 = get_initial_conditions4,
#'   func_for_solution       = func_for_solution
#' )




DAISIE_DE_trait_logpNE <- function(brts,
                                   status,
                                   trait,
                                   trait_mainland_ancestor = FALSE,
                                   num_observed_states,
                                   num_hidden_states,
                                   parameter,
                                   atol  = 1e-15,
                                   rtol  = 1e-15,
                                   methode                 = "ode45",
                                   get_initial_conditions2 = get_initial_conditions2,
                                   get_initial_conditions3 = get_initial_conditions3,
                                   get_initial_conditions4 = get_initial_conditions4,
                                   func_for_solution       = func_for_solution) {

  # Unpack times from brts
  t0   <- brts[1]
  tmax <- brts[2]
  t1   <- brts[2]

  tp   <- 0


  # Time intervals

  time2 <<- c(tp, t1)
  time3 <<- c(tp, tmax)
  time4 <<- c(tmax, t0)

  # Number of states in the system
  num_observed_states <<- num_observed_states
  num_hidden_states <<- num_hidden_states
  n <<- num_observed_states * num_hidden_states
  parameter <<- parameter
  methode <<- methode
  atol <<- atol
  rtol <<- rtol

  # Solve for interval [tp, t2] (stem phase)


  # Run appropriate sequence of intervals
  if (status == 4) {
    interval2 <- get_func_interval(interval = "interval2")
    initial_conditions2 <<- get_initial_conditions2(status = status, res = res, trait = trait, num_observed_states = num_observed_states, num_hidden_states = num_hidden_states)
    solution2 <- func_for_solution(interval = "interval2",initial_conditions = initial_conditions2, time = time2, parameter = parameter, methode = methode, atol = atol, rtol = rtol)

    interval4 <- get_func_interval(interval = "interval4")
    initial_conditions4 <<- get_initial_conditions4(status = status, solution = solution2, parameter = parameter, trait_mainland_ancestor = trait_mainland_ancestor, num_observed_states = num_observed_states, num_hidden_states = num_hidden_states)
    solution4 <- func_for_solution(interval = "interval4", initial_conditions = initial_conditions4, time = time4, parameter = parameter, methode = methode, atol = atol, rtol = rtol)
  }

  if (status == 1) {
    interval3 <- get_func_interval(interval = "interval3")
    initial_conditions3 <<- get_initial_conditions3(status = status, res = res,  num_observed_states = num_observed_states, num_hidden_states = num_hidden_states, trait = trait)
    solution3 <- func_for_solution("interval3", initial_conditions = initial_conditions3, time = time3, parameter = parameter, methode = methode, atol = atol, rtol = rtol)

    interval4 <- get_func_interval(interval = "interval4")
    initial_conditions4 <<- get_initial_conditions4(status = status, solution = solution3, parameter = parameter, trait_mainland_ancestor = trait_mainland_ancestor,  num_observed_states = num_observed_states, num_hidden_states = num_hidden_states)
    solution4 <- func_for_solution(interval ="interval4", initial_conditions = initial_conditions4, time = time4, parameter = parameter, methode = methode, atol = atol, rtol = rtol)
  }

  # Extract log-likelihood from final solution
  Lk <- solution4[2, length(solution4[2, ])]
  logLkb <- log(Lk)
  return(logLkb)
}
