#' testing function, for comparison with DAISIE
#' @description
#' this function calculates the likelihood of observing a singleton endemic species on an island
#' with the trait state `i`, and for which only the estimated maximum and minimum ages of colonization are known.
#' @export
#' @param brts branching times
#' @param parameter parameters
#' @param num_observed_states number of observed traits
#' @param num_hidden_states number of hidden traits
#' @param trait trait state of the species at the tip
#' @param atol absolute tolerance
#' @param rtol relative tolerance
#' @param methode method of integration
#' @examples
#' library(DAISIE)
#' data("Biwa_datalist")
#' datalist <- Biwa_datalist
#' sf <- 1
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
#' DAISIE_DE_trait_logpES_max_min_age_hidden(
#'   brts                  = c(4, 3.9999, 0.001),
#'   trait                 = 0,
#'   status                = 9,
#'   parameter             = parameter,
#'   num_observed_states   = 1,
#'   num_hidden_states     = 1,
#'   atol                  = 1e-10,
#'   rtol                  = 1e-10,
#'   methode               = "ode45"
#' )
DAISIE_DE_trait_logpES_max_min_age_hidden <- function(brts,
                                                      trait,
                                                      status,
                                                      sf = 1,
                                                      parameter,
                                                      num_observed_states,
                                                      num_hidden_states,
                                                      atol = 1e-10,
                                                      rtol = 1e-10,
                                                      methode = "ode45") {
  t0   <- brts[1]
  tmax <- brts[2]
  tmin <- brts[3]
  tp   <- 0

  # number of unique state
  n <- num_observed_states * num_hidden_states

  #########interval2 [t_p, tmin]

  m = length(parameter[[1]])


  ## SOLVED: can't we call 'get_initial_conditions' here? //NO, because brts > 2
  initial_conditions2 <- get_initial_conditions2(status = status,
                                                 num_observed_states = num_observed_states,
                                                 num_hidden_states = num_hidden_states,
                                                 trait = trait,
                                                 brts = brts,
                                                 sf = sf)



  # Time sequence for interval [tp, tmin]
  time2 <- c(tp, tmin)

  # Solve the system for interval [tp, tmin]
  solution2 <- deSolve::ode(y = initial_conditions2,
                            times = time2,
                            func = interval2,
                            parms = parameter,
                            method = methode,
                            atol = atol,
                            rtol = rtol)

  solution2 <- matrix(solution2[,-1], nrow = 2) # remove the time from the result

  #########interval3 [tmin, tmax]

  # Initial conditions

  # only use second row, because the first row of solution3 is the initial state
  initial_conditions3_max_min <- c(solution2[2,][1:n],                                             ### DE: select DE in solution2
                                   rep(0, n),                                                      ### DM1: select DE in solution2
                                   solution2[2,][(n + 1):(n + n)],                         ### DM2: select DM2 in solution2
                                   solution2[2,][(n + n + 1):(n + n + n)],                 ### DM3: select DM3 in solution2
                                   solution2[2,][(n + n + n + 1):(n + n + n + n)],         ### E: select E in solution2
                                   0,                                                              ### DA2
                                   solution2[2,][length(solution2[2,])])                           ### DA3: select DA3 in solution2

  initial_conditions3_max_min <- matrix(initial_conditions3_max_min, nrow = 1)

  # Time sequence for interval [tmin, tmax]
  time3 <- c(tmin, tmax)

  # Solve the system for interval [tp, tmax]
  solution3 <- deSolve::ode(y = initial_conditions3_max_min,
                            times = time3,
                            func = interval3,
                            parms = parameter,
                            method = methode,
                            atol = atol,
                            rtol = rtol)


  solution3 <- matrix(solution3[,-1], nrow = 2) # remove the time from the result

  #########interval4 [tmax, t0]

  # Initial conditions

  # only use second row, because the first row of solution3 is the initial state
  initial_conditions4_max_min <- c(solution3[2,][(n + n + 1):(n + n + n)],                         ### DM1: select DM2 in solution3
                                   solution3[2,][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution3
                                   solution3[2,][length(solution3[2,]) - 1])

  initial_conditions4_max_min <- matrix(initial_conditions4_max_min, nrow = 1)

  # Time sequence for interval [tmax, t0]
  time4 <- c(tmax, t0)

  # Solve the system for interval [tmax, t0]
  solution4 <- deSolve::ode(y = initial_conditions4_max_min,
                            times = time4,
                            func = interval4,
                            parms = parameter,
                            method = methode,
                            atol = atol,
                            rtol = rtol)

  solution4 <- matrix(solution4[,-1], nrow = 2)

  # Extract log-likelihood
  Lk <- solution4[2,][length(solution4[2,])]
  logLkb <- log(Lk)
  return(logLkb)
}
