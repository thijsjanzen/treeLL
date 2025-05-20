#' testing fuction, for comparison with DAISIE
#' @description
#' this function calculates the likelihood of observing a singleton species on an island
#' with the trait state `i`, either non-endemic or rendered endemic by a trait change, and
#' for which only the estimated maximum age of colonization is known.
#' @export
#' @param brts branching times
#' @param parameter parameters
#' @param num_observed_states number of observed traits
#' @param num_hidden_states number of hidden traits
#' @param trait trait
#' @param cond conditioning, default = "proper_cond"
#' @param root_state_weight root weight, default = "proper_weights"
#' @param setting_calculation used in ML
#' @param see_ancestral_states recover the ancestral states
#' @param atol absolute tolerance
#' @param rtol relative tolerance
#' @param methode method of integration
#' @examples
#' # load DAISIE package and data
#' library(DAISIE)
#' parameter <- list(
#'   c(2.546591, 1.2, 1, 0.2),
#'   c(2.678781, 2, 1.9, 3),
#'   c(0.009326754, 0.003, 0.002, 0.2),
#'   c(1.008583, 1, 2, 1.5),
#'   matrix(c(0, 1, 0.5, 0,
#'            0, 0, 0.002, 0.005,
#'            rep(0, 8)),
#'          nrow = 4),
#'   0
#' )
#'
#' DAISIE_DE_trait_logpNE_max_min_age_hidden(
#'   brts                  = c(4, 3.999, 0.0001),
#'   trait                 = 0,
#'   parameter             = parameter,
#'   num_observed_states   = 2,
#'   num_hidden_states     = 2,
#'   cond                  = "proper_cond",
#'   root_state_weight     = "proper_weights",
#'   see_ancestral_states  = TRUE,
#'   atol                  = 1e-15,
#'   rtol                  = 1e-15,
#'   methode               = "ode45"
#' )

DAISIE_DE_trait_logpNE_max_min_age_hidden <- function(brts,
                                                  trait,
                                                  parameter,
                                                  num_observed_states,
                                                  num_hidden_states,
                                                  cond = "proper_cond",
                                                  root_state_weight = "proper_weights",
                                                  see_ancestral_states = TRUE,
                                                  atol = 1e-10,
                                                  rtol = 1e-10,
                                                  methode = "ode45") {
  t0   <- brts[1]
  tmax <- brts[2]
  tmin <- brts[3]
  tp   <- 0

  # number of unique state
  n <- num_observed_states * num_hidden_states

  #########Interval1 [t_p, tmin]

  interval1 <- function(t, state, parameter) {
    with(as.list(c(state, parameter)), {

      lambdac <- parameter[[1]]
      mu      <- parameter[[2]]
      gamma   <- parameter[[3]]
      lambdaa <- parameter[[4]]
      q       <- parameter[[5]]
      p       <- parameter[[6]]


      # n <- (length(state) - 1) / 4
      n <- num_observed_states * num_hidden_states

      dDE     <- numeric(n)
      dDM2    <- numeric(n)
      dDM3    <- numeric(n)
      dE      <- numeric(n)

      t_vec   <- rowSums(q)

      DE  <- state[1:n]
      DM2 <- state[(n + 1):(n + n)]
      DM3 <- state[(n + n + 1):(n + n + n)]
      E   <- state[(n + n + n + 1):(n + n + n + n)]
      DA3 <- state[length(state)]

      gamma_matrix <- matrix(gamma, nrow = n, ncol = length(gamma), byrow = TRUE)
      gamma_nonself <- rowSums(gamma_matrix - diag(gamma))

      q_mult_E   <- t(q %*% E)
      q_mult_DE  <- t(q %*% DE)
      q_mult_DM2 <- t(q %*% DM2)
      q_mult_DM3 <- t(q %*% DM3)


      dDE <- -(lambdac + mu + t_vec) * DE +
        2 * lambdac * DE * E +
        q_mult_DE


      dDM2 <- -(lambdac + mu + gamma + lambdaa + t_vec) * DM2 +
        (lambdaa * DE + 2 * lambdac * DE * E + p * q_mult_DE) * DA3 +
        (1 - p) * q_mult_DM2 + gamma_nonself * DM2



      dDM3 <-  -(lambdac + mu + gamma_nonself + lambdaa + t_vec) * DM3 +
        (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA3 +
        (1 - p) * q_mult_DM3 +
        gamma_nonself * DM3

      dE <- mu - (mu + lambdac + t_vec) * E +
        lambdac * E * E +
        q_mult_E

      dDA3 <- -sum(gamma) * DA3 + sum(gamma * DM3)

      return(list(c(dDE, dDM2, dDM3, dE, dDA3)))
    })
  }


  m = length(parameter[[1]])


  calc_init_state_hidden <- function(trait,
                                     num_unique_states,
                                     num_hidden_states) {

    DE  <- rep(0, num_unique_states)
    DM2 <- rep(0, num_unique_states)
    DM3 <- rep(0, num_unique_states)
    E   <- rep(0, num_unique_states)
    DA3 <- 1

    #for (i in 1:num_hidden_states) {
    # assuming the traits start counting at 0 !!!!
    # DM2[(1 + trait) + (i - 1) * num_hidden_states] <- 1
    #}
    DM2[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- 1

    return( c(DE, DM2, DM3, E, DA3))
  }





  num_unique_states <- length(parameter[[1]])
  initial_conditions1 <-   calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  initial_conditions1 <- matrix(initial_conditions1, nrow = 1)




  # Time sequence for interval [tp, tmin]
  time1 <- c(tp, tmin)

  # Solve the system for interval [tp, tmin]
  solution1 <- deSolve::ode(y = initial_conditions1,
                            times = time1,
                            func = interval1,
                            parms = parameter,
                            method = methode,
                            atol = atol,
                            rtol = rtol)


  solution1 <- matrix(solution1[,-1], nrow = 2) # remove the time from the result

  #########Interval2 [tmin, tmax]


  interval2 <- function(t, state, parameter) {
    with(as.list(c(state, parameter)), {

      lambdac <- parameter[[1]]
      mu      <- parameter[[2]]
      gamma   <- parameter[[3]]
      lambdaa <- parameter[[4]]
      q       <- parameter[[5]]
      p       <- parameter[[6]]


      # n <- (length(state) - 1) / 4
      n <- num_observed_states*num_hidden_states
      dDE     <- numeric(n)
      dDM1    <- numeric(n)
      dDM2    <- numeric(n)
      dDM3    <- numeric(n)
      dE      <- numeric(n)

      t_vec   <- rowSums(q)

      DE  <- state[1:n]
      DM1 <- state[(n + 1):(n + n)]
      DM2 <- state[(n + n + 1):(n + n + n)]
      DM3 <- state[(n + n + n + 1):(n + n + n + n)]
      E   <- state[(n + n + n + n + 1):(n + n + n + n + n)]
      DA2 <- state[length(state) - 1]
      DA3 <- state[length(state)]

      gamma_matrix <- matrix(gamma, nrow = n, ncol = length(gamma), byrow = TRUE)
      gamma_nonself <- rowSums(gamma_matrix - diag(gamma))

      q_mult_E   <- t(q %*% E)
      q_mult_DE  <- t(q %*% DE)
      q_mult_DM1 <- t(q %*% DM1)
      q_mult_DM2 <- t(q %*% DM2)
      q_mult_DM3 <- t(q %*% DM3)


      dDE <- -(lambdac + mu + t_vec) * DE +
        2 * lambdac * DE * E +
        q_mult_DE

      dDM1 <- -(lambdac + mu + gamma + lambdaa + t_vec) * DM1 +
        (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA2 +
        (1 - p) * q_mult_DM1 + gamma * DM2


      dDM2 <- -(lambdac + mu + gamma_nonself + lambdaa + t_vec) * DM2 +
        (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA2 +
        (lambdaa * DE + 2*lambdac * DE + p * q_mult_DE) * DA3 +
        (1 - p) * q_mult_DM2 + gamma_nonself * DM2



      dDM3 <-  -(lambdac + mu + gamma_nonself + lambdaa + t_vec) * DM3 +
        (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA3 +
        (1 - p) * q_mult_DM3 +
        gamma_nonself * DM3

      dE <- mu - (mu + lambdac + t_vec) * E +
        lambdac * E * E +
        q_mult_E

      dDA2 <- -sum(gamma) * DA2 + sum(gamma * DM2)

      dDA3 <- -sum(gamma) * DA3 + sum(gamma * DM3)

      return(list(c(dDE, dDM1, dDM2, dDM3, dE, dDA2, dDA3)))
    })
  }

  # Initial conditions

  # only use second row, because the first row of solution1 is the initial state
  initial_conditions2 <- c(solution1[2,][1:n],
                           rep(0, n),       ### DE: select DE in solution1
                           solution1[2,][(n + 1):(n + n)],         ### DM2: select DM2 in solution1
                           solution1[2,][(n + n + 1):(n + n + n)],         ### DM3: select DM3 in solution1
                           solution1[2,][(n + n + n + 1):(n + n + n + n)],         ### E: select E in solution1
                           0,
                           solution1[2,][length(solution1[2,])])                       ### DA3: select DA3 in solution1

  initial_conditions2 <- matrix(initial_conditions2, nrow = 1)

  # Time sequence for interval [tmin, tmax]
  time1 <- c(tmin, tmax)

  # Solve the system for interval [tp, tmax]
  solution2 <- deSolve::ode(y = initial_conditions2,
                            times = time2,
                            func = interval2,
                            parms = parameter,
                            method = methode,
                            atol = atol,
                            rtol = rtol)


  solution2 <- matrix(solution2[,-1], nrow = 2) # remove the time from the result

  #########Interval3 [tmax, t0]

  interval3 <- function(t, state, parameter) {
    with(as.list(c(state, parameter)), {

      lambdac <- parameter[[1]]
      mu      <- parameter[[2]]
      gamma   <- parameter[[3]]
      lambdaa <- parameter[[4]]
      q       <- parameter[[5]]
      p       <- parameter[[6]]


      #n <- (length(state) - 1) / 2
      n <- num_observed_states * num_hidden_states

      dDM1    <- numeric(n)
      dDE     <- numeric(n)


      t_vec   <- rowSums(q)

      DM1  <- state[1:n]
      E    <- state[(n + 1):(n + n)]
      DA1  <- state[length(state)]

      gamma_matrix <- matrix(gamma, nrow = n, ncol = length(gamma), byrow = TRUE)
      gamma_nonself <- rowSums(gamma_matrix - diag(gamma))

      q_mult_E   <- t(q %*% E)
      q_mult_DM1 <- t(q %*% DM1)

      dDM1 <-  -(lambdac + mu + gamma_nonself + lambdaa + t_vec) * DM1 +
        (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA1 +
        (1 - p) * q_mult_DM1 + gamma_nonself * DM1

      dE <- mu - (mu + lambdac + t_vec) * E +
        lambdac * E * E +
        q_mult_E

      dDA1 <- -sum(gamma) * DA1 + sum(gamma * DM1)

      return(list(c(dDM1, dE, dDA1)))
    })
  }


  # Initial conditions

  # only use second row, because the first row of solution2 is the initial state
  initial_conditions3 <- c(solution2[2,][(n + 1):(n + n)],                                 ### DM1: select DM2 in solution2
                           solution2[2,][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution2
                           solution2[2,][length(solution2[2,]) - 1])                       ### DA1: select DA2 in solution2

  initial_conditions3 <- matrix(initial_conditions3, nrow = 1)

  # Time sequence for interval [tmax, t0]
  time3 <- c(tmax, t0)

  # Solve the system for interval [tmax, t0]
  solution3 <- deSolve::ode(y = initial_conditions3,
                            times = time3,
                            func = interval3,
                            parms = parameter,
                            method = methode,
                            atol = atol,
                            rtol = rtol)

  solution3 <- matrix(solution3[,-1], nrow = 2)

  # Extract log-likelihood
  Lk <- solution3[2,][length(solution3[2,])]
  logLkb <- log(Lk)
  return(logLkb)
}
