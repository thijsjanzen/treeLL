#' testing fuction, for comparison with DAISIE
#' @description
#' this function calculates the likelihood of observing a singleton endemic species on an island
#' with the trait state `i`, and for which the estimated colonization time is known.
#' @export
#' @param brts branching times
#' @param missnumspec number of missing species
#' @param parameter parameters
#' @param num_observed_states number of observed traits
#' @param num_hidden_states number of hidden traits
#' @param trait trait state of the species at the tip
#' @param atol absolute tolerance
#' @param rtol relative tolerance
#' @param methode method of integration
#' @param sf sampling fraction
#' @examples
#' library(DAISIE)
#' data("Biwa_datalist")
#' datalist <- Biwa_datalist
#' i <- 49
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
#' DAISIE_DE_trait_logpES_hidden(
#'   brts                   = datalist[[i]]$branching_times,
#'   trait                  = 0,
#'   trait_mainland_ancestor= 0,
#'   parameter              = parameter,
#'   num_observed_states    = 2,
#'   num_hidden_states      = 2,
#'   atol                   = 1e-10,
#'   rtol                   = 1e-10,
#'   methode                = "ode45"
#' )


DAISIE_DE_trait_logpES_hidden <- function(brts,
                                          trait,
                                          sf = 1,
                                          trait_mainland_ancestor = "FALSE",
                                          parameter,
                                          num_observed_states,
                                          num_hidden_states,
                                          atol = 1e-10,
                                          rtol = 1e-10,
                                          methode = "ode45") {
  t0 <- brts[1]
  t1 <- brts[2]
  tp <- 0
  # number of unique state
  n <- num_observed_states * num_hidden_states

  #########interval2 [t_p, t_1]

  interval2 <- function(t, state, parameter) {
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

      gamma_matrix <- matrix(parameter[[3]], nrow = n, ncol = length(parameter[[3]]), byrow = TRUE)

      if (nrow(gamma_matrix) == 1) {
        # when there's only one row, there is no “self” element to subtract
        gamma_nonself <- 0
      } else {
        # for n > 1, subtract the diagonal (self‐effects) as before
        gamma_nonself <- rowSums(gamma_matrix - diag(parameter[[3]]))
      }

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



    DE[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- sf
    E[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- 1 - sf

    return( c(DE, DM2, DM3, E, DA3))
  }





  num_unique_states <- length(parameter[[1]])
  initial_conditions2 <-   calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  initial_conditions2 <- matrix(initial_conditions2, nrow = 1)




  # Time sequence for interval [tp, t1]
  time2 <- c(tp, t1)

  # Solve the system for interval [tp, t1]
  solution2 <- deSolve::ode(y = initial_conditions2,
                            times = time2,
                            func = interval2,
                            parms = parameter,
                            method = methode,
                            atol = atol,
                            rtol = rtol)


  solution2 <- matrix(solution2[,-1], nrow = 2) # remove the time from the result

  #########interval4 [t1, t0]

  interval4 <- function(t, state, parameter) {
    with(as.list(c(state, parameter)), {

      lambdac <- parameter[[1]]
      mu      <- parameter[[2]]
      gamma   <- parameter[[3]]
      lambdaa <- parameter[[4]]
      q       <- parameter[[5]]
      p       <- parameter[[6]]


     # n <- (length(state) - 1) / 2
      n <- num_observed_states * num_hidden_states

      dDM1    <- numeric(n)
      dDE     <- numeric(n)


      t_vec   <- rowSums(q)

      DM1  <- state[1:n]
      E    <- state[(n + 1):(n + n)]
      DA1  <- state[length(state)]

      gamma_matrix <- matrix(parameter[[3]], nrow = n, ncol = length(parameter[[3]]), byrow = TRUE)

      if (nrow(gamma_matrix) == 1) {
        # when there's only one row, there is no “self” element to subtract
        gamma_nonself <- 0
      } else {
        # for n > 1, subtract the diagonal (self‐effects) as before
        gamma_nonself <- rowSums(gamma_matrix - diag(parameter[[3]]))
      }

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

  gamma <- parameter[[3]]

  # only use second row, because the first row of solution4 is the initial state

  #if the trait state of the species at the stem is unknown
  if (trait_mainland_ancestor == "FALSE")
  {
    initial_conditions4 <- c(rep( sum(gamma * (solution2[2,][(n + 1):(n + n)])), n), ### DM1: select DM2 in solution2
                             solution2[2,][(n + n + n + 1):(n + n + n + n)],         ### E: select E in solution2
                             sum(gamma * (solution2[2,][(n + 1):(n + n)])))          ### DA1: select DM2 in solution2

  }
  #if the trait state of the species at the stem is known
  else if(trait_mainland_ancestor == trait_mainland_ancestor)
  {
    pos <- c((num_hidden_states*trait_mainland_ancestor + 1), num_hidden_states + trait_mainland_ancestor* num_hidden_states)
    initial_conditions4 <- c(rep (sum(parameter[[3]][pos] * (solution2[2,][(n + 1):(n + n)])[pos])/2,n ), ### DM1: select DM2 in solution2
                             solution2[2,][(n + n + n + 1):(n + n + n + n)],                                                                       ### E: select E in solution2
                             sum(parameter[[3]][pos] * (solution2[2,][(n + 1):(n + n)])[pos]/2))          ### DA1: select DM2 in solution2

  }


  initial_conditions4 <- matrix(initial_conditions4, nrow = 1)

  # Time sequence for interval [t1, t0]
  time4 <- c(t1, t0)

  # Solve the system for interval [t1, t0]
  solution4 <- deSolve::ode(y = initial_conditions4,
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
