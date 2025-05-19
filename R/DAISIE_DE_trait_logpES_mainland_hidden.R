#' testing fuction, for comparison with DAISIE
#' @description
#' this function calculates the likelihood of observing a singleton endemic species on an island
#' with the trait state `i`, that coexist on the island with its mainland ancestors.
#' @export
#' @param brts branching times
#' @param parameter parameters
#' @param num_observed_states number of observed traits
#' @param num_hidden_states number of hidden traits
#' @param trait trait of the species at the tip
#' @param trait_mainland_ancestor trait of the mainland ancestors
#' @param cond conditioning, default = "proper_cond"
#' @param root_state_weight root weight, default = "proper_weights"
#' @param setting_calculation used in ML
#' @param see_ancestral_states recover the ancestral states
#' @param atol absolute tolerance
#' @param rtol relative tolerance
#' @param methode method of integration
DAISIE_DE_trait_logpES_mainland_hidden <- function(brts,
                                                  trait,
                                                  parameter,
                                                  num_observed_states,
                                                  num_hidden_states,
                                                  trait_mainland_ancestor = trait_mainland_ancestor,
                                                  cond = "proper_cond",
                                                  root_state_weight = "proper_weights",
                                                  see_ancestral_states = TRUE,
                                                  atol = 1e-10,
                                                  rtol = 1e-10,
                                                  methode = "ode45") {
  t0 <- brts[1]
  t1 <- brts[2]
  # number of unique state
  n <- num_observed_states * num_hidden_states

  #########Interval1 [t_p, t_1]

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




  calc_init_state_hidden <- function(trait,
                                     num_unique_states,
                                     num_hidden_states) {

    DE  <- rep(0, num_unique_states)
    DM2 <- rep(0, num_unique_states)
    DM3 <- rep(0, num_unique_states)
    E   <- rep(0, num_unique_states)
    DA3 <- 0

    #for (i in 1:num_hidden_states) {
    #assuming the traits start counting at 0 !!!!
    #DM2[(1 + trait) + (i - 1) * num_hidden_states] <- 1
    #}
    DE[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- 1
    DM3[c((num_hidden_states*trait_mainland_ancestor + 1), num_hidden_states + trait_mainland_ancestor* num_hidden_states)] <- 1

    return( c(DE, DM2, DM3, E, DA3))
  }





  num_unique_states <- length(parameter[[1]])
  initial_conditions1 <-   calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  initial_conditions1 <- matrix(initial_conditions1, nrow = 1)




  # Time sequence for interval [tp, t1]
  time1 <- c(tp, t1)

  # Solve the system for interval [tp, t1]
  solution1 <- deSolve::ode(y = initial_conditions1,
                            times = time1,
                            func = interval1,
                            parms = parameter,
                            method = methode,
                            atol = atol,
                            rtol = rtol)


  solution1 <- matrix(solution1[,-1], nrow = 2) # remove the time from the result

  #########Interval2 [t1, t0]

  interval2 <- function(t, state, parameter) {
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
  gamma <- parameter[[3]]

  # only use second row, because the first row of solution2 is the initial state

  #if the trait state of the species at the stem is unknown
  if (trait_mainland_ancestor == "FALSE")
  {
    initial_conditions2 <- c(rep( sum(gamma * (solution1[2,][(n + 1):(n + n)])), n), ### DM1: select DM2 in solution1
                             solution1[2,][(n + n + n + 1):(n + n + n + n)],         ### E: select E in solution1
                             sum(gamma * (solution1[2,][(n + 1):(n + n)])))          ### DA1: select DA3 in solution1

  }
  #if the trait state of the species at the stem is known
  else if(trait_mainland_ancestor == trait_mainland_ancestor)
  {
    initial_conditions2 <- c(rep (parameter[[3]][trait_mainland_ancestor + 1] * (solution1[2,][(n + 1):(n + n)])[trait_mainland_ancestor + 1], n), ### DM1: select DM2 in solution1
                             solution1[2,][(n + n + n + 1):(n + n + n + n)],         ### E: select E in solution1
                             parameter[[3]][trait_mainland_ancestor + 1] * (solution1[2,][(n + 1):(n + n)])[trait_mainland_ancestor + 1])          ### DA1: select DA3 in solution1

  }


  initial_conditions2 <- matrix(initial_conditions2, nrow = 1)

  # Time sequence for interval [t1, t0]
  time2 <- c(t1, t0)

  # Solve the system for interval [t1, t0]
  solution2 <- deSolve::ode(y = initial_conditions2,
                            times = time2,
                            func = interval2,
                            parms = parameter,
                            method = methode,
                            atol = atol,
                            rtol = rtol)

  solution2 <- matrix(solution2[,-1], nrow = 2)

  # Extract log-likelihood
  Lk <- solution2[2,][length(solution2[2,])]
  logLkb <- log(Lk)
  return(logLkb)
}
