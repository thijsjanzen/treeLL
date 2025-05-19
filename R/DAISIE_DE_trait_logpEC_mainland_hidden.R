#' testing fuction, for comparison with DAISIE
#' @description
#' this function calculates the likelihood of observing a clade on an island
#' that coexists with its mainland ancestors.
#' @export
#' @param brts branching times
#' @param missnumspec number of missing species
#' @param parameter parameters
#' @param num_observed_states number of observed traits
#' @param num_hidden_states number of hidden traits
#' @param trait_mainland_ancestor trait of the species at the stem
#' @param phy phy
#' @param traits traits
#' @param cond conditioning, default = "proper_cond"
#' @param root_state_weight root weight, default = "proper_weights"
#' @param setting_calculation used in ML
#' @param see_ancestral_states recover the ancestral states
#' @param atol absolute tolerance
#' @param rtol relative tolerance
#' @param methode method of integration
DAISIE_DE_trait_logpEC_mainland_hidden <- function(brts,
                                           missnumspec,
                                           parameter,
                                           phy,
                                           traits,
                                           trait_mainland_ancestor = trait_mainland_ancestor,
                                           num_observed_states,
                                           num_hidden_states,
                                           cond = "proper_cond",
                                           root_state_weight = "proper_weights",
                                           see_ancestral_states = TRUE,
                                           atol = 1e-10,
                                           rtol = 1e-10,
                                           methode = "ode45",
                                           rhs_func = loglik_hidden_rhs) {
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0

  # number of unique state
  n  <- num_observed_states * num_hidden_states

  #########Interval2 [t_2, t_1]

  interval2 <- function(t, state, parameter) {
    with(as.list(c(state, parameter)), {

      lambdac <- parameter[[1]]
      mu      <- parameter[[2]]
      gamma   <- parameter[[3]]
      lambdaa <- parameter[[4]]
      q       <- parameter[[5]]
      p       <- parameter[[6]]


      n <- num_observed_states*num_hidden_states

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


  # Solve the system for interval [tp, t2]
  res <- treeLL::loglik_R_mainland_hidden(parameter,
                                 phy,
                                 traits,
                                 num_hidden_traits = num_hidden_traits,
                                 see_ancestral_states = TRUE,
                                 atol = atol,
                                 rtol = rtol)


  initial_conditions2 <- c(res[1:n],                      ## DE
                           (res[1:n]) * res[length(res)], ## DM2
                           res[(n + 1):(n + n)],          ## DM3
                           res[(n + n + 1):(n + n + n)],  ## E
                           res[length(res)])              ## DA3

  initial_conditions2 <- matrix(initial_conditions2, nrow = 1)




  # Time sequence for interval [t2, t1]
  time2 <- c(t2, t1)

  # Solve the system for interval [t2, t1]
  solution2 <- deSolve::ode(y = initial_conditions2,
                            times = time2,
                            func = interval2,
                            parms = parameter,
                            method = methode,
                            atol = atol,
                            rtol = rtol)


  solution2 <- matrix(solution2[,-1], nrow = 2) # remove the time from the result

  #########Interval3 [t1, t0]

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
  gamma <- parameter[[3]]

  # only use second row, because the first row of solution2 is the initial state

  #if the trait state of the species at the stem is unknown
  if (trait_mainland_ancestor == "FALSE")
  {  initial_conditions3 <- c(rep( sum(gamma * (solution2[2,][(m + 1):(m + m)])), m), ### DM1: select DM2 in solution2
                              solution2[2,][(m + m + m + 1):(m + m + m + m)],         ### E: select E in solution2
                              sum(gamma * (solution2[2,][(m + 1):(m + m)])))          ### DA1: select DA3 in solution2

  }
  #if the trait state of the species at the stem is known
  else if(trait_mainland_ancestor == trait_mainland_ancestor)
  {
    initial_conditions3 <- c(rep (parameter[[3]][trait_mainland_ancestor + 1] * (solution2[2,][(m + 1):(m + m)])[trait_mainland_ancestor + 1], m), ### DM1: select DM2 in solution2
                             solution2[2,][(m + m + m + 1):(m + m + m + m)],                                                                       ### E: select E in solution2
                             parameter[[3]][trait_mainland_ancestor + 1] * (solution2[2,][(m + 1):(m + m)])[trait_mainland_ancestor + 1])          ### DA1: select DA3 in solution2

  }

  initial_conditions3 <- matrix(initial_conditions3, nrow = 1)

  # Time sequence for interval [t1, t0]
  time3 <- c(t1, t0)

  # Solve the system for interval [t1, t0]
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

