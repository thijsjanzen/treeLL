

get_initial_conditions2 <- function(status, res, trait, num_observed_states, num_hidden_states)

{
  n <- num_observed_states * num_hidden_states
  num_unique_states <- n

  if (status == 2 && length(brts)>2)
  {

    initial_conditions2 <- c(res[1:n],                      ## DE
                             (res[1:n]) * res[length(res)],  ## DM2
                             res[(n + 1):(n + n)],          ## DM3
                             res[(n + n + 1):(n + n + n)],  ## E
                             res[length(res)])              ## DA3

  }
  else if (status == 2 && length(brts) == 2)
  {

    calc_init_state_hidden <- function(trait,
                                       num_unique_states,
                                       num_hidden_states) {

      DE  <- rep(0, num_unique_states)
      DM2 <- rep(0, num_unique_states)
      DM3 <- rep(0, num_unique_states)
      E   <- rep(0, num_unique_states)
      DA3 <- 1


      DE[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- 1

      return( c(DE, DM2, DM3, E, DA3))
    }


    initial_conditions2 <- calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  }


  else if (status == 3 && length(brts) == 2 )
  {

    calc_init_state_hidden <- function(trait,
                                       num_unique_states,
                                       num_hidden_states) {

      DE  <- rep(0, num_unique_states)
      DM2 <- rep(0, num_unique_states)
      DM3 <- rep(0, num_unique_states)
      E   <- rep(0, num_unique_states)
      DA3 <- 0


      DE[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- 1
      DM3[c((num_hidden_states*trait_mainland_ancestor + 1), num_hidden_states + trait_mainland_ancestor* num_hidden_states)] <- 1

      return( c(DE, DM2, DM3, E, DA3))
    }

    initial_conditions2 <- calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  }

  else if (status == 4)
  {
    calc_init_state_hidden <- function(trait,
                                       num_unique_states,
                                       num_hidden_states) {

      DE  <- rep(0, num_unique_states)
      DM2 <- rep(0, num_unique_states)
      DM3 <- rep(0, num_unique_states)
      E   <- rep(0, num_unique_states)
      DA3 <- 1


      DM2[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- 1

      return( c(DE, DM2, DM3, E, DA3))
    }



    initial_conditions2 <- calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  }

  else if(status == 8)

  {
    calc_init_state_hidden <- function(trait,
                                       num_unique_states,
                                       num_hidden_states) {

      DE  <- rep(0, num_unique_states)
      DM2 <- rep(0, num_unique_states)
      DM3 <- rep(0, num_unique_states)
      E   <- rep(0, num_unique_states)
      DA3 <- 1


      DM2[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- 1

      return( c(DE, DM2, DM3, E, DA3))
    }

    initial_conditions2 <-   calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  }
  else if(status == 9)

  {
    calc_init_state_hidden <- function(trait,
                                       num_unique_states,
                                       num_hidden_states) {

      DE  <- rep(0, num_unique_states)
      DM2 <- rep(0, num_unique_states)
      DM3 <- rep(0, num_unique_states)
      E   <- rep(0, num_unique_states)
      DA3 <- 1


      DE[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- 1

      return( c(DE, DM2, DM3, E, DA3))
    }

    initial_conditions2 <-   calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  }
  return(matrix(initial_conditions2, nrow = 1))

}

get_initial_conditions3 <- function(status, res, num_observed_states, num_hidden_states, trait) {

  n <- num_observed_states * num_hidden_states
  num_unique_states <- n

  if (status == 1) {
    calc_init_state_hidden <- function(trait,
                                       num_unique_states,
                                       num_hidden_states) {

      DE  <- rep(0, num_unique_states)
      DM2 <- rep(0, num_unique_states)
      DM1 <- rep(0, num_unique_states)
      DM3 <- rep(0, num_unique_states)
      E   <- rep(0, num_unique_states)
      DA2 <- 0
      DA3 <- 1

      DM2[c((num_hidden_states * trait + 1), num_hidden_states + trait * num_hidden_states)] <- 1

      return(c(DE, DM1, DM2, DM3, E, DA2, DA3))
    }

    initial_conditions3 <- calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  }


  else if (status == 6) {
    initial_conditions3 <- c(res[1:n],                                              ## DE
                             rep(0, n),                                             ## DM1
                             (res[1:n]) * res[length(res)],                         ## DM2
                             res[(n + 1):(n + n)],                                  ## DM3
                             res[(n + n + 1):(n + n + n)],                          ## E
                             0,                                                     ## DA2
                             res[length(res)])                                      ## DA3

  }
  else if (status == 5) {
    calc_init_state_hidden <- function(trait,
                                       num_unique_states,
                                       num_hidden_states) {

      DE  <- rep(0, num_unique_states)
      DM2 <- rep(0, num_unique_states)
      DM1 <- rep(0, num_unique_states)
      DM3 <- rep(0, num_unique_states)
      E   <- rep(0, num_unique_states)
      DA2 <- 0
      DA3 <- 1

      DE[c((num_hidden_states * trait + 1), num_hidden_states + trait * num_hidden_states)] <- 1

      return(c(DE, DM1, DM2, DM3, E, DA2, DA3))
    }

    initial_conditions3 <- calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  } else if (status == 8 || status == 9) {
    initial_conditions3 <- c(solution[2, ][1:n],
                             rep(0, n),                                             ### DE: select DE in solution2
                             solution[2, ][(n + 1):(n + n)],                        ### DM2: select DM2 in solution2
                             solution[2, ][(n + n + 1):(n + n + n)],                ### DM3: select DM3 in solution2
                             solution[2, ][(n + n + n + 1):(n + n + n + n)],        ### E: select E in solution2
                             0,
                             solution[2, ][length(solution[2, ])])                 ### DA3: select DA3 in solution2

  }

  return(matrix(initial_conditions3, nrow = 1))
}




get_initial_conditions4 <- function(status, solution, parameter, trait_mainland_ancestor, num_observed_states, num_hidden_states)

{
  n <- num_observed_states * num_hidden_states
  num_unique_states <- n
  if (status == 2 || status == 3 || status == 4)
  {

    if (trait_mainland_ancestor == "FALSE")
    {  initial_conditions4 <- c(rep( sum(parameter[[3]] * (solution[2,][(n + 1):(n + n)])), n), ### DM1: select DM2 in solution2
                                solution[2,][(n + n + n + 1):(n + n + n + n)],         ### E: select E in solution2
                                sum(parameter[[3]] * (solution[2,][(n + 1):(n + n)])))          ### DA1: select DM2 in solution2

    }
    #if the trait state of the species at the sten is known
    else if(trait_mainland_ancestor == trait_mainland_ancestor)
    {
      pos <- c((num_hidden_states*trait_mainland_ancestor + 1), num_hidden_states + trait_mainland_ancestor* num_hidden_states)
      initial_conditions4 <- c(rep (sum(parameter[[3]][pos] * (solution[2,][(n + 1):(n + n)])[pos])/2,n ), ### DM1: select DM2 in solution2
                               solution[2,][(n + n + n + 1):(n + n + n + n)],                                                                       ### E: select E in solution2
                               sum(parameter[[3]][pos] * (solution[2,][(n + 1):(n + n)])[pos]/2))          ### DA1: select DM2 in solution2

    }
  }
  else if(status == 1 || status == 5 || status == 6)

  {
    initial_conditions4 <- c(solution[2,][(n + 1):(n + n)],                                 ### DM1: select DM1 in solution1
                                 solution[2,][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution1
                                 solution[2,][length(solution[2,]) - 1])                       ### DA1: select DA2 in solution1

  }

  else if(status == 8 || status == 9)

  {
    initial_conditions4 <- c(solution[2,][(n + 1):(n + n)],                                 ### DM1: select DM2 in solution3
                                     solution[2,][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution3
                                     solution[2,][length(solution[2,]) - 1])                       ### DA1: select DA2 in solution3

  }
  return(matrix(initial_conditions4, nrow = 1))

}


func_for_solution <- function(interval,
                              initial_conditions,
                              time,
                              parameter,
                              methode,
                              atol,
                              rtol) {
  interval_func <- get_func_interval(interval)

  if (interval == "interval2") {
    solution <- deSolve::ode(
      y = initial_conditions,
      times = time,
      func = interval_func,
      parms = parameter,
      method = methode,
      atol = atol,
      rtol = rtol
    )
  } else if (interval == "interval3") {
    solution <- deSolve::ode(
      y = initial_conditions,
      times = time,
      func = interval_func,
      parms = parameter,
      method = methode,
      atol = atol,
      rtol = rtol
    )
  } else if (interval == "interval4") {
    solution <- deSolve::ode(
      y = initial_conditions,
      times = time,
      func = interval_func,
      parms = parameter,
      method = methode,
      atol = atol,
      rtol = rtol
    )
  } else {
    stop("Unknown interval name passed.")
  }

  return(matrix(solution[, -1], nrow = 2))
}




get_func_interval <- function(interval) {
  interval_fun <- NULL

  if (interval == "interval2") {
    interval_fun <- function(t, state, parameter) {
      with(as.list(c(state, parameter)), {
        lambdac <- parameter[[1]]
        mu      <- parameter[[2]]
        gamma   <- parameter[[3]]
        lambdaa <- parameter[[4]]
        q       <- parameter[[5]]
        p       <- parameter[[6]]

        n <- num_observed_states * num_hidden_states

        dDE     <- numeric(n)
        dDM2    <- numeric(n)
        dDM3    <- numeric(n)
        dE      <- numeric(n)

        t_vec <- rowSums(q)

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

        dDM3 <- -(lambdac + mu + gamma_nonself + lambdaa + t_vec) * DM3 +
          (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA3 +
          (1 - p) * q_mult_DM3 + gamma_nonself * DM3

        dE <- mu - (mu + lambdac + t_vec) * E +
          lambdac * E * E +
          q_mult_E

        dDA3 <- -sum(gamma) * DA3 + sum(gamma * DM3)

        return(list(c(dDE, dDM2, dDM3, dE, dDA3)))
      })
    }
  } else if (interval == "interval3") {
    interval_fun <- function(t, state, parameter) {
      with(as.list(c(state, parameter)), {
        lambdac <- parameter[[1]]
        mu      <- parameter[[2]]
        gamma   <- parameter[[3]]
        lambdaa <- parameter[[4]]
        q       <- parameter[[5]]
        p       <- parameter[[6]]

        n <- num_observed_states * num_hidden_states

        dDE     <- numeric(n)
        dDM1    <- numeric(n)
        dDM2    <- numeric(n)
        dDM3    <- numeric(n)
        dE      <- numeric(n)

        t_vec <- rowSums(q)

        DE  <- state[1:n]
        DM1 <- state[(n + 1):(n + n)]
        DM2 <- state[(n + n + 1):(n + n + n)]
        DM3 <- state[(n + n + n + 1):(n + n + n + n)]
        E   <- state[(n + n + n + n + 1):(n + n + n + n + n)]
        DA2 <- state[length(state) - 1]
        DA3 <- state[length(state)]

        gamma_matrix  <- matrix(gamma, nrow = n, ncol = length(gamma), byrow = TRUE)
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
          (lambdaa * DE + 2 * lambdac * DE + p * q_mult_DE) * DA3 +
          (1 - p) * q_mult_DM2 + gamma_nonself * DM2

        dDM3 <- -(lambdac + mu + gamma_nonself + lambdaa + t_vec) * DM3 +
          (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA3 +
          (1 - p) * q_mult_DM3 + gamma_nonself * DM3

        dE <- mu - (mu + lambdac + t_vec) * E +
          lambdac * E * E +
          q_mult_E

        dDA2 <- -sum(gamma) * DA2 + sum(gamma * DM2)
        dDA3 <- -sum(gamma) * DA3 + sum(gamma * DM3)

        return(list(c(dDE, dDM1, dDM2, dDM3, dE, dDA2, dDA3)))
      })
    }
  } else if (interval == "interval4") {
    interval_fun <- function(t, state, parameter) {
      with(as.list(c(state, parameter)), {
        lambdac <- parameter[[1]]
        mu      <- parameter[[2]]
        gamma   <- parameter[[3]]
        lambdaa <- parameter[[4]]
        q       <- parameter[[5]]
        p       <- parameter[[6]]

        n <- num_observed_states * num_hidden_states

        dDM1 <- numeric(n)
        dDE  <- numeric(n)

        t_vec <- rowSums(q)

        DM1 <- state[1:n]
        E   <- state[(n + 1):(n + n)]
        DA1 <- state[length(state)]

        gamma_matrix <- matrix(gamma, nrow = n, ncol = length(gamma), byrow = TRUE)
        gamma_nonself <- rowSums(gamma_matrix - diag(gamma))

        q_mult_E   <- t(q %*% E)
        q_mult_DM1 <- t(q %*% DM1)

        dDM1 <- -(lambdac + mu + gamma_nonself + lambdaa + t_vec) * DM1 +
          (mu + lambdaa * E + lambdac * E * E + p * q_mult_E) * DA1 +
          (1 - p) * q_mult_DM1 + gamma_nonself * DM1

        dE <- mu - (mu + lambdac + t_vec) * E +
          lambdac * E * E +
          q_mult_E

        dDA1 <- -sum(gamma) * DA1 + sum(gamma * DM1)

        return(list(c(dDM1, dE, dDA1)))
      })
    }
  }

  return(interval_fun)
}

