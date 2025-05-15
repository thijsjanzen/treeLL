#' testing fuction, for comparison with DAISIE
#' @description
#' This function calculates something we can verify with DAISIE
#' @export
#' @param brts branching times
#' @param missnumspec number of missing species
#' @param parameter parameters
#' @param phy phy
#' @param traits traits
#' @param cond conditioning, default = "proper_cond"
#' @param root_state_weight root weight, default = "proper_weights"
#' @param setting_calculation used in ML
#' @param see_ancestral_states recover the ancestral states
#' @param atol absolute tolerance
#' @param rtol relative tolerance
#' @param methode method of integration
DAISIE_DE_logp0_trait <- function(parameter,
                                  cond = "proper_cond",
                                  root_state_weight = "proper_weights",
                                  see_ancestral_states = TRUE,
                                  atol = 1e-10,
                                  rtol = 1e-10,
                                  methode = "ode45") {



  #########Interval1 [t_p, t_0]
  interval1 <- function(t, state, parameter) {
    with(as.list(c(state, parameter)), {

      lambdac <- parameter[[1]]
      mu      <- parameter[[2]]
      gamma   <- parameter[[3]]
      lambdaa <- parameter[[4]]
      q       <- parameter[[5]]
      p       <- parameter[[6]]


      n <- (length(state) - 1) / 2


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

  initial_conditions1 <- c(rep(0, num_unique_states), ### DM1
                           rep(0, num_unique_states), ### E
                           1)                         ### DA1



  # Time sequence for interval [tp, t0]
  time1 <- c(tp, t0)

  # Solve the system for interval [tp, t1]
  solution1 <- deSolve::ode(y = initial_conditions1,
                            times = time1,
                            func = interval1,
                            parms = parameter,
                            method = methode,
                            atol = 1e-10,
                            rtol = 1e-10)



  solution1 <- matrix(solution1[,-1], nrow = 2)

  # Extract log-likelihood
  Lk <- solution1[2,][length(solution1[2,])]
  logLkb <- log(Lk)
  return(logLkb)
}

library (DAISIE)
data("Galapagos_datalist")
datalist <- Galapagos_datalist
i <- 3
parameter <- list( c(2.546591,0), c(2.678781, 0), c(0.009326754, 0), c(1.008583, 0), matrix(rep(0,4), nrow = 2), 0)

DAISIE_DE_logp0_trait(parameter = parameter,
                      cond = "proper_cond",
                      root_state_weight = "proper_weights",
                      see_ancestral_states = TRUE,
                      atol = 1e-10,
                      rtol = 1e-10,
                      methode = "ode45")

