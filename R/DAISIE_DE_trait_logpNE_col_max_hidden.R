#' testing function, for comparison with DAISIE
#' @description
#' this function calculates the likelihood of observing a singleton species on an island
#' with the trait state `i`, either non-endemic or rendered endemic by a trait change, and
#' for which only the estimated colonization time is unknown.
#' @export
#' @param brts branching times
#' @param missnumspec number of missing species
#' @param parameter parameters
#' @param num_observed_states number of observed traits
#' @param num_hidden_states number of hidden traits
#' @param phy phy
#' @param traits traits
#' @param cond conditioning, default = "proper_cond"
#' @param root_state_weight root weight, default = "proper_weights"
#' @param setting_calculation used in ML
#' @param see_ancestral_states recover the ancestral states
#' @param atol absolute tolerance
#' @param rtol relative tolerance
#' @param methode method of integration
#' @examples
#' library(DAISIE)
#' data("Galapagos_datalist")
#' datalist <- Galapagos_datalist
#' i <- 2
#'
#' parameter <- list(
#'   c(2.546591, 1.2, 1, 0.2),
#'   c(2.678781, 2, 1.9, 3),
#'   c(0.009326754, 0.003, 0.002, 0.2),
#'   c(1.008583, 1, 2, 1.5),
#'   matrix(
#'     c(0,   1,    0.5, 0,
#'       0,   0,    0,   0,
#'       0.002, 0.005, rep(0, 8)),
#'     nrow = 4
#'   ),
#'   0
#' )
#'
#' DAISIE_DE_trait_logpNE_col_max_hidden(
#'   brts                  = datalist[[i]]$branching_times,
#'   trait                 = 0,
#'   parameter <- parameter,
#'   num_observed_states   = 2,
#'   num_hidden_states     = 2,
#'   cond                  = "proper_cond",
#'   root_state_weight     = "proper_weights",
#'   see_ancestral_states  = TRUE,
#'   atol                  = 1e-10,
#'   rtol                  = 1e-10,
#'   methode               = "ode45"
#' )


DAISIE_DE_trait_logpNE_col_max_hidden <- function(brts,
                                                  trait,
                                                  parameter,
                                                  num_hidden_states,
                                                  num_observed_states = 1,
                                                  cond = "proper_cond",
                                                  root_state_weight = "proper_weights",
                                                  see_ancestral_states = TRUE,
                                                  atol = 1e-10,
                                                  rtol = 1e-10,
                                                  methode = "ode45") {


  t0 <- brts[1]
  tp <- 0
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

  calc_init_state_hidden <- function(trait,
                                     num_unique_states,
                                     num_hidden_states) {

    DM1 <- rep(0, num_unique_states)
    E   <- rep(0, num_unique_states)
    DA1 <- 0

    #for (i in 1:num_hidden_states) {
      # assuming the traits start counting at 0 !!!!
#DM1[(1 + trait) + (i - 1) * num_hidden_states] <- 1
   # }

    DM1[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- 1

    return( c(DM1, E, DA1))
  }


  num_unique_states <- length(parameter[[1]])
  initial_conditions1 <- calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  initial_conditions1 <- matrix(initial_conditions1, nrow = 1)



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


