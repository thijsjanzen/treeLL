#' testing fuction, for comparison with DAISIE
#' @description
#' This function compute the likelihood that all species that colonize the island
#' have gone extinct prior to the present.
#' @export
#' @param cond conditioning, default = "proper_cond"
#' @param root_state_weight root weight, default = "proper_weights"
#' @param setting_calculation used in ML
#' @param see_ancestral_states recover the ancestral states
#' @param atol absolute tolerance
#' @param rtol relative tolerance
#' @param methode method of integration
#' @examples
#' Load DAISIE package and data
#' library(DAISIE)
#' data("Galapagos_datalist")
#'
#' datalist <- Galapagos_datalist
#' parameter <- list(
#'   c(2.546591, 0),        # cladogenesis rate & zero for hidden state
#'   c(2.678781, 0),        # extinction rate & zero for hidden state
#'   c(0.009326754, 0),     # trait‐change rate & zero for hidden state
#'   c(1.008583, 0),        # anagenesis rate & zero for hidden state
#'   matrix(rep(0, 4), nrow = 2),  # transition matrix Q
#'   0                      # probability p
#' )
#'
#' # Compute log‐likelihood under the DE‐trait model
#' DAISIE_DE_trait_logp0(
#'   datalist            = datalist,
#'   parameter           = parameter,
#'   cond                = "proper_cond",
#'   root_state_weight   = "proper_weights",
#'   see_ancestral_states= TRUE,
#'   atol                = 1e-10,
#'   rtol                = 1e-10,
#'   methode             = "ode45"
#' )


DAISIE_DE_trait_logp0 <- function(datalist,
                                  parameter,
                                  cond = "proper_cond",
                                  root_state_weight = "proper_weights",
                                  see_ancestral_states = TRUE,
                                  atol = 1e-10,
                                  rtol = 1e-10,
                                  methode = "ode45") {


  t0 <- datalist[[1]]$island_age
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
      # number of unique state


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
                            atol = atol,
                            rtol = rtol)



  solution1 <- matrix(solution1[,-1], nrow = 2)

  # Extract log-likelihood
  Lk <- solution1[2,][length(solution1[2,])]
  logLkb <- log(Lk)
  return(logLkb)
}
