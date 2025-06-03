#' testing fuction, for comparison with DAISIE
#' @description
#' This function compute the likelihood that all species that colonize the island
#' have gone extinct prior to the present.
#' @export
#' @inheritParams default_params_doc
#' @examples
#' #Load DAISIE package and data
#' library(DAISIE)
#' data("Galapagos_datalist")
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
#' # Compute log‐likelihood under the DE‐trait model
#' DAISIE_DE_trait_logp0(
#'   datalist            = datalist,
#'   parameter           = parameter,
#'   num_observed_states     = 2,
#'   num_hidden_states       = 2,
#'   atol                = 1e-10,
#'   rtol                = 1e-10,
#'   methode             = "ode45"
#' )


DAISIE_DE_trait_logp0 <- function(datalist,
                                  parameter,
                                  atol = 1e-10,
                                  rtol = 1e-10,
                                  num_observed_states,
                                  num_hidden_states,
                                  methode = "ode45") {

  n <- num_observed_states * num_hidden_states
  t0 <- datalist[[1]]$island_age
  tp <- 0
  #########interval4 [t_p, t_0]

  interval4_local <- function(t, state, parameter) {
    with(as.list(c(state, parameter)), {

      lambdac <- parameter[[1]]
      mu      <- parameter[[2]]
      gamma   <- parameter[[3]]
      lambdaa <- parameter[[4]]
      q       <- parameter[[5]]
      p       <- parameter[[6]]


      n <- num_observed_states * num_hidden_states
      # number of unique state


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

  initial_conditions40 <- c(rep(0, n), ### DM1
                           rep(0, n), ### E
                           1)                         ### DA1



  # Time sequence for interval [tp, t0]
  time4 <- c(tp, t0)

  # Solve the system for interval [tp, t1]
  solution4 <- deSolve::ode(y = initial_conditions40,
                            times = time4,
                            func = interval4_local,
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
