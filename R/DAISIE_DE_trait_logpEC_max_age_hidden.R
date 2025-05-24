#' Testing function for comparison with DAISIE
#'
#' @description
#' This function calculates the likelihood of observing a clade with specified species trait states,
#' for which the estimated maximum age of colonization is known.
#'
#' @param brts Branching times.
#' @param missnumspec Number of missing species.
#' @param parameter Model parameters.
#' @param num_observed_states Number of observed trait states.
#' @param num_hidden_states Number of hidden trait states.
#' @param phy Phylogeny (of class `phylo`).
#' @param traits Trait states for the species.
#' @param cond Conditioning scheme (default = "proper_cond").
#' @param root_state_weight Root state weighting method (default = "proper_weights").
#' @param setting_calculation Setting used in maximum likelihood estimation.
#' @param see_ancestral_states Logical; whether to recover ancestral states.
#' @param atol Absolute tolerance for ODE integration.
#' @param rtol Relative tolerance for ODE integration.
#' @param methode Method used for numerical integration (e.g., "ode45").
#'
#' @export
#'
#' @examples
#' library(DAISIE)
#' data("NewZealand_birds_datalist")
#' datalist <- NewZealand_birds_datalist
#' i <- 23
#' phy <- DDD::brts2phylo(datalist[[i]]$branching_times[-c(1, 2)])
#' traits <- sample(c(0, 1), length(phy$tip.label), replace = TRUE)
#' sampling_fraction <- sample(c(1, 1), length(phy$tip.label), replace = TRUE)
#'
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
#' DAISIE_DE_trait_logpEC_max_age_hidden(
#'   brts                  = datalist[[i]]$branching_times,
#'   phy                   = phy,
#'   traits                = traits,
#'   sampling_fraction     = sampling_fraction,
#'   parameter             = parameter,
#'   num_observed_states   = 2,
#'   num_hidden_states     = 2,
#'   cond                  = "proper_cond",
#'   root_state_weight     = "proper_weights",
#'   see_ancestral_states  = TRUE,
#'   atol                  = 1e-10,
#'   rtol                  = 1e-10,
#'   methode               = "ode45",
#'   rhs_func              = loglik_hidden_rhs
#' )

DAISIE_DE_trait_logpEC_max_age_hidden <- function(brts,
                                                   parameter,
                                                   phy,
                                                   traits,
                                                   num_observed_states,
                                                   num_hidden_states,
                                                   sampling_fraction,
                                                   cond = "proper_cond",
                                                   root_state_weight = "proper_weights",
                                                   see_ancestral_states = TRUE,
                                                   atol = 1e-10,
                                                   rtol = 1e-10,
                                                   methode = "ode45",
                                                   rhs_func = loglik_hidden_rhs) {
  t0   <- brts[1]
  tmax <- brts[2]
  t2   <- brts[3]
  tp   <- 0

  # number of unique state
  n <- num_observed_states*num_hidden_states

  #########interval3 [t_2, tmax]

  interval3 <- function(t, state, parameter) {
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


  # Solve the system for interval [tp, t2]
  res <- treeLL::loglik_R_hidden(parameter,
                                 phy,
                                 traits,
                                 sampling_fraction,
                                 num_hidden_states = num_hidden_states,
                                 see_ancestral_states = TRUE,
                                 atol = atol,
                                 rtol = rtol)



  initial_conditions3 <- c(res[1:n],                                              ## DE
                           rep(0,n),                                              ## DM1
                           (res[1:n]) * res[length(res)],                         ## DM2
                           res[(n + 1):(n + n)],                                  ## DM3
                           res[(n + n + 1):(n + n + n)],                          ## E
                           0,                                                     ## DA2
                           res[length(res)])                                      ## DA3

  initial_conditions3 <- matrix(initial_conditions3, nrow = 1)




  # Time sequence for interval [t2, tmax]
  time3 <- c(t2, tmax)

  # Solve the system for interval [t2, tmax]
  solution3 <- deSolve::ode(y = initial_conditions3,
                            times = time3,
                            func = interval3,
                            parms = parameter,
                            method = methode,
                            atol = atol,
                            rtol = rtol)


  solution3 <- matrix(solution3[,-1], nrow = 2) # remove the time from the result

  #########interval4 [tmax, t0]

  interval4 <- function(t, state, parameter) {
    with(as.list(c(state, parameter)), {

      lambdac <- parameter[[1]]
      mu      <- parameter[[2]]
      gamma   <- parameter[[3]]
      lambdaa <- parameter[[4]]
      q       <- parameter[[5]]
      p       <- parameter[[6]]


      #n <- (length(state) - 1) / 2
      n <- num_observed_states*num_hidden_states

      dDM1    <- numeric(n)
      dDE     <- numeric(n)


      t_vec   <- rowSums(q)

      DM1  <- state[1:n]
      E   <- state[(n + 1):(n + n)]
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


  # only use second row, because the first row of solution3 is the initial state
 initial_conditions4_max <- c(solution3[2,][(n + 1):(n + n)],                                 ### DM1: select DM1 in solution1
                           solution3[2,][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution1
                           solution3[2,][length(solution3[2,]) - 1])                       ### DA1: select DA2 in solution1


 initial_conditions4_max <- matrix(initial_conditions4_max, nrow = 1)

  # Time sequence for interval [tmax, t0]
  time4 <- c(tmax, t0)

  # Solve the system for interval [tmax, t0]
  solution4 <- deSolve::ode(y =initial_conditions4_max,
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

