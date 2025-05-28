#' @keywords internal
interval2 <- function(t, state, parameter) {
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

#' @keywords internal
interval3 <- function(t, state, parameter) {
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

#' @keywords internal
interval4 <- function(t, state, parameter) {
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
