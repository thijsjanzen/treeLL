count_species_per_trait <- function(traits, num_observed_states) {
  tabulate(traits + 1, nbins = num_observed_states)
}


# Endemic singleton
make_log_f_ES <- function(brts,
                          trait,
                          status,
                          parameter,
                          num_observed_states,
                          num_hidden_states,
                          S,
                          trait_mainland_ancestor) {

  function(x) {
    rho <- 1 - x

    loglik <- DAISIE_DE_trait_logpES(
      brts = brts,
      sampling_fraction = rho,
      trait = trait,
      status = status,
      parameter = parameter,
      num_observed_states = num_observed_states,
      num_hidden_states = num_hidden_states,
      trait_mainland_ancestor = trait_mainland_ancestor
    )$loglik

    loglik - sum(S * log(rho))
  }
}


# Non-endemic lineage
make_log_f_NE <- function(brts,
                          trait,
                          status,
                          parameter,
                          num_observed_states,
                          num_hidden_states,
                          S,
                          trait_mainland_ancestor) {

  function(x) {
    rho <- 1 - x

    loglik <- DAISIE_DE_trait_logpNE(
      brts = brts,
      sampling_fraction = rho,
      trait = trait,
      status = status,
      parameter = parameter,
      num_observed_states = num_observed_states,
      num_hidden_states = num_hidden_states,
      trait_mainland_ancestor = trait_mainland_ancestor
    )$loglik

    loglik - sum(S * log(rho))
  }
}


# Endemic clade
make_log_f_EC <- function(brts,
                          traits,
                          status,
                          parameter,
                          num_observed_states,
                          num_hidden_states,
                          S,
                          trait_mainland_ancestor,
                          phy,
                          num_threads = 1) {

  function(x) {
    rho <- 1 - x

    loglik <- DAISIE_DE_trait_logpEC(
      brts = brts,
      sampling_fraction = rho,
      traits = traits,
      status = status,
      parameter = parameter,
      num_observed_states = num_observed_states,
      num_hidden_states = num_hidden_states,
      trait_mainland_ancestor = trait_mainland_ancestor,
      phy = phy,
      num_threads = num_threads
    )$loglik

    loglik - sum(S * log(rho))
  }
}


all_multi_indices <- function(maximum) {
  unname(as.matrix(
    expand.grid(lapply(maximum, function(x) 0:x))
  ))
}


mixed_partial_derivative <- function(f, orders, step_size = NULL) {

  derivative_order <- sum(orders)

  if (derivative_order == 0) {
    return(f(rep(0, length(orders))))
  }

  if (is.null(step_size)) {
    step_size <- .Machine$double.eps^(
      1 / (derivative_order + 2)
    )
  }

  indices <- all_multi_indices(orders)

  coefficients <- apply(
    indices,
    1,
    function(index) prod((-1)^index * choose(orders, index))
  )

  points <- sweep(-indices, 2, orders / 2, "+") * step_size
  values <- apply(points, 1, f)

  sum(coefficients * values) / step_size^derivative_order
}


bell_polynomial_derivative <- function(orders, log_derivatives) {

  bell <- array(0, dim = orders + 1)
  bell[matrix(rep(1, length(orders)), nrow = 1)] <- 1

  indices <- all_multi_indices(orders)

  for (row in seq_len(nrow(indices))) {
    k <- indices[row, ]

    if (all(k == 0)) {
      next
    }

    dimension <- which(k > 0)[1]

    previous <- k
    previous[dimension] <- previous[dimension] - 1

    j_indices <- all_multi_indices(previous)
    terms <- numeric(nrow(j_indices))

    for (i in seq_len(nrow(j_indices))) {
      j <- j_indices[i, ]

      log_index <- j
      log_index[dimension] <- log_index[dimension] + 1

      terms[i] <-
        prod(choose(previous, j)) *
        bell[matrix(previous - j + 1, nrow = 1)] *
        log_derivatives[matrix(log_index + 1, nrow = 1)]
    }

    bell[matrix(k + 1, nrow = 1)] <- sum(terms)
  }

  bell[matrix(orders + 1, nrow = 1)]
}


n_sampling_derivative <- function(log_f, missnumspec, step_size = NULL) {

  log_derivatives <- array(0, dim = missnumspec + 1)
  indices <- all_multi_indices(missnumspec)

  for (row in seq_len(nrow(indices))) {
    orders <- indices[row, ]

    if (!all(orders == 0)) {
      log_derivatives[matrix(orders + 1, nrow = 1)] <-
        mixed_partial_derivative(
          f = log_f,
          orders = orders,
          step_size = step_size
        )
    }
  }

  list(
    log_f_at_zero = log_f(rep(0, length(missnumspec))),
    derivative = bell_polynomial_derivative(
      orders = missnumspec,
      log_derivatives = log_derivatives
    ),
    log_derivatives = log_derivatives
  )
}


#' Convert one rho-sampling lineage likelihood into an n-sampling likelihood
#'
#' @name DAISIE_DE_trait_n
#' @export
DAISIE_DE_trait_n <- function(log_f,
                              missnumspec,
                              S,
                              step_size = NULL) {

  missnumspec <- as.integer(round(missnumspec))

  if (length(missnumspec) != length(S)) {
    stop("missnumspec and S must have the same length.")
  }

  if (any(missnumspec < 0)) {
    stop("missnumspec cannot be negative.")
  }

  result <- n_sampling_derivative(
    log_f = log_f,
    missnumspec = missnumspec,
    step_size = step_size
  )

  if (is.na(result$derivative) || result$derivative <= 0) {
    warning("Numerical derivative is not positive. Returning -Inf.")

    return(list(
      loglik = -Inf,
      deriv = NA,
      g_derivs = result$log_derivatives,
      S = S
    ))
  }

  loglikelihood <-
    result$log_f_at_zero +
    log(result$derivative) +
    sum(lfactorial(S)) -
    sum(lfactorial(S + missnumspec))

  list(
    loglik = loglikelihood,
    deriv = exp(result$log_f_at_zero) * result$derivative,
    g_derivs = result$log_derivatives,
    S = S
  )
}
