#' @name DAISIE_DE_trait_n
#' @export


trait_likelihood_arguments <- function(DAISIE_DE_trait_function,
                                       brts,
                                       traits,
                                       status,
                                       parameter,
                                       num_observed_states,
                                       num_hidden_states) {

  fun <- match.fun(DAISIE_DE_trait_function)
  accepted_arguments <- names(formals(fun))

  arguments <- list(
    brts = brts,
    status = status,
    parameter = parameter,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states
  )

  # logpEC uses `traits`; logpES and logpNE use `trait`
  if ("traits" %in% accepted_arguments) {
    arguments$traits <- traits
  } else {
    arguments$trait <- traits
  }

  arguments <- c(arguments, list(...))

  if (!"..." %in% accepted_arguments) {
    arguments <- arguments[names(arguments) %in% accepted_arguments]
  }

  list(fun = fun, arguments = arguments)
}


all_indices <- function(maximum) {
  unname(as.matrix(
    expand.grid(lapply(maximum, function(x) 0:x))
  ))
}


mixed_derivative <- function(f, orders, step_size = NULL) {

  derivative_order <- sum(orders)

  if (derivative_order == 0) {
    return(f(rep(0, length(orders))))
  }

  if (is.null(step_size)) {
    step_size <- .Machine$double.eps^(1 / (derivative_order + 2))
  }

  indices <- all_indices(orders)

  coefficients <- apply(
    indices,
    1,
    function(index) prod((-1)^index * choose(orders, index))
  )

  points <- sweep(-indices, 2, orders / 2, "+") * step_size
  values <- apply(points, 1, f)
  terms <- coefficients * values

  sum(terms[order(abs(terms))]) / step_size^derivative_order
}


# Calculate derivatives of exp(g) from derivatives of g.
# This is the multivariate Bell-polynomial recursion.
derivative_exp_from_log <- function(orders, log_derivatives) {

  bell <- array(0, dim = orders + 1)

  # exp(g(0) - g(0)) = 1
  bell[matrix(rep(1, length(orders)), nrow = 1)] <- 1

  indices <- all_indices(orders)

  for (row in seq_len(nrow(indices))) {
    k <- indices[row, ]

    if (all(k == 0)) {
      next
    }

    dimension <- which(k > 0)[1]

    previous <- k
    previous[dimension] <- previous[dimension] - 1

    j_indices <- all_indices(previous)
    terms <- numeric(nrow(j_indices))

    for (i in seq_len(nrow(j_indices))) {
      j <- j_indices[i, ]

      log_derivative_index <- j
      log_derivative_index[dimension] <-
        log_derivative_index[dimension] + 1

      terms[i] <-
        prod(choose(previous, j)) *
        bell[matrix(previous - j + 1, nrow = 1)] *
        log_derivatives[
          matrix(log_derivative_index + 1, nrow = 1)
        ]
    }

    bell[matrix(k + 1, nrow = 1)] <- sum(terms[order(abs(terms))])
  }

  bell[matrix(orders + 1, nrow = 1)]
}


DAISIE_DE_trait_n <- function(DAISIE_DE_trait_function,
                              brts,
                              missnumspec,
                              traits,
                              status,
                              parameter,
                              num_observed_states,
                              num_hidden_states,
                              S = NULL,
                              step_size = NULL,
                              ...) {

  if (length(missnumspec) != num_observed_states) {
    stop("missnumspec must have one value per observed trait state.")
  }

  missnumspec <- as.integer(round(missnumspec))

  if (any(missnumspec < 0)) {
    stop("missnumspec cannot be negative.")
  }

  if (anyNA(traits) && sum(missnumspec) > 0) {
    stop("Known trait states are required when species are missing.")
  }

  # S_i = number of sampled species in state i
  if (is.null(S)) {
    S <- tabulate(traits + 1, nbins = num_observed_states)
  }

  likelihood <- trait_likelihood_arguments(
    DAISIE_DE_trait_function = DAISIE_DE_trait_function,
    brts = brts,
    traits = traits,
    status = status,
    parameter = parameter,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states,
    ...
  )

  # g(x) = log f(x), where x_i = 1 - rho_i
  log_f <- function(x) {
    sampling_fraction <- 1 - x

    result <- do.call(
      likelihood$fun,
      c(
        list(sampling_fraction = sampling_fraction),
        likelihood$arguments
      )
    )

    loglik <- if (is.list(result)) result$loglik else result

    # Remove the rho_i^S_i factors
    loglik - sum(S * log(sampling_fraction))
  }

  # Avoid solving the same ODE more than once
  cache <- new.env(hash = TRUE, parent = emptyenv())

  log_f_cached <- function(x) {
    key <- paste(sprintf("%.17g", x), collapse = "|")

    if (!exists(key, envir = cache, inherits = FALSE)) {
      assign(key, log_f(x), envir = cache)
    }

    get(key, envir = cache, inherits = FALSE)
  }

  # Derivatives of g = log(f)
  log_derivatives <- array(0, dim = missnumspec + 1)

  indices <- all_indices(missnumspec)

  for (row in seq_len(nrow(indices))) {
    orders <- indices[row, ]

    if (!all(orders == 0)) {
      log_derivatives[matrix(orders + 1, nrow = 1)] <-
        suppressWarnings(
          mixed_derivative(
            f = log_f_cached,
            orders = orders,
            step_size = step_size
          )
        )
    }
  }

  # f(0) is kept in log space
  log_f_at_zero <- log_f_cached(rep(0, num_observed_states))

  # Bell polynomial: derivative of exp(log_f - log_f_at_zero)
  scaled_derivative <- derivative_exp_from_log(
    orders = missnumspec,
    log_derivatives = log_derivatives
  )

  if (is.na(scaled_derivative) || scaled_derivative <= 0) {
    warning("Numerical derivative is not positive. Returning -Inf.")

    return(list(
      loglik = -Inf,
      deriv = NA,
      g_derivs = log_derivatives,
      S = S
    ))
  }

  loglikelihood <-
    log_f_at_zero +
    log(scaled_derivative) +
    sum(lfactorial(S)) -
    sum(lfactorial(S + missnumspec))

  list(
    loglik = loglikelihood,
    deriv = exp(log_f_at_zero) * scaled_derivative,
    g_derivs = log_derivatives,
    S = S
  )
}
