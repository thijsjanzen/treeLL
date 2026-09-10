#' @name DAISIE_DE_trait_n
#' @export
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
    stop(
      "missnumspec must have one value per observed trait state."
    )
  }

  missnumspec <- as.integer(round(missnumspec))


  # We can relax this later.
  if (anyNA(traits) && sum(missnumspec) > 0) {
    stop(
      "DAISIE_DE_trait_n requires known trait states ",
      "when species are missing."
    )
  }

  # Number of sampled species in each observed trait state
  if (is.null(S)) {
    S <- tabulate(
      traits + 1,
      nbins = num_observed_states
    )
  }

  if (length(S) != num_observed_states) {
    stop("S must have one value per observed trait state.")
  }

  # Prepare arguments for logpEC, logpES, or logpNE
  fun <- match.fun(DAISIE_DE_trait_function)
  fun_args <- names(formals(fun))
  has_dots <- "..." %in% fun_args

  fixed_args <- list(
    brts = brts,
    status = status,
    parameter = parameter,
    num_observed_states = num_observed_states,
    num_hidden_states = num_hidden_states
  )

  # logpEC uses `traits`; logpES and logpNE use `trait`
  if ("traits" %in% fun_args) {
    fixed_args$traits <- traits
  } else {
    fixed_args$trait <- traits
  }

  fixed_args <- c(fixed_args, list(...))

  if (!has_dots) {
    fixed_args <- fixed_args[names(fixed_args) %in% fun_args]
  }

  # x_i = 1 - rho_i
  log_f <- function(x) {

    sampling_fraction <- 1 - x

    result <- do.call(
      fun,
      c(
        list(sampling_fraction = sampling_fraction),
        fixed_args
      )
    )

    loglik <- if (is.list(result)) result$loglik else result

    # Remove the rho_i^S_i factors
    loglik - sum(S * log(sampling_fraction))
  }

  # All vectors k such that 0 <= k_i <= upper_i
  multi_indices <- function(upper) {
    unname(as.matrix(
      expand.grid(lapply(upper, function(x) 0:x))
    ))
  }

  # Numerical mixed partial derivative at x = 0
  mixed_partial <- function(f, orders, step_size = NULL) {

    derivative_order <- sum(orders)

    if (derivative_order == 0) {
      return(f(rep(0, length(orders))))
    }

    if (is.null(step_size)) {
      step_size <- .Machine$double.eps^(
        1 / (derivative_order + 2)
      )
    }

    indices <- multi_indices(orders)

    coefficients <- apply(
      indices,
      1,
      function(index) {
        prod((-1)^index * choose(orders, index))
      }
    )

    points <- sweep(-indices, 2, orders / 2, "+") * step_size
    values <- apply(points, 1, f)
    terms <- coefficients * values

    sum(terms[order(abs(terms))]) / step_size^derivative_order
  }

  # Bell-polynomial recursion:
  # derivatives of f = exp(g), given derivatives of g = log(f)
  bell_polynomials <- function(orders, f_at_zero, log_derivatives) {

    bell <- array(0, dim = orders + 1)

    # B_0 = f(0)
    bell[matrix(rep(1, length(orders)), nrow = 1)] <- f_at_zero

    indices <- multi_indices(orders)

    for (row in seq_len(nrow(indices))) {
      k <- indices[row, ]

      if (all(k == 0)) {
        next
      }

      dimension <- which(k > 0)[1]

      previous <- k
      previous[dimension] <- previous[dimension] - 1

      j_indices <- multi_indices(previous)
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

      bell[matrix(k + 1, nrow = 1)] <-
        sum(terms[order(abs(terms))])
    }

    bell
  }

  # Derivatives of g = log(f) at x = 0
  log_derivatives <- array(0, dim = missnumspec + 1)

  indices <- multi_indices(missnumspec)

  for (row in seq_len(nrow(indices))) {
    orders <- indices[row, ]

    if (!all(orders == 0)) {
      log_derivatives[matrix(orders + 1, nrow = 1)] <-
        suppressWarnings(
          mixed_partial(
            f = log_f,
            orders = orders,
            step_size = step_size
          )
        )
    }
  }

  # f(0), kept on the log scale
  log_f_at_zero <- log_f(rep(0, num_observed_states))

  # Since f(0) is factored out, start the Bell recursion at 1
  bell <- bell_polynomials(
    orders = missnumspec,
    f_at_zero = 1,
    log_derivatives = log_derivatives
  )

  derivative <- bell[matrix(missnumspec + 1, nrow = 1)]

  if (is.na(derivative) || derivative <= 0) {
    warning(
      "The mixed derivative is not positive; ",
      "finite differences may be numerically unstable. Returning -Inf."
    )

    loglikelihood <- -Inf
  } else {
    loglikelihood <-
      log_f_at_zero +
      log(derivative) +
      sum(lfactorial(S)) -
      sum(lfactorial(S + missnumspec))
  }

  list(
    loglik = loglikelihood,
    deriv = exp(log_f_at_zero) * derivative,
    g_derivs = log_derivatives,
    S = S
  )
}

