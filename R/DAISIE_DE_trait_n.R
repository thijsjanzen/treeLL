#' @name DAISIE_DE_trait_n
#' @title Convert a trait-dependent rho-sampling likelihood into an n-sampling
#' likelihood
#' @description The trait-dependent DAISIE_DE likelihood functions
#' (\code{DAISIE_DE_trait_logpES}, \code{DAISIE_DE_trait_logpEC}, ...) are
#' parameterised by a vector of sampling probabilities
#' \code{sampling_fraction} \eqn{= (\rho_1, \dots, \rho_d)}, one per observed
#' trait state, i.e. they perform rho-sampling. This function converts such a
#' likelihood into the n-sampling likelihood, in which the numbers of unsampled
#' species per observed trait state \eqn{(n_1, \dots, n_d)} are specified
#' instead.
#'
#' With \eqn{S_i} the number of sampled species in observed state \eqn{i}, the
#' two are related by
#' \deqn{P_{\rho_1,\rho_2}(clade) = \sum_{n_1}\sum_{n_2}
#'       P_{n_1,n_2}(clade)
#'       \binom{S_1+n_1}{n_1}\rho_1^{S_1}(1-\rho_1)^{n_1}
#'       \binom{S_2+n_2}{n_2}\rho_2^{S_2}(1-\rho_2)^{n_2}.}
#' Writing \eqn{x_i = 1 - \rho_i} and
#' \deqn{f(x_1,x_2) = \frac{P_{1-x_1,1-x_2}(clade)}
#'                         {(1-x_1)^{S_1}(1-x_2)^{S_2}},}
#' the right-hand side is the two-dimensional Taylor expansion of \eqn{f} about
#' \eqn{(0,0)}, so that
#' \deqn{P_{n_1,n_2}(clade) = \frac{S_1! S_2!}{(S_1+n_1)! (S_2+n_2)!}
#'       \left.\frac{\partial^{n_1+n_2} f(x_1,x_2)}
#'                  {\partial x_1^{n_1}\partial x_2^{n_2}}\right|_{(0,0)}.}
#' The implementation below is written for an arbitrary number \eqn{d} of
#' observed trait states; \eqn{d = 2} reproduces the equations above and
#' \eqn{d = 1} reproduces the trait-free \code{DAISIE_DE_n}.
#'
#' As in the trait-free case the mixed partial derivatives are obtained from
#' the derivatives of \eqn{g = \log f} rather than from \eqn{f} itself, because
#' it is \eqn{\log P} that the likelihood functions return. In one dimension
#' \eqn{f^{(n)} = f B_n(g', \dots, g^{(n)})} with \eqn{B_n} the Bell
#' polynomials; the multivariate analogue follows from applying Leibniz's rule
#' to \eqn{\partial_r f = f \partial_r g},
#' \deqn{\partial^k f = \sum_{j \le k - e_r} \binom{k - e_r}{j}
#'       (\partial^{k - e_r - j} f)(\partial^{j + e_r} g),}
#' where \eqn{k} is a multi-index, \eqn{r} is any direction with \eqn{k_r > 0}
#' and \eqn{e_r} is the corresponding unit multi-index. For \eqn{d = 1} this is
#' exactly the recursion \eqn{B_n = \sum_k \binom{n-1}{k-1} B_{n-k} g^{(k)}}.
#'
#' @param DAISIE_DE_trait_function The rho-sampling likelihood function to
#' convert, e.g. \code{DAISIE_DE_trait_logpES} or
#' \code{DAISIE_DE_trait_logpEC}. It must take a \code{sampling_fraction}
#' argument and return either a numeric loglikelihood or a list with a
#' \code{loglik} element.
#' @param brts Branching times of the clade, as passed to
#' \code{DAISIE_DE_trait_function}.
#' @param missnumspec Numbers of unsampled (missing) species per observed trait
#' state, i.e. \eqn{(n_1, \dots, n_d)}. Must have length
#' \code{num_observed_states}.
#' @param traits Trait states of the sampled species, coded
#' \code{0, ..., num_observed_states - 1} (a single value for a singleton
#' lineage). Passed on to \code{DAISIE_DE_trait_function} as \code{traits} or
#' \code{trait}, whichever that function accepts.
#' @param status The stac / status of the colonist.
#' @param parameter The model parameters, as a list.
#' @param num_observed_states Number of observed trait states.
#' @param num_hidden_states Number of hidden trait states.
#' @param S Numbers of sampled species per observed trait state,
#' \eqn{(S_1, \dots, S_d)}. Defaults to the counts derived from \code{traits}.
#' @param h Step size of the finite differences. The default,
#' \code{.Machine$double.eps^(1 / (sum(missnumspec) + 2))}, balances truncation
#' against round-off error for a central difference of total order
#' \code{sum(missnumspec)} and matches \code{pracma::fderiv} in one dimension.
#' @param ... Further arguments passed on to \code{DAISIE_DE_trait_function},
#' e.g. \code{datalist}, \code{phy}, \code{trait_mainland_ancestor},
#' \code{atol}, \code{rtol}, \code{methode}, \code{rcpp_methode},
#' \code{use_Rcpp}, \code{num_threads}. Arguments the function does not accept
#' are silently dropped.
#' @return A list with elements \code{loglik} (the n-sampling loglikelihood),
#' \code{deriv} (the mixed partial derivative of \eqn{f}), \code{g_derivs} (the
#' array of mixed partial derivatives of \eqn{g = \log f}) and \code{S}.
#' @author Rampal S. Etienne & Bart Haegeman
#' @seealso \code{\link{DAISIE_DE_trait_logpES}},
#' \code{\link{DAISIE_DE_trait_logpEC}}
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
                              h = NULL,
                              ...) {

  d <- num_observed_states

  if (length(missnumspec) != d) {
    stop("missnumspec must give the number of missing species per observed ",
         "trait state, so it must have length num_observed_states (", d, ").")
  }
  missnumspec <- as.integer(round(missnumspec))
  if (any(missnumspec < 0)) {
    stop("missnumspec cannot be negative.")
  }

  # The rho-sampling likelihood is discontinuous at rho = 1 when the trait
  # state of a species is unknown: the is.na(trait) branches of
  # get_initial_conditions2/3 and calc_init_state_hidden zero the DE entries of
  # the fully sampled states as soon as any sampling fraction drops below 1, so
  # the Taylor expansion around rho = 1 is meaningless in that case. With no
  # missing species nothing is expanded, though - the answer is the likelihood
  # at rho = 1 itself, which is perfectly well defined - so the restriction
  # only bites once species are actually missing.
  if (anyNA(traits) && sum(missnumspec) > 0) {
    stop("DAISIE_DE_trait_n requires known trait states when there are ",
         "missing species: the rho-sampling likelihood is not continuous at ",
         "rho = 1 when traits are NA. This clade has ", sum(missnumspec),
         " missing species and ", sum(is.na(traits)), " sampled species of ",
         "unknown trait state.")
  }

  # given S_i, the observed number of species in observed trait state i.
  # Traits are coded 0, ..., d - 1, as in sampling_fraction[1 + trait].
  if (is.null(S)) {
    S <- tabulate(traits + 1, nbins = d)
  }
  if (length(S) != d) {
    stop("S must have length num_observed_states (", d, ").")
  }

  # Assemble the arguments of the rho-sampling likelihood, keeping only those
  # it actually accepts: logpEC takes 'traits' (and 'phy', 'num_threads'),
  # while logpES/logpNE take a single 'trait'.
  fun <- match.fun(DAISIE_DE_trait_function)
  fun_args <- names(formals(fun))
  has_dots <- "..." %in% fun_args

  fixed_args <- list(brts = brts,
                     status = status,
                     parameter = parameter,
                     num_observed_states = num_observed_states,
                     num_hidden_states = num_hidden_states)
  if ("traits" %in% fun_args) {
    fixed_args$traits <- traits
  } else if ("trait" %in% fun_args) {
    fixed_args$trait <- traits
  } else if (has_dots) {
    fixed_args$traits <- traits
  }
  fixed_args <- c(fixed_args, list(...))
  if (!has_dots) {
    fixed_args <- fixed_args[names(fixed_args) %in% fun_args]
  }

  # f(x) and g(x) = log f(x), with x = 1 - sampling_fraction. The trait
  # likelihoods include the rho_i^{S_i} tip factors (DE_i = rho_i at every
  # tip), so these have to be divided out again to obtain f.
  log_f <- function(x) {
    sampling_fraction <- 1 - x
    loglik <- do.call(fun, c(list(sampling_fraction = sampling_fraction),
                             fixed_args))
    if (is.list(loglik)) loglik <- loglik$loglik
    return(loglik - sum(S * log(sampling_fraction)))
  }

  # Every evaluation of log_f is a full ODE solve and neighbouring stencils
  # share nodes, so memoise on the grid point.
  cache <- new.env(hash = TRUE, parent = emptyenv())
  log_f_cached <- function(x) {
    key <- paste(sprintf("%.17g", x), collapse = "|")
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    value <- log_f(x)
    assign(key, value, envir = cache)
    return(value)
  }

  # All multi-indices 0 <= k <= upper, in an order in which every k' <= k
  # (componentwise) comes before k.
  multi_indices <- function(upper) {
    unname(as.matrix(expand.grid(lapply(upper, function(m) 0:m))))
  }

  # Mixed partial derivative of g at 0, as a tensor product of the central
  # difference formula used by pracma::fderiv:
  #   d^k g(0) = h^-|k| sum_{j <= k} (-1)^|j| choose(k, j) g((k/2 - j) h)
  mixed_partial <- function(g, k, h = NULL) {
    order_k <- sum(k)
    if (order_k == 0) return(g(rep(0, length(k))))
    if (is.null(h)) h <- .Machine$double.eps^(1 / (order_k + 2))
    js <- multi_indices(k)
    coefs <- apply(js, 1, function(j) prod((-1)^j * choose(k, j)))
    nodes <- sweep(-js, 2, k / 2, "+") * h
    values <- apply(nodes, 1, g)
    terms <- coefs * values
    o <- order(abs(terms))
    return(sum(terms[o]) / h^order_k)
  }

  # Forward-difference alternative, the analogue of nth_deriv_richardson
  # mixed_partial_forward <- function(g, k, h = 1e-5) {
  #   js <- multi_indices(k)
  #   coefs <- apply(js, 1, function(j) prod((-1)^j * choose(k, j)))
  #   nodes <- sweep(-js, 2, k, "+") * h
  #   values <- apply(nodes, 1, g)
  #   return(sum(coefs * values) / h^sum(k))  # simply inaccurate
  # }

  # Multivariate Bell polynomials: the array of mixed partial derivatives of
  # f = exp(g) at 0, given f(0) and the mixed partial derivatives of g at 0.
  bell_polynomials_up_to_n <- function(nvec, f_val, g_derivs) {
    B <- array(0, dim = nvec + 1)
    B[matrix(rep(1, length(nvec)), nrow = 1)] <- f_val  # B_0
    ks <- multi_indices(nvec)

    for (row in seq_len(nrow(ks))) {
      k <- ks[row, ]
      if (all(k == 0)) next
      r <- which(k > 0)[1]
      km <- k
      km[r] <- km[r] - 1                                 # k - e_r

      js <- multi_indices(km)
      tmp <- numeric(nrow(js))
      for (p in seq_len(nrow(js))) {
        j <- js[p, ]
        kg <- j
        kg[r] <- kg[r] + 1                               # j + e_r
        tmp[p] <- prod(choose(km, j)) *
          B[matrix(km - j + 1, nrow = 1)] *
          g_derivs[matrix(kg + 1, nrow = 1)]
      }
      o <- order(abs(tmp))
      B[matrix(k + 1, nrow = 1)] <- sum(tmp[o])
    }
    return(B)  # B[1, ..., 1] = B_0, ..., B[nvec + 1] = B_nvec
  }

  nth_derivative_from_log <- function(nvec, f_val, g_derivs) {
    B <- bell_polynomials_up_to_n(nvec, f_val, g_derivs)
    return(B[matrix(nvec + 1, nrow = 1)])
  }

  # Multivariate version of Cauchy's integral formula, the analogue of the
  # one-dimensional contour integral. Needs f for complex sampling fractions.
  # f <- function(x) exp(log_f(x))
  # integrand <- function(t) {
  #   z <- exp(1i * t)
  #   fz <- 1 / (2 * pi * 1i)^d * f(z) / prod(z^(missnumspec + 1))
  #   dz_dt <- prod(1i * z)
  #   return(Re(fz * dz_dt))
  # }

  lderiv <- array(0, dim = missnumspec + 1)
  ks <- multi_indices(missnumspec)
  for (row in seq_len(nrow(ks))) {
    k <- ks[row, ]
    if (all(k == 0)) next
    lderiv[matrix(k + 1, nrow = 1)] <-
      suppressWarnings(mixed_partial(log_f_cached, k = k, h = h))
  }

  # f(0) is factored out in log space instead of being passed in as f_val, so
  # that a very small likelihood cannot underflow before the Bell recursion.
  g0 <- log_f_cached(rep(0, d))
  deriv <- nth_derivative_from_log(nvec = missnumspec,
                                   f_val = 1,
                                   g_derivs = lderiv)

  if (is.na(deriv) || deriv <= 0) {
    warning("The mixed partial derivative of f is not positive (", deriv,
            "); the finite differences are probably dominated by round-off. ",
            "Returning -Inf.")
    loglikelihood <- -Inf
  } else {
    loglikelihood <- g0 + log(deriv) +
      sum(lfactorial(S)) - sum(lfactorial(S + missnumspec))
  }

  return(list(loglik = loglikelihood,
              deriv = exp(g0) * deriv,
              g_derivs = lderiv,
              S = S))
}
