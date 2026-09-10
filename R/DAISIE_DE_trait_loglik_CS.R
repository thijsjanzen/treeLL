make_log_f_from_clade <- function(logp_function, clade_args, S) {

  accepted_arguments <- names(formals(logp_function))

  clade_args <- clade_args[
    names(clade_args) %in% accepted_arguments
  ]

  function(x) {
    rho <- 1 - x

    result <- do.call(
      logp_function,
      c(
        list(sampling_fraction = rho),
        clade_args
      )
    )

    loglik <- if (is.list(result)) result$loglik else result

    loglik - sum(S * log(rho))
  }
}


DAISIE_DE_trait_loglik_CS <- function(
    parameter,
    datalist,
    methode = "lsodes",
    rcpp_methode = "odeint::runge_kutta_cash_karp54",
    atol = 1e-15,
    rtol = 1e-15,
    num_observed_states,
    num_hidden_states,
    cond = 1,
    num_threads = 1,
    verbose = FALSE,
    use_Rcpp = 2,
    sampling = "rho"
) {

  sampling <- match.arg(sampling, c("rho", "n"))

  loglik <- 0
  logcond <- 0

  # Contribution of mainland species that are not present on the island
  if (length(parameter) >= 6) {

    logp0 <- DAISIE_DE_trait_logp0(
      datalist = datalist,
      parameter = parameter,
      atol = atol,
      rtol = rtol,
      num_observed_states = num_observed_states,
      num_hidden_states = num_hidden_states,
      trait_mainland_ancestor = NA,
      methode = methode,
      use_Rcpp = use_Rcpp
    )$loglik

    if (is.null(datalist[[1]]$not_present)) {

      not_present <-
        datalist[[1]]$not_present_type1 +
        datalist[[1]]$not_present_type2

      loglik <- not_present * logp0

      numimm <- not_present + length(datalist) - 1

      logp0_pooled <- logp0

    } else {

      not_present_by_state <- datalist[[1]]$not_present_by_state
      not_present_NA <- datalist[[1]]$not_present_NA

      state_proportions <-
        not_present_by_state / sum(not_present_by_state)

      effective_counts_observed <-
        not_present_by_state +
        not_present_NA * state_proportions

      effective_counts_full <- rep(
        effective_counts_observed / num_hidden_states,
        each = num_hidden_states
      )

      loglik <- sum(effective_counts_full * logp0)

      numimm <-
        sum(not_present_by_state) +
        not_present_NA +
        length(datalist) - 1

      weights <- effective_counts_full / sum(effective_counts_full)

      logp0_pooled <- log(sum(weights * exp(logp0)))
    }

    if (cond == 1) {
      logcond <- log(1 - exp(numimm * logp0_pooled))
    }

    for (i in 2:length(datalist)) {
      datalist[[i]]$type1or2 <- 1
    }
  }

  loglik <- loglik - logcond

  clade_loglikelihoods <- numeric(length(datalist) - 1)

  # Compute the likelihood of every colonist lineage
  for (i in 2:length(datalist)) {

    stac <- datalist[[i]]$stac
    brts <- datalist[[i]]$branching_times
    traits <- datalist[[i]]$traits
    phy <- datalist[[i]]$phylogeny
    root_state <- datalist[[i]]$root_state

    clade_methode <- methode

    # Select the appropriate rho-sampling likelihood
    if (stac %in% c(1, 4)) {

      logp_function <- DAISIE_DE_trait_logpNE

    } else if (stac %in% c(2, 3, 5)) {

      if (length(brts) == 2) {
        logp_function <- DAISIE_DE_trait_logpES
      } else {
        logp_function <- DAISIE_DE_trait_logpEC
      }

    } else if (stac == 6) {

      logp_function <- DAISIE_DE_trait_logpEC

    } else if (stac == 8) {

      logp_function <- DAISIE_DE_trait_logpNE_max_min_age_hidden
      clade_methode <- "ode45"

    } else if (stac == 9) {

      logp_function <- DAISIE_DE_trait_logpES_max_min_age_hidden
      clade_methode <- "ode45"

    } else {
      stop("Unknown stac value: ", stac)
    }

    clade_args <- list(
      datalist = datalist,
      brts = brts,
      parameter = parameter,
      phy = phy,
      traits = traits,
      trait = traits,
      num_observed_states = num_observed_states,
      num_hidden_states = num_hidden_states,
      trait_mainland_ancestor = root_state,
      status = stac,
      atol = atol,
      rtol = rtol,
      methode = clade_methode,
      rcpp_methode = rcpp_methode,
      use_Rcpp = use_Rcpp,
      num_threads = num_threads
    )

    if (sampling == "rho") {

      clade_args$sampling_fraction <- datalist[[i]]$sampling_fraction

      accepted_arguments <- names(formals(logp_function))

      clade_loglikelihood <- do.call(
        logp_function,
        clade_args[names(clade_args) %in% accepted_arguments]
      )

    } else {

      missnumspec <- datalist[[i]]$missing_species

      if (anyNA(traits) && sum(missnumspec) > 0) {
        stop(
          "Known trait states are required when species are missing."
        )
      }

      S <- count_species_per_trait(
        traits = traits,
        num_observed_states = num_observed_states
      )

      log_f <- make_log_f_from_clade(
        logp_function = logp_function,
        clade_args = clade_args,
        S = S
      )

      clade_loglikelihood <- DAISIE_DE_trait_n(
        log_f = log_f,
        missnumspec = missnumspec,
        S = S
      )
    }

    clade_loglikelihoods[i - 1] <- clade_loglikelihood$loglik
  }

  sum(clade_loglikelihoods) + loglik
}
