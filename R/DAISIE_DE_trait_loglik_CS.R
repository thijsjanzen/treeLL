DAISIE_DE_trait_loglik_CS <- function( parameter,
                                       datalist,
                                       methode = "lsodes",
                                       rcpp_methode =
                                         "odeint::runge_kutta_cash_karp54",
                                       atol = 1e-15,
                                       rtol = 1e-15,
                                       num_observed_states,
                                       num_hidden_states,
                                       cond = 1,
                                       num_threads = 1,
                                       verbose = FALSE,
                                       use_Rcpp = use_Rcpp,
                                       sampling = "rho")

{
  # "rho" uses the sampling fractions stored in the datalist, "n" the numbers
  # of unsampled species per observed trait state, taken from
  # datalist[[i]]$missing_species
  sampling <- match.arg(sampling, c("rho", "n"))

  logcond <- 0 # default value gives no effect

  if (length(parameter) >= 6) {


    logp0 <- DAISIE_DE_trait_logp0(datalist = datalist,
                                   parameter = parameter,
                                   atol = atol,
                                   rtol = rtol,
                                   num_observed_states = num_observed_states,
                                   num_hidden_states = num_hidden_states,
                                   trait_mainland_ancestor =  NA,
                                   methode = methode,
                                   use_Rcpp = use_Rcpp)

    if (is.null(datalist[[1]]$not_present)) {
      loglik <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0$loglik
      numimm <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      not_present_by_state <- datalist[[1]]$not_present_by_state
      not_present_NA <- datalist[[1]]$not_present_NA

      # Empirical distribution across OBSERVED states
      p <- not_present_by_state / sum(not_present_by_state)

      # Add the NA species according to the observed-state proportions
      effective_counts_obs <-
        not_present_by_state + not_present_NA * p

      # Expand each observed state over its hidden states
      # assuming hidden states are equally likely
      effective_counts_full <- rep(
        effective_counts_obs / num_hidden_states,
        each = num_hidden_states
      )

      # logp0 has length:
      # num_observed_states * num_hidden_states
      loglik <- sum(effective_counts_full * logp0)

      numimm <-
        sum(not_present_by_state) +
        not_present_NA +
        length(datalist) - 1
    }

    ### pool p0 over the mainland state distribution -> scalar
    w <- effective_counts_full / sum(effective_counts_full)
    logp0_pooled <- log(sum(w * exp(logp0)))

    ### condition on at least one successful colonization
    logcond <- (cond == 1) * log(1 - exp(numimm * logp0_pooled))
    for (i in 2:length(datalist)) {
      datalist[[i]]$type1or2 <- 1
    }
  }

  loglik <- loglik - logcond

  vec_loglikelihood <- rep(NA, length(datalist) - 1) # first entry is not data
  for (i in 2:length(datalist)) {
    stac <- datalist[[i]]$stac
    brts <- datalist[[i]]$branching_times
    traits <- datalist[[i]]$traits
    trait <- datalist[[i]]$traits

    sampling_fraction <- datalist[[i]]$sampling_fraction

    phy <- datalist[[i]]$phylogeny

    trait_mainland_ancestor <- datalist[[i]]$root_state

    # Select the likelihood function for this colonist. stac 8 and 9 are
    # integrated with ode45 regardless of the requested method.
    clade_methode <- methode
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

    # The likelihood functions differ in which of these they accept (logpEC
    # takes 'traits' and 'phy', logpES and logpNE a single 'trait'), so only
    # the arguments a function declares are passed on to it.
    clade_args <- list(datalist                = datalist,
                       brts                    = brts,
                       parameter               = parameter,
                       phy                     = phy,
                       traits                  = traits,
                       trait                   = trait,
                       num_observed_states     = num_observed_states,
                       num_hidden_states       = num_hidden_states,
                       trait_mainland_ancestor = trait_mainland_ancestor,
                       status                  = stac,
                       atol                    = atol,
                       rtol                    = rtol,
                       methode                 = clade_methode,
                       rcpp_methode            = rcpp_methode,
                       use_Rcpp                = use_Rcpp,
                       num_threads             = num_threads)

    if (sampling == "rho") {
      clade_args$sampling_fraction <- sampling_fraction
      loglikelihood <- do.call(
        logp_function,
        clade_args[names(clade_args) %in% names(formals(logp_function))])
    } else {
      # n-sampling: the numbers of unsampled species per observed trait state
      # come straight from the datalist and replace the sampling fractions.
      # DAISIE_DE_trait_n differentiates the rho-sampling likelihood to get
      # there, and supplies 'traits' itself under whichever name
      # logp_function uses - hence dropping the singular 'trait' first, so it
      # is not passed twice.
      loglikelihood <- do.call(
        DAISIE_DE_trait_n,
        c(list(DAISIE_DE_trait_function = logp_function,
               missnumspec = datalist[[i]]$missing_species),
          clade_args[names(clade_args) != "trait"]))
    }

    vec_loglikelihood[i - 1] <- loglikelihood$loglik
  }

  loglik <- sum(vec_loglikelihood) + loglik
  return(loglik)

}
