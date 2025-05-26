
DAISIE_DE_trait_loglik_CS <- function( parameter,
                                 datalist,
                                 methode = "lsodes",
                                 abstolint = 1e-15,
                                 reltolint = 1e-15)

{


  if (length(parameter) == 6) {
    logp0 <- DAISIE_DE_trait_logp0(brts,
                                   parameter,
                                   atol = 1e-10,
                                   rtol = 1e-10,
                                   methode = "ode45")
    if (is.null(datalist[[1]]$not_present)) {
      loglik <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0
      numimm <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      loglik <- datalist[[1]]$not_present * logp0
      numimm <- datalist[[1]]$not_present + length(datalist) - 1
    }
    logcond <- (cond == 1) * log(1 - exp(numimm * logp0))
    for (i in 2:length(datalist)) {
      datalist[[i]]$type1or2 <- 1
    }
  } else {
    numimm <- datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2 + length(datalist) - 1
    numimm_type2 <- length(which(unlist(datalist)[which(names(unlist(datalist)) == "type1or2")] == 2))
    numimm_type1 <- length(datalist) - 1 - numimm_type2
    if (!is.na(parameter[11])) {
      if (parameter[11] < numimm_type2 / numimm | parameter[11] > (1 - numimm_type1 / numimm)) {
        return(-Inf)
      }
      datalist[[1]]$not_present_type2 <- max(0, round(parameter[11] * numimm) - numimm_type2)
      datalist[[1]]$not_present_type1 <- numimm - (length(datalist) - 1) - datalist[[1]]$not_present_type2
    }
    logp0_type1 <- DAISIE_DE_trait_logp0(brts,
                                                      parameter,
                                                      atol = 1e-10,
                                                      rtol = 1e-10,
                                                      methode = "ode45")
    logp0_type2 <- DAISIE_DE_trait_logp0(brts,
                                                      parameter,
                                                      atol = 1e-10,
                                                      rtol = 1e-10,
                                                      methode = "ode45")
    loglik <- datalist[[1]]$not_present_type1 * logp0_type1 + datalist[[1]]$not_present_type2 * logp0_type2
    logcond <- (cond == 1) * log(1 - exp((datalist[[1]]$not_present_type1 + numimm_type1) * logp0_type1 +
                                           (datalist[[1]]$not_present_type2 + numimm_type2) * logp0_type2))
  }

  loglik <- loglik - logcond
  vec_loglikelihood <- c()

  for (i in 2:length(datalist)) {

    stac <- datalist[[i]]$stac
    brts <- datalist[[i]]$branching_times
    missnumspec <- datalist[[i]]$missing_species

    if (stac == 1) {
      loglikelihood <- DAISIE_DE_trait_logpNE_hidden(brts, status = 1, trait, trait_mainland_ancestor = FALSE,
                                                             num_observed_states,
                                                             num_hidden_states,
                                                             parameter,
                                                             atol  = 1e-15,
                                                             rtol  = 1e-15,
                                                             methode                 = "ode45",
                                                             get_initial_conditions2 = get_initial_conditions2,
                                                             get_initial_conditions3 = get_initial_conditions3,
                                                             get_initial_conditions4 = get_initial_conditions4,
                                                             func_for_solution       = func_for_solution)
    } else if (stac == 2) {
      if (length(brts) == 2)
        loglikelihood <- DAISIE_DE_trait_logpES(brts,
                                                 status = 2,
                                                 trait,
                                                 sf = 1,
                                                 trait_mainland_ancestor = FALSE,
                                                 num_observed_states,
                                                 num_hidden_states,
                                                 parameter,
                                                 atol  = 1e-10,
                                                 rtol  = 1e-10,
                                                 methode                 = "ode45",
                                                 get_initial_conditions2 = get_initial_conditions2,
                                                 get_initial_conditions4 = get_initial_conditions4,
                                                 get_initial_conditions3 = get_initial_conditions3,
                                                 func_for_solution       = func_for_solution)
      else
        loglikelihood <- DAISIE_DE_trait_logpEC(
                                                            brts,
                                                            parameter,
                                                            phy,
                                                            traits,
                                                            num_observed_states,
                                                            num_hidden_states,
                                                            trait_mainland_ancestor = FALSE,
                                                            status = 2,
                                                            sampling_fraction,
                                                            see_ancestral_states = TRUE,
                                                            atol = 1e-10,
                                                            rtol = 1e-10,
                                                            func_for_solution,
                                                            get_initial_conditions2,
                                                            get_initial_conditions3,
                                                            get_initial_conditions4,
                                                            rhs_func = loglik_hidden_rhs,
                                                            get_func_interval,
                                                            methode = "ode45"
                                                                )
    } else if (stac == 3) {
      if (length(brts) == 2)
        loglikelihood <- DAISIE_DE_trait_logpES(brts,
                                                status = 3,
                                                trait,
                                                sf = 1,
                                                trait_mainland_ancestor = FALSE,
                                                num_observed_states,
                                                num_hidden_states,
                                                parameter,
                                                atol  = 1e-10,
                                                rtol  = 1e-10,
                                                methode                 = "ode45",
                                                get_initial_conditions2 = get_initial_conditions2,
                                                get_initial_conditions4 = get_initial_conditions4,
                                                get_initial_conditions3 = get_initial_conditions3,
                                                func_for_solution       = func_for_solution)
      else
        loglikelihood <- DAISIE_DE_trait_logpEC(
                                                  brts,
                                                  parameter,
                                                  phy,
                                                  traits,
                                                  num_observed_states,
                                                  num_hidden_states,
                                                  trait_mainland_ancestor = FALSE,
                                                  status = 3,
                                                  sampling_fraction,
                                                  see_ancestral_states = TRUE,
                                                  atol = 1e-10,
                                                  rtol = 1e-10,
                                                  func_for_solution,
                                                  get_initial_conditions2,
                                                  get_initial_conditions3,
                                                  get_initial_conditions4,
                                                  rhs_func = loglik_hidden_rhs,
                                                  get_func_interval,
                                                  methode = "ode45"
                                                )
    } else if (stac == 4) {
      loglikelihood <- DAISIE_DE_trait_logpNE_hidden(brts, status = 4, trait, trait_mainland_ancestor = FALSE,
                                              num_observed_states,
                                              num_hidden_states,
                                              parameter,
                                              atol  = 1e-15,
                                              rtol  = 1e-15,
                                              methode                 = "ode45",
                                              get_initial_conditions2 = get_initial_conditions2,
                                              get_initial_conditions3 = get_initial_conditions3,
                                              get_initial_conditions4 = get_initial_conditions4,
                                              func_for_solution       = func_for_solution)
    } else if (stac == 5) {
      loglikelihood <- DAISIE_DE_trait_logpES(brts,
                                                          status,
                                                          trait,
                                                          sf = 5,
                                                          trait_mainland_ancestor = FALSE,
                                                          num_observed_states,
                                                          num_hidden_states,
                                                          parameter,
                                                          atol  = 1e-10,
                                                          rtol  = 1e-10,
                                                          methode                 = "ode45",
                                                          get_initial_conditions2 = get_initial_conditions2,
                                                          get_initial_conditions4 = get_initial_conditions4,
                                                          get_initial_conditions3 = get_initial_conditions3,
                                                          func_for_solution       = func_for_solution)
    } else if (stac == 6) {
      loglikelihood <- DAISIE_DE_trait_logpEC(
        brts,
        parameter,
        phy,
        traits,
        num_observed_states,
        num_hidden_states,
        trait_mainland_ancestor = FALSE,
        status = 6,
        sampling_fraction,
        see_ancestral_states = TRUE,
        atol = 1e-10,
        rtol = 1e-10,
        func_for_solution,
        get_initial_conditions2,
        get_initial_conditions3,
        get_initial_conditions4,
        rhs_func = loglik_hidden_rhs,
        get_func_interval,
        methode = "ode45"
      )
    } else if (stac == 7) {
      if (length(brts) == 2)
        loglikelihood <- DAISIE_DE_trait_logpES_max_age_coltime_and_mainland_hidden(brts,parameter,methode,reltolint,abstolint)
      else
        loglikelihood <- DAISIE_DE_trait_logpEC_max_age_coltime_and_mainland_hidden(brts,parameter,methode,reltolint,abstolint)
    } else if (stac == 8) {
      loglikelihood <- DAISIE_DE_trait_logpNE_max_min_age_coltime_hidden(brts, trait, parameter, num_observed_states, num_hidden_states, atol = 1e-10, rtol = 1e-10, methode = "ode45")
    } else if (stac == 9) {
      loglikelihood <- DAISIE_DE_trait_logpES_max_min_age_coltime_hidden(brts, trait, parameter, num_observed_states, num_hidden_states, atol = 1e-10, rtol = 1e-10, methode = "ode45")
    } else {
      stop("Unknown stac value: ", stac)
    }

    vec_loglikelihood <- c(vec_loglikelihood, loglikelihood)


  }
  loglik <- sum(vec_loglikelihood) + loglik
  return(loglik)
}
