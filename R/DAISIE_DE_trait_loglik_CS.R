DAISIE_DE_trait_loglik_CS <- function( parameter,
                                       datalist,
                                       methode = "lsodes",
                                       atol = 1e-15,
                                       rtol = 1e-15,
                                       num_observed_states,
                                       num_hidden_states,
                                       get_initial_conditions2 = get_initial_conditions2,
                                       get_initial_conditions3 = get_initial_conditions3,
                                       get_initial_conditions4 = get_initial_conditions4,
                                       func_for_solution       = func_for_solution,
                                       cond =1)

{

  if (length(parameter) == 6) {
    logp0 <- DAISIE_DE_trait_logp0(data_list,
                                   parameter,
                                   atol = 1e-10,
                                   rtol = 1e-10,
                                   num_observed_states,
                                   num_hidden_states,
                                   methode = methode)
    if (is.null(datalist[[1]]$not_present)) {
      loglik <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) * logp0
      numimm <- (datalist[[1]]$not_present_type1 + datalist[[1]]$not_present_type2) + length(datalist) - 1
    } else {
      loglik <- datalist[[1]]$not_present * logp0
      numimm <- datalist[[1]]$not_present + length(datalist) - 1
    }

    ### condition on at least one successful colonization
    logcond <- (cond == 1) * log(1 - exp(numimm * logp0))
    for (i in 2:length(datalist)) {
      datalist[[i]]$type1or2 <- 1
    }
  }

  loglik <- loglik - logcond
  vec_loglikelihood <- c()

  for (i in 2:length(datalist)) {

    stac <- data_list[[i]]$stac
    brts <- data_list[[i]]$branching_times
    traits <- data_list[[i]]$traits
    trait <- data_list[[i]]$traits
    trait_mainland_ancestor <- data_list[[i]]$root_state[2]
    sf <- data_list[[i]]$sampling_fraction[2]
    sf0 <- data_list[[i]]$sampling_fraction[2]
    sf1 <- data_list[[i]]$sampling_fraction[2]
    phy <- data_list[[i]]$phylogeny
    if (stac %in% c(1, 4)) {
      loglikelihood <- DAISIE_DE_trait_logpNE(brts, status = stac, trait, trait_mainland_ancestor = FALSE,
                                              num_observed_states,
                                              num_hidden_states,
                                              parameter,
                                              atol  = 1e-15,
                                              rtol  = 1e-15,
                                              methode                 = methode,
                                              get_initial_conditions2 = get_initial_conditions2,
                                              get_initial_conditions3 = get_initial_conditions3,
                                              get_initial_conditions4 = get_initial_conditions4,
                                              func_for_solution       = func_for_solution)
    } else if (stac %in% c(2, 5)) {
      if (length(brts) == 2)
        loglikelihood <- DAISIE_DE_trait_logpES(brts,
                                                status = stac,
                                                trait,
                                                sf = 1,
                                                trait_mainland_ancestor = FALSE,
                                                num_observed_states,
                                                num_hidden_states,
                                                parameter,
                                                atol  = 1e-10,
                                                rtol  = 1e-10,
                                                methode                 = methode,
                                                get_initial_conditions2 = get_initial_conditions2,
                                                get_initial_conditions4 = get_initial_conditions4,
                                                get_initial_conditions3 = get_initial_conditions3,
                                                func_for_solution       = func_for_solution)
      else
        loglikelihood <- DAISIE_DE_trait_logpEC(brts,
                                                parameter,
                                                phy,
                                                traits = traits,
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
                                                methode = methode)

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
                                                methode = methode,
                                                get_initial_conditions2 = get_initial_conditions2,
                                                get_initial_conditions4 = get_initial_conditions4,
                                                get_initial_conditions3 = get_initial_conditions3,
                                                func_for_solution       = func_for_solution)
      else
        loglikelihood <- DAISIE_DE_trait_logpEC( brts,
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
                                                 methode = methode
        )
    }
    else if (stac == 6) {
      loglikelihood <- DAISIE_DE_trait_logpEC(brts,
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
                                              methode = methode
      )
    } else if (stac == 7) {
      if (length(brts) == 2)
        loglikelihood <- DAISIE_DE_trait_logpES_max_age_coltime_and_mainland_hidden(brts,parameter,methode,rtol,atol)
      else
        loglikelihood <- DAISIE_DE_trait_logpEC_max_age_coltime_and_mainland_hidden(brts,parameter,methode,rtol,atol)
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


load("C:/Users/P306709/Desktop/chapter3/Canary_angiosperms_HW.RData")

DAISIE_DE_trait_loglik_CS (parameter,
                           datalist = data_list[12:20],
                           methode = "lsodes",
                           atol = 1e-10,
                           rtol = 1e-10,
                           num_observed_states = 2,
                           num_hidden_states = 2,
                           get_initial_conditions2 = get_initial_conditions2,
                           get_initial_conditions3 = get_initial_conditions3,
                           get_initial_conditions4 = get_initial_conditions4,
                           func_for_solution       = func_for_solution)

datalist <- data_list[[366]]
i <- 3

parameter <- list(
   c(2.546591, 1.2, 1, 0.2),
   c(2.678781, 2, 1.9, 3),
   c(0.009326754, 0.003, 0.002, 0.2),
   c(1.008583, 1, 2, 1.5),
   matrix(c(
     0,    1,    0.5,  0,
     0,    0,    0.002,0.005,
    rep(0, 8)
   ), nrow = 4),
  0)
