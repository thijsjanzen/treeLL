
library(DAISIE)
library(ape)
library(secsse)

set.seed(1)



parameter <- list( c(2.546591, 2.546591), c(2.678781, 2.678781), c(0.009326754, 0.009326754),
                   c(1.008583, 1.008583), matrix(rep(0, 4), nrow = 2), 0)

data("Galapagos_datalist")
datalist <- Galapagos_datalist
i <- 4
phy <- DDD::brts2phylo(datalist[[i]]$branching_times[-c(1, 2)])
traits <- sample(c(0, 1),ape::Ntip(phy),replace=TRUE)
num_hidden_traits <- 1


DAISIE_DE_logpEC_trait_loglik <- DAISIE_DE_logpEC_trait1_hidden(brts = datalist[[i]]$branching_times,
                                                                missnumspec = datalist[[i]]$missing_species,
                                                                parameter = parameter,
                                                                phy,
                                                                traits,
                                                                num_hidden_traits,
                                                                cond = "proper_cond",
                                                                root_state_weight = "proper_weights",
                                                                see_ancestral_states = TRUE,
                                                                atol = 1e-10,
                                                                rtol = 1e-10,
                                                                methode = "ode45",
                                                                rhs_func = loglik_hidden_rhs)
DAISIE_DE_logpEC_trait_loglik
