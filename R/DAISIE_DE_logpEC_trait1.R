### function used
calcThruNodes <- function(
    ances,
    states,
    loglik,
    forTime,
    parameter,
    methode,
    phy,
    rhs_func,
    rtol,
    atol
) {

  nb_node <- phy$Nnode
  focal <- ances
  desRows <- which(phy$edge[, 1] == focal)
  desNodes <- phy$edge[desRows, 2]

  nodeM <- numeric()
  nodeN <- numeric()

  for (desIndex in 1:2) {
    y <- states[desNodes[desIndex],]
    #
    timeInte <- forTime[which(forTime[, 2] == desNodes[desIndex]), 3]
    ##  To do the calculation in both lineages

    nodeMN <- deSolve::ode(y = y,
                           func = rhs_func,
                           times = c(0, timeInte),
                           parms = parameter,
                           rtol = rtol,
                           atol = atol,
                           method = methode)

    if (desIndex == 1) {
      nodeN <- nodeMN
    }
    if (desIndex == 2) {
      nodeM <- nodeMN
    }
  }
  ## At the node
  nodeM <- as.numeric(nodeM[2,-1])
  nodeN <- as.numeric(nodeN[2,-1])

  combined_state <- nodeN
  lambdac <- parameter[[1]]
  combined_state[1] <- lambdac[1] * nodeN[1] * nodeM[1];
  combined_state[2] <- lambdac[2] * nodeN[2] * nodeM[2];

  states[focal,] <- combined_state
  return(list(states = states,
              loglik = loglik,
              combined_state = combined_state,
              nodeM = nodeM,
              nodeN = nodeN))
}



#' @keywords internal
master_loglik_R <- function(parameter,
                            phy,
                            traits,
                            cond = "proper_cond",
                            root_state_weight = "proper_weights",
                            setting_calculation = NULL,
                            see_ancestral_states = TRUE,
                            atol = 1e-10,
                            rtol = 1e-10,
                            methode = "ode45",
                            use_normalization = TRUE) {

  number_of_lineages <- length(phy$tip.label)

  states <- matrix(nrow = number_of_lineages + phy$Nnode,
                   ncol = 7,
                   data = NA)

  for (i in 1:length(traits)) {
    states[i, ] <- calc_init_state(traits[i])
  }
  colnames(states) <- c("DE_0", "DE_1", "DM3_0", "DM3_1", "E_0", "E_1", "DA3")
  phy$node.label <- NULL
  split_times <- sort(event_times(phy), decreasing = FALSE)
  ances <- as.numeric(names(split_times))
  forTime <- cbind(phy$edge, phy$edge.length)

  d <- ncol(states) / 2
  loglik <- 0

  for(i in 1:length(ances)){
    calcul <- calcThruNodes(ances = ances[i],
                            states = states,
                            loglik = loglik,
                            forTime = forTime,
                            parameter = parameter,
                            methode = methode,
                            phy = phy,
                            rhs_func = loglik_rhs,
                            rtol = rtol,
                            atol = atol)
    states <- calcul$states
    loglik <- calcul$loglik
    nodeN <- calcul$nodeN
  }

  prob_states <- calcul$combined_state
  prob_states <- matrix(prob_states, nrow = 1,
                        dimnames = list(NULL, c("DE_0", "DE_1", "DM3_0", "DM3_1", "E_0", "E_1", "DA3")))


  return(prob_states)
}


DAISIE_DE_logpEC_trait1 <- function (brts,
                                     missnumspec,
                                     parameter,
                                     phy,
                                     traits,
                                     cond = "proper_cond",
                                     root_state_weight = "proper_weights",
                                     setting_calculation = NULL,
                                     see_ancestral_states = TRUE,
                                     atol = 1e-10,
                                     rtol = 1e-10,
                                     methode = "ode45",
                                     use_normalization = TRUE)

{




  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  ti <- sort(brts)
  ti <- ti[1:(length(ti)-2)]

  #########Ininitial conditions [tp, t2]
  calc_init_state <- function(trait) {

    out <- rep(7, 0)
    if (trait == 0) {
      out <- c(DE_0 = 1, DE_1 = 0, DM3_0 = 0, DM3_1 = 0, E_0 = 0, E_1 = 0, DA3 = 1)
    }
    if (trait == 1) {
      out <- c(DE_0 = 0, DE_1 = 1, DM3_0 = 0, DM3_1 = 0, E_0 = 0, E_1 = 0, DA3 = 1)
    }
    return(out)
  }

  #########Interval1 [t_p, t_2]

  loglik_rhs <- function(t, state, parameter) {
    with(as.list(c(state, parameter)), {

      lambdac_0 <- parameter[[1]][1]; lambdac_1 <- parameter[[1]][2]
      mu_0 <- parameter[[2]][1]; mu_1 <- parameter[[2]][2]
      gamma_0 <- parameter[[3]][1]; gamma_1 <- parameter[[3]][2]
      lambdaa_0 <- parameter[[4]][1]; lambdaa_1 <- parameter[[4]][2]
      q_01 <- parameter[[5]][1]; q_10 <- parameter[[5]][2]
      p <- parameter[[6]]

      dDE_0 <- -(lambdac_0 + mu_0 + q_01) * DE_0 + 2 * lambdac_0 * DE_0 * E_0 + q_01 * DE_1

      dDE_1 <- -(lambdac_1 + mu_1 + q_10) * DE_1 + 2 * lambdac_1 * DE_1 * E_1 + q_10 * DE_0


      dDM3_0 <- -(lambdac_0 + mu_0 + gamma_1 + lambdaa_0 + q_01) * DM3_0 +
                 (mu_0 + lambdaa_0 * E_0 + lambdac_0 * E_0^2 + p * q_01 * E_1) * DA3 +
                 (1 - p) * q_01 * DM3_1 + gamma_1 * DM3_1

      dDM3_1 <- -(lambdac_1 + mu_1 + gamma_0 + lambdaa_1 + q_10) * DM3_1 +
                 (mu_1 + lambdaa_1 * E_1 + lambdac_1 * E_1^2 + p * q_10 * E_0) * DA3 +
                 (1 - p) * q_10 * DM3_0 + gamma_0 * DM3_0



      dE_0 <- mu_0 - (mu_0 + lambdac_0 + q_01) * E_0 + lambdac_0 * E_0^2 + q_01 * E_1
      dE_1 <- mu_1 - (mu_1 + lambdac_1 + q_10) * E_1 + lambdac_1 * E_1^2 + q_10 * E_0

      dDA3 <- - (gamma_0 + gamma_1) * DA3 + gamma_0 * DM3_0 + gamma_1 * DM3_1

      list(c(dDE_0, dDE_1, dDM3_0, dDM3_1, dE_0, dE_1, dDA3))
    })
  }



  #########Interval2 [t_2, t_1]

  interval2 <- function(t, state, parameter) {
    with(as.list(c(state, parameter)), {

      lambdac_0 <- parameter[[1]][1]; lambdac_1 <- parameter[[1]][2]
      mu_0 <- parameter[[2]][1]; mu_1 <- parameter[[2]][2]
      gamma_0 <- parameter[[3]][1]; gamma_1 <- parameter[[3]][2]
      lambdaa_0 <- parameter[[4]][1]; lambdaa_1 <- parameter[[4]][2]
      q_01 <- parameter[[5]][1]; q_10 <- parameter[[5]][2]
      p <- parameter[[6]]

      dDE_0  <- -(lambdac_0 + mu_0 + q_01) * DE_0 + 2 * lambdac_0 * DE_0 * E_0 + q_01 * DE_1

      dDE_1  <- -(lambdac_1 + mu_1 + q_10) * DE_1 + 2 * lambdac_1 * DE_1 * E_1 + q_10 * DE_0

      dDM3_0 <- -(lambdac_0 + mu_0 + gamma_1 + lambdaa_0 + q_01) * DM3_0 +
                 (mu_0 + lambdaa_0 * E_0 + lambdac_0 * E_0^2 + p * q_01 * E_1) * DA3 +
                 (1 - p) * q_01 * DM3_1 + gamma_1 * DM3_1

      dDM3_1 <- -(lambdac_1 + mu_1 + gamma_0 + lambdaa_1 + q_10) * DM3_1 +
        (mu_1 + lambdaa_1 * E_1 + lambdac_1 * E_1^2 + p * q_10 * E_0) * DA3 +
        (1 - p) * q_10 * DM3_0 + gamma_0 * DM3_0

      dDM2_0 <- -(lambdac_0 + mu_0 + gamma_1 + gamma_0 + lambdaa_0 + q_01) * DM2_0 +
                 (lambdaa_0 * DE_0 + 2*lambdac_0 * DE_0*E_0 + p * q_01 * DE_1) * DA3 +
                 (1 - p) * q_01 * DM2_1 + gamma_1 * DM2_1


      dDM2_1 <- -(lambdac_1 + mu_1 + gamma_0 + gamma_1 + lambdaa_1 + q_10) * DM2_1 +
                 ( lambdaa_1 * DE_1 + 2*lambdac_1 * DE_1*E_1 + p * q_10 * DE_0) * DA3 +
                 (1 - p) * q_10 * DM2_0 + gamma_0 * DM2_0

      dE_0  <- mu_0 - (mu_0 + lambdac_0 + q_01) * E_0 + lambdac_0 * E_0^2 + q_01 * E_1

      dE_1  <- mu_1 - (mu_1 + lambdac_1 + q_10) * E_1 + lambdac_1 * E_1^2 + q_10 * E_0

      dDA3  <- - (gamma_0 + gamma_1) * DA3 + gamma_0 * DM3_0 + gamma_1 * DM3_1

      list(c(dDE_0, dDE_1, dDM3_0, dDM3_1, dDM2_0, dDM2_1, dE_0, dE_1, dDA3))
    })
  }


  # Solve the system for interval [tp, t2]
  res <- master_loglik_R(parameter,
                         phy,
                         traits,
                         cond = "proper_cond",
                         root_state_weight = "proper_weights",
                         setting_calculation = NULL,
                         see_ancestral_states = TRUE,
                         atol = 1e-10,
                         rtol = 1e-10,
                         methode = "ode45",
                         use_normalization = TRUE)



  # Initial conditions
  initial_conditions2 <- c(DE_0 = res[,'DE_0'][[1]],
                           DE_1 = res[,'DE_1'][[1]],
                           DM3_0 = res[,'DM3_0'][[1]],
                           DM3_1 = res[,'DM3_1'][[1]],
                           DM2_0 = res[[1]] * res[,'DA3'][[1]],
                           DM2_1 = res[[2]] * res[,'DA3'][[1]],
                           E_0 = res[,'E_0'][[1]],
                           E_1 = res[,'E_1'][[1]],
                           DA3 = res[,'DA3'][[1]])


  # Time sequence for interval [t2, t1]
  time2 <- c(t2, t1)

  # Solve the system for interval [t2, t1]
  solution2 <- deSolve::ode(y = initial_conditions2,
                            times = time2,
                            func = interval2,
                            parms = parameter,
                            method = methode,
                            atol = 1e-10,
                            rtol = 1e-10)



  #########Interval3 [t1, t0]


  interval3 <- function(t, state, parameter) {
    with(as.list(c(state)), {

      lambdac_0 <- parameter[[1]][1]; lambdac_1 <- parameter[[1]][2]
      mu_0 <- parameter[[2]][1]; mu_1 <- parameter[[2]][2]
      gamma_0 <- parameter[[3]][1]; gamma_1 <- parameter[[3]][2]
      lambdaa_0 <- parameter[[4]][1]; lambdaa_1 <- parameter[[4]][2]
      q_01 <- parameter[[5]][1]; q_10 <- parameter[[5]][2]
      p <- parameter[[6]]

      dDM1_0 <- -(lambdac_0 + mu_0 + gamma_1 + lambdaa_0 + q_01) * DM1_0 +
                 (mu_0 + lambdaa_0 * E_0 + lambdac_0 * E_0^2 + p * q_01 * E_1) * DA1 +
                 (1 - p) * q_01 * DM1_1 + gamma_1 * DM1_1

      dDM1_1 <- -(lambdac_1 + mu_1 + gamma_0 + lambdaa_1 + q_10) * DM1_1 +
                 (mu_1 + lambdaa_1 * E_1 + lambdac_1 * E_1^2 + p * q_10 * E_0) * DA1 +
                 (1 - p) * q_10 * DM1_0 + gamma_0 * DM1_0


      dE_0 <- mu_0 - (mu_0 + lambdac_0 + q_01) * E_0 + lambdac_0 * E_0^2 + q_01 * E_1
      dE_1 <- mu_1 - (mu_1 + lambdac_1 + q_10) * E_1 + lambdac_1 * E_1^2 + q_10 * E_0

      dDA1 <- - (gamma_0 + gamma_1) * DA1 + gamma_0 * DM1_0 + gamma_1 * DM1_1


      list(c(dDM1_0, dDM1_1, dE_0, dE_1, dDA1))
    })
  }






  # Initial conditions
  initial_conditions3 <- c(
    DM1_0 = parameter[[3]][1]*solution2[, "DM2_0"][[2]] + parameter[[3]][2]*solution2[, "DM2_1"][[2]],
    DM1_1 = parameter[[3]][1]*solution2[, "DM2_0"][[2]] + parameter[[3]][2]*solution2[, "DM2_1"][[2]],
    E_0 = solution2[, "E_0"][[2]],
    E_1 = solution2[, "E_1"][[2]],
    DA1 = parameter[[3]][1]*solution2[, "DM2_0"][[2]] + parameter[[3]][2]*solution2[, "DM2_1"][[2]])


  # Time sequence for interval [t1, t0]
  time3 <- c(t1, t0)

  # Solve the system for interval [t1, t0]
  solution3 <- deSolve::ode(y = initial_conditions3,
                            times = time3,
                            func = interval3,
                            parms = parameter,
                            method = methode,
                            atol = 1e-10,
                            rtol = 1e-10)


  # Extract log-likelihood
  Lk <- solution3[, "DA1"][[2]]
  logLkb <- log(Lk)
  return(logLkb)

}



library(DAISIE)
library(ape)
library(secsse)

set.seed(1)

data("Galapagos_datalist")
datalist <- Galapagos_datalist

## select a clade in the datalist
i <- 5
phy <- DDD::brts2phylo(datalist[[i]]$branching_times[-c(1, 2)])
traits <- sample(c(0),ape::Ntip(phy),replace=TRUE)

#parameter <- list( c(0.550682, 1.550682), c (1.007281, 2.007281), c (1.683818, 0.683818),
#                  c (0.003344, 0.009344), c(0.003,0.001), 0 )

parameter <- list( c(2.546591, 0), c (2.678781, 0), c(0.009326754,0),
                  c(1.008583, 0), c(0.000,0.000), 0 )

loglik_DAISIE_trait <- DAISIE_DE_logpEC_trait1 (brts = datalist[[i]]$branching_times,
                         missnumspec = datalist[[i]]$missing_species,
                         parameter = parameter,
                         phy,
                         traits,
                         cond = "proper_cond",
                         root_state_weight = "proper_weights",
                         setting_calculation = NULL,
                         see_ancestral_states = TRUE,
                         atol = 1e-10,
                         rtol = 1e-10,
                         methode = "lsodes",
                         use_normalization = TRUE)




### This function calculate the likelihood of a single clade with a given colonization time of the island
pars2 <- c(100, 11,0,2)
pars <- c(2.546591, 2.678781, Inf, 0.009326754, 1.008583)

loglik_DAISIE <- DAISIE:::DAISIE_loglik_CS_choice( pars1 = pars,
                                  pars2 = pars2,
                                  datalist = NULL,
                                  brts = datalist[[i]]$branching_times,
                                  stac = datalist[[i]]$stac,
                                  missnumspec = datalist[[i]]$missing_species,
                                  methode = "lsodes",
                                  CS_version = 1,
                                  abstolint = 1E-16,
                                  reltolint = 1E-10,
                                  verbose = FALSE)


loglik_DAISIE_trait
loglik_DAISIE
