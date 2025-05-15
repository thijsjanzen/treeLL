#' @keywords internal
loglik_rhs <- function(t, state, parameter) {
  with(as.list(c(state, parameter)), {

    DE_0 = state[1]
    DE_1 = state[2]
    DM3_0 = state[3]
    DM3_1 = state[4]
    E_0 = state[5]
    E_1 = state[6]
    DA3 = state[7]

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


#' @keywords internal
calcThruNodes <- function(
    ances,
    states,
    loglik,
    forTime,
    parameter,
    methode,
    phy,
    rhs_func,
    reltol,
    abstol
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
                           rtol = reltol,
                           atol = abstol,
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
  lambda_c <- parameter[[1]]
  combined_state[1] <- lambda_c[1] * nodeN[1] * nodeM[1];
  combined_state[2] <- lambda_c[2] * nodeN[2] * nodeM[2];

  states[focal,] <- combined_state
  return(list(states = states,
              loglik = loglik,
              combined_state = combined_state,
              nodeM = nodeM,
              nodeN = nodeN))
}

#' @keywords internal
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


#' @keywords internal
simple_loglik_R <- function(parameter,
                            phy,
                            traits,
                            cond = "proper_cond",
                            root_state_weight = "proper_weights",
                            setting_calculation = NULL,
                            see_ancestral_states = TRUE,
                            atol = 1e-8,
                            rtol = 1e-7,
                            methode = "ode45",
                            use_normalization = TRUE,
                            rhs_func = loglik_rhs) {

  number_of_lineages <- length(phy$tip.label)

  states <- matrix(nrow = number_of_lineages + phy$Nnode,
                   ncol = 7,
                   data = NA)

  for (i in 1:length(traits)) {
    states[i, ] <- calc_init_state(traits[i])
  }

  phy$node.label <- NULL
  split_times <- sort(event_times(phy), decreasing = FALSE)
  ances <- as.numeric(names(split_times))
  forTime <- cbind(phy$edge, phy$edge.length)

  d <- ncol(states) / 2
  loglik <- 0

  for (i in 1:length(ances)){
    calcul <- calcThruNodes(ances = ances[i],
                            states = states,
                            loglik = loglik,
                            forTime = forTime,
                            parameter = parameter,
                            methode = methode,
                            phy = phy,
                            rhs_func = rhs_func,
                            reltol = rtol,
                            abstol = atol)
    states <- calcul$states
    loglik <- calcul$loglik
    nodeN <- calcul$nodeN
  }

  prob_states <- calcul$combined_state
  prob_states <- matrix(prob_states, nrow = 1,
                        dimnames = list(NULL, c("DE_0", "DE_1", "DM3_0", "DM3_1", "E_0", "E_1", "DA3")))
  return(prob_states)
}
