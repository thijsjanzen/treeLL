loglik_rhs <- function(t, y, parameter) {
  lambda_c <- parameter[[1]]
  lambda_a <- parameter[[2]]
  mus <- parameter[[3]]
  gammas <- parameter[[4]]
  qs <- parameter[[5]]
  p <- parameter[[6]]


  dxdt <- rep(NA, 7)

  # 1 = DE_0, 2 = DE_1
  # 3 = DM3_0, 4 = DM3_1
  # 5 =  E_0, 6 = E_1
  # 7 =  DA3


  #DE_0
  dxdt[1] <- -(lambda_c[1] + mus[1] + qs[1]) * y[1] + 2 * lambda_c[1] * y[1] * y[5] + qs[1] * y[2]

  #DE_1
  dxdt[2] <-  -(lambda_c[2] + mus[2] + qs[2]) * y[2] + 2 * lambda_c[2] * y[2] * y[6] + qs[2] * y[1]

  #DM3_0
  dxdt[3] <-  -(lambda_a[1] + lambda_c[1] + mus[1] + gammas[2] + qs[1]) * y[3]
  + (lambda_a[0] * y[5]  + lambda_c[1] * y[5] * y[5] + mus[1]
     +  p * qs[1] * y[6]) * y[7]
  +  (1 - p) * qs[1] * y[4] + gammas[2] * y[4]

  # DM3_1
  dxdt[4] =  -(lambda_c[2] + lambda_a[2] + mus[2] + gammas[1] + qs[2]) * y[4]
  + (lambda_a[2] * y[6]  + lambda_c[2] * y[6] * y[6] + mus[2]
     + p * qs[2] * y[5]) * y[7]
  + (1 - p) * qs[2] * y[3] + gammas[1] * y[3];


  # E_0
  dxdt[5]   = mus[1] - (lambda_c[1] + mus[1] + qs[1]) * y[5] +  lambda_c[1] * y[5] * y[5] + qs[1] * y[6];

  # E_1
  dxdt[6]   = mus[2] - (lambda_c[2] + mus[2] + qs[2]) * y[6] +  lambda_c[2] * y[6] * y[6] + qs[1] * y[5];

  # DA_3
  dxdt[7]  = -(gammas[1] + gammas[2]) * y[7] + gammas[1] * y[3] + gammas[2] * y[4];

  return(list(c(dxdt)))
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
    out <- c(1, 0, 0, 0, 0, 0, 1)
  }
  if (trait == 1) {
    out <- c(0, 1, 0, 0, 0, 0, 1)
  }
  return(out)
}


#' @keywords internal
master_loglik_R <- function(parameter,
                          phy,
                          traits,
                          cond = "proper_cond",
                          root_state_weight = "proper_weights",
                          setting_calculation = NULL,
                          see_ancestral_states = TRUE,
                          atol = 1e-8,
                          rtol = 1e-7,
                          methode = "ode45",
                          use_normalization = TRUE) {

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

  for(i in 1:length(ances)){
    calcul <- calcThruNodes(ances = ances[i],
                            states = states,
                            loglik = loglik,
                            forTime = forTime,
                            parameter = parameter,
                            methode = methode,
                            phy = phy,
                            rhs_func = loglik_rhs,
                            reltol = rtol,
                            abstol = atol)
    states <- calcul$states
    loglik <- calcul$loglik
    nodeN <- calcul$nodeN
  }

  loglik <- calcul$loglik
  nodeM <- calcul$node_M
  mergeBranch <- calcul$merge_branch

  if (length(nodeM) > 2 * d) nodeM <- nodeM[1:(2 * d)]

  LL <- loglik

  if (see_ancestral_states == TRUE) {
    states <- calcul$states
    numustips <- ape::Ntip(phy)
    colnames(states) <- c("DE_0", "DE_1", "DM3_0", "DM3_1", "E_0", "E_1", "DA_3")
    return(list(LL = LL, states = states))
  } else {
    return(LL)
  }
}
