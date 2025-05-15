#' testing fuction, for comparison with DAISIE
#' @description
#' This function calculates something we can verify with DAISIE
#' @export
#' @param brts branching times
#' @param missnumspec number of missing species
#' @param parameter parameters
#' @param phy phy
#' @param traits traits
#' @param cond conditioning, default = "proper_cond"
#' @param root_state_weight root weight, default = "proper_weights"
#' @param setting_calculation used in ML
#' @param see_ancestral_states recover the ancestral states
#' @param atol absolute tolerance
#' @param rtol relative tolerance
#' @param methode method of integration
DAISIE_DE_logpEC_trait1_hidden_old <- function (brts,
                                           missnumspec,
                                           parameter,
                                           phy,
                                           traits,
                                           num_hidden_traits,
                                           cond = "proper_cond",
                                           root_state_weight = "proper_weights",
                                           see_ancestral_states = TRUE,
                                           atol = 1e-10,
                                           rtol = 1e-10,
                                           methode = "ode45",
                                           rhs_func = loglik_hidden_rhs) {
  t0 <- brts[1]
  t1 <- brts[2]
  t2 <- brts[3]
  tp <- 0
  ti <- sort(brts)
  ti <- ti[1:(length(ti)-2)]

  #########Initial conditions [tp, t2]
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
  res <- treeLL::loglik_R_hidden(parameter,
                                 phy,
                                 traits,
                                 num_hidden_traits = num_hidden_traits,
                                 cond = "proper_cond",
                                 root_state_weight = "proper_weights",
                                 see_ancestral_states = TRUE,
                                 atol = 1e-10,
                                 rtol = 1e-10,
                                 use_R_version = TRUE,
                                 rhs_func = loglik_hidden_rhs)

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
