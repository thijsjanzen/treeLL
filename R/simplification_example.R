#' testing fuction, for comparison with DAISIE
#' @description
#' this function calculate the likelihood of observing a clade with specified species trait states,
#' for which the estimated maximum age of colonization is known.
#' @export
#' @param brts branching times
#' @param missnumspec number of missing species
#' @param parameter parameters
#' @param num_observed_states number of observed traits
#' @param num_hidden_states number of hidden traits
#' @param phy phy
#' @param traits traits
#' @param cond conditioning, default = "proper_cond"
#' @param root_state_weight root weight, default = "proper_weights"
#' @param setting_calculation used in ML
#' @param see_ancestral_states recover the ancestral states
#' @param atol absolute tolerance
#' @param rtol relative tolerance
#' @param methode method of integration
#'
#' library(DAISIE)
#' data("NewZealand_birds_datalist")
#' datalist <- NewZealand_birds_datalist
#' i <- 23
#' phy <- DDD::brts2phylo(datalist[[i]]$branching_times[-c(1, 2)])
#' traits <- sample(c(0,1), length(phy$tip.label), replace = TRUE)
#' sampling_fraction <- sample(c(1,1), length(phy$tip.label), replace = TRUE)
#' parameter <- list(
#'   c(2.546591, 1.2, 1, 0.2),
#'   c(2.678781, 2, 1.9, 3),
#'   c(0.009326754, 0.003, 0.002, 0.2),
#'   c(1.008583, 1, 2, 1.5),
#'   matrix(c(
#'     0,    1,    0.5,  0,
#'     0,    0,    0.002,0.005,
#'     rep(0, 8)
#'   ), nrow = 4),
#'   0
#' )
#'
#' DAISIE_DE_trait_logpEC_max_age_hidden(
#'   brts                  = datalist[[i]]$branching_times,
#'   phy                   = phy,
#'   traits                = traits,
#'   sampling_fraction     = sampling_fraction,
#'   parameter             = parameter,
#'   num_observed_states   = 2,
#'   num_hidden_traits     = 2,
#'   cond                  = "proper_cond",
#'   root_state_weight     = "proper_weights",
#'   see_ancestral_states  = TRUE,
#'   atol                  = 1e-10,
#'   rtol                  = 1e-10,
#'   methode               = "ode45",
#'   rhs_func              = loglik_hidden_rhs
#' )
DAISIE_DE_trait_logpEC_max_age_hidden_example <- function(brts,
                                                  parameter,
                                                  phy,
                                                  traits,
                                                  num_hidden_traits,
                                                  type = "max_age_hidden", # this is new
                                                  cond = "proper_cond",
                                                  root_state_weight = "proper_weights",
                                                  see_ancestral_states = TRUE,
                                                  atol = 1e-10,
                                                  rtol = 1e-10,
                                                  methode = "ode45",
                                                  rhs_func = loglik_hidden_rhs) {
  t0   <- brts[1]
  tmax <- brts[2]
  t2   <- brts[3]
  tp   <- 0

  # number of unique state
  n <- num_observed_states * num_hidden_states

  # Solve the system for interval [tp, t2]
  res <- treeLL::loglik_R_hidden(parameter,
                                 phy,
                                 traits,
                                 num_hidden_traits = num_hidden_traits,
                                 see_ancestral_states = TRUE,
                                 atol = atol,
                                 rtol = rtol)

  initial_conditions2 <- get_initial_conditions2(type, res)  # TODO: this needs to be defined

  initial_conditions2 <- matrix(initial_conditions2, nrow = 1)

  # Time sequence for interval [t2, tmax]
  time2 <- c(t2, tmax)

  # Solve the system for interval [t2, tmax]

  func_for_solution2 <- interval2      # TODO: make interval1, interval2, interval3 and interval4 in separate R file
  if (type == "max_age_hidden") func_for_solution2 <- interval3

  solution2 <- deSolve::ode(y = initial_conditions2,
                            times = time2,
                            func = func_for_solution2,
                            parms = parameter,
                            method = methode,
                            atol = atol,
                            rtol = rtol)


  solution2 <- matrix(solution2[,-1], nrow = 2) # remove the time from the result

  #########Interval3 [tmax, t0]

  initial_conditions3 <- get_initial_conditions3(type, solution2) # TODO: define this function

  # Time sequence for interval [tmax, t0]
  time3 <- c(tmax, t0)

  # Solve the system for interval [tmax, t0]
  solution3 <- deSolve::ode(y = initial_conditions3,
                            times = time3,
                            func = interval4,
                            parms = parameter,
                            method = methode,
                            atol = atol,
                            rtol = rtol)

  solution3 <- matrix(solution3[,-1], nrow = 2)

  # Extract log-likelihood
  Lk <- solution3[2,][length(solution3[2,])]
  logLkb <- log(Lk)
  return(logLkb)
}

