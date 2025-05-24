func_for_solution <- function(interval,
                              initial_conditions,
                              time,
                              parameter,
                              methode,
                              atol,
                              rtol) {
  interval_func <- get_func_interval(interval)

  if (interval == "interval2") {
    solution <- deSolve::ode(
      y = initial_conditions,
      times = time,
      func = interval_func,
      parms = parameter,
      method = methode,
      atol = atol,
      rtol = rtol
    )
  } else if (interval == "interval3") {
    solution <- deSolve::ode(
      y = initial_conditions,
      times = time,
      func = interval_func,
      parms = parameter,
      method = methode,
      atol = atol,
      rtol = rtol
    )
  } else if (interval == "interval4") {
    solution <- deSolve::ode(
      y = initial_conditions,
      times = time,
      func = interval_func,
      parms = parameter,
      method = methode,
      atol = atol,
      rtol = rtol
    )
  } else {
    stop("Unknown interval name passed.")
  }

  return(matrix(solution[, -1], nrow = 2))
}
