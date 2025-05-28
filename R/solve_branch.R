#' @keywords internal
solve_branch <- function(interval_func,
                         initial_conditions,
                         time,
                         parameter,
                         methode,
                         atol,
                         rtol) {

  solution <- deSolve::ode(
      y = initial_conditions,
      times = time,
      func = interval_func,
      parms = parameter,
      method = methode,
      atol = atol,
      rtol = rtol
  )

  return(matrix(solution[, -1], nrow = 2))
}
