#' branch solving
#' @description
#' solve along branch
#' @inheritParams default_params_doc
#' @export
solve_branch <- function(interval_func,
                         initial_conditions,
                         time,
                         parameter,
                         methode,
                         atol,
                         rtol,
                         use_R = TRUE) {
  solution <- c()
  if (use_R == TRUE) {

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
  } else {
    interval_name <- as.character(substitute(interval_func))
    solution <- solve_branch_cpp(interval_name,
                                 initial_conditions,
                                 time,
                                 parameter,
                                 methode,
                                 atol,
                                 rtol)
    return(solution)
  }


}
