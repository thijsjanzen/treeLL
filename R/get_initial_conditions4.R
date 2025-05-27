get_initial_conditions4 <- function(status, solution, parameter, trait_mainland_ancestor, num_observed_states, num_hidden_states)

{
  n <- num_observed_states * num_hidden_states
  num_unique_states <- n
  if (status == 2 || status == 3 || status == 4)
  {

    if (trait_mainland_ancestor == "FALSE")
    {  initial_conditions4 <- c(rep( sum(parameter[[3]] * (solution[2,][(n + 1):(n + n)])), n), ### DM1: select DM2 in solution2
                                solution[2,][(n + n + n + 1):(n + n + n + n)],         ### E: select E in solution2
                                sum(parameter[[3]] * (solution[2,][(n + 1):(n + n)])))          ### DA1: select DM2 in solution2

    }
    #if the trait state of the species at the sten is known
    else if(trait_mainland_ancestor == trait_mainland_ancestor)
    {
      pos <- c((num_hidden_states*trait_mainland_ancestor + 1), num_hidden_states + trait_mainland_ancestor* num_hidden_states)
      initial_conditions4 <- c(rep (sum(parameter[[3]][pos] * (solution[2,][(n + 1):(n + n)])[pos])/2,n ), ### DM1: select DM2 in solution2
                               solution[2,][(n + n + n + 1):(n + n + n + n)],                                                                       ### E: select E in solution2
                               sum(parameter[[3]][pos] * (solution[2,][(n + 1):(n + n)])[pos]/2))          ### DA1: select DM2 in solution2

    }
  }
  else if(status == 1 || status == 5 || status == 6)

  {
    initial_conditions4 <- c(solution[2,][(n + 1):(n + n)],                                 ### DM1: select DM1 in solution1
                             solution[2,][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution1
                             solution[2,][length(solution[2,]) - 1])                        ### DA1: select DA2 in solution1

  }

  else if(status == 8 || status == 9)

  {
    initial_conditions4 <- c(solution[2,][(n + 1):(n + n)],                                 ### DM1: select DM2 in solution3
                             solution[2,][(n + n + n + n + 1):(n + n + n + n + n)],         ### E: select E in solution3
                             solution[2,][length(solution[2,]) - 1])                        ### DA1: select DA2 in solution3

  }
  return(matrix(initial_conditions4, nrow = 1))

}
