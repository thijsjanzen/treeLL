get_initial_conditions2 <- function(status, res, trait, num_observed_states, num_hidden_states)

{
  n <- num_observed_states * num_hidden_states
  num_unique_states <- n

  if (status == 2 && length(brts)>2)
  {

    initial_conditions2 <- c(res[1:n],                      ## DE
                             (res[1:n]) * res[length(res)], ## DM2
                             res[(n + 1):(n + n)],          ## DM3
                             res[(n + n + 1):(n + n + n)],  ## E
                             res[length(res)])              ## DA3

  }
  else if (status == 2 && length(brts) == 2)
  {

    calc_init_state_hidden <- function(trait,
                                       num_unique_states,
                                       num_hidden_states) {

      DE  <- rep(0, num_unique_states)
      DM2 <- rep(0, num_unique_states)
      DM3 <- rep(0, num_unique_states)
      E   <- rep(0, num_unique_states)
      DA3 <- 1


      DE[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- sf
      E[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)]  <- sf

      return( c(DE, DM2, DM3, E, DA3))
    }


    initial_conditions2 <- calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  }


  else if (status == 3 && length(brts) == 2 )
  {

    calc_init_state_hidden <- function(trait,
                                       num_unique_states,
                                       num_hidden_states) {

      DE  <- rep(0, num_unique_states)
      DM2 <- rep(0, num_unique_states)
      DM3 <- rep(0, num_unique_states)
      E   <- rep(0, num_unique_states)
      DA3 <- 0


      DE[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- sf
       E[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- sf
      DM3[c((num_hidden_states*trait_mainland_ancestor + 1), num_hidden_states + trait_mainland_ancestor* num_hidden_states)] <- 1

      return( c(DE, DM2, DM3, E, DA3))
    }

    initial_conditions2 <- calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  }

  else if (status == 4)
  {
    calc_init_state_hidden <- function(trait,
                                       num_unique_states,
                                       num_hidden_states) {

      DE  <- rep(0, num_unique_states)
      DM2 <- rep(0, num_unique_states)
      DM3 <- rep(0, num_unique_states)
      E   <- rep(0, num_unique_states)
      DA3 <- 1


      DM2[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- 1

      return( c(DE, DM2, DM3, E, DA3))
    }



    initial_conditions2 <- calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  }

  else if(status == 8)

  {
    calc_init_state_hidden <- function(trait,
                                       num_unique_states,
                                       num_hidden_states) {

      DE  <- rep(0, num_unique_states)
      DM2 <- rep(0, num_unique_states)
      DM3 <- rep(0, num_unique_states)
      E   <- rep(0, num_unique_states)
      DA3 <- 1


      DM2[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- 1

      return( c(DE, DM2, DM3, E, DA3))
    }

    initial_conditions2 <-   calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  }
  else if(status == 9)

  {
    calc_init_state_hidden <- function(trait,
                                       num_unique_states,
                                       num_hidden_states) {

      DE  <- rep(0, num_unique_states)
      DM2 <- rep(0, num_unique_states)
      DM3 <- rep(0, num_unique_states)
      E   <- rep(0, num_unique_states)
      DA3 <- 1


      DE[c((num_hidden_states*trait + 1), num_hidden_states + trait* num_hidden_states)] <- 1

      return( c(DE, DM2, DM3, E, DA3))
    }

    initial_conditions2 <-   calc_init_state_hidden(trait, num_unique_states, num_hidden_states)

  }
  return(matrix(initial_conditions2, nrow = 1))

}
