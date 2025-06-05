initvals <- c(2.546591, 2.678781, 0.009326754, 1.008583, runif(4))
idparslist <- list()
idparslist[[1]] <- c(1, 1, 1, 1) # lambda_c
idparslist[[2]] <- c(2, 2, 2, 2) # mu
idparslist[[3]] <- c(3, 3, 3, 3) # gamma
idparslist[[4]] <- c(4, 4, 4, 4) # lambda_a

idparslist[[5]] <- matrix(0, 4, 4)
# observed state transitions
idparslist[[5]][1, 3] <- 5 # 0 -> 1
idparslist[[5]][2, 4] <- 5 # 0 -> 1
idparslist[[5]][3, 1] <- 6 # 1 -> 0
idparslist[[5]][4, 2] <- 6 # 1 -> 0
colnames(idparslist[[5]]) <- c("0A", "0B", "1A", "1B")
rownames(idparslist[[5]]) <- colnames(idparslist[[5]])
# hidden state transitions
idparslist[[5]][1, 2] <- 7 # A -> B
idparslist[[5]][2, 1] <- 8 # B -> A
idparslist[[5]][3, 4] <- 7 # A -> B
idparslist[[5]][4, 3] <- 8 # B -> A

idparslist[[6]] <- 0 #p

parameter <- list()
parameter[[1]] <- secsse::fill_in(idparslist[[1]], initvals)
parameter[[2]] <- secsse::fill_in(idparslist[[2]], initvals)
parameter[[3]] <- secsse::fill_in(idparslist[[3]], initvals)
parameter[[4]] <- secsse::fill_in(idparslist[[4]], initvals)
parameter[[5]] <- secsse::fill_in(idparslist[[5]], initvals)
parameter[[6]] <- 0





initial_conditions2 <- c(0,    0,    0 ,   0,    1,    1,    0,    0 ,   0,     0,     0,     0,     0 ,    0,     0,     0 ,    1)
time2 <- c(0.00, 0.34)

num_hidden_states <- 2


test_func <- function(rcpp_flag) {
    treeLL::solve_branch(interval_func = interval2,
                           initial_conditions = initial_conditions2,
                            time = time2,
                            parameter = parameter,
                            methode = "ode45",
                            atol = 1e-9,
                            rtol = 1e-9,
                            use_Rcpp = rcpp_flag)
}

test_func(0)
