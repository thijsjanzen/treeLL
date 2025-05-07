
all_trees <- readRDS("/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/Ornela/multi_trees_woody_16042025.rds")
all_traits <- readRDS("/Users/thijsjanzen/Library/CloudStorage/Dropbox/projects/Ornela/multi_traits_woody_16042025.rds")

focal_tree_size <- 3

focal_index <- 1
for (i in 1:length(all_traits)) {
  if (length(all_traits[[i]]) == focal_tree_size) {
    focal_index <- i
    break
  }
}

phy <- all_trees[[focal_index]]
traits <- all_traits[[focal_index]]
require(ape)
plot(phy)

lambda_c <- c(1, 1)
lambda_a <- c(0.0, 0.0)
mus <- c(0, 0)
gammas <- c(0.01, 0.01)
qs <- c(0, 0)
p <- 0

parameters <- list()
parameters[[1]] <- lambda_c
parameters[[2]] <- lambda_a
parameters[[3]] <- mus
parameters[[4]] <- gammas
parameters[[5]] <- qs
parameters[[6]] <- p


res <- calc_loglik(parameter = parameters,
            phy = phy,
            traits = traits,
            sampling_fraction = c(1, 1),
            see_ancestral_states = TRUE,
            use_normalization = FALSE)
res






