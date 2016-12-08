fission_in_ensemble <- function( frac, bdp, N, nm, ng, nuc = c('m','a') ){
  ensemble = sapply( 1:nm, function(x) sample(nuc, N, replace = T, prob = c(frac, 1 - frac) ) )
  ent <- matrix(0, nrow = ng, ncol = nm)
  for( i in 1:ng ){
    ent[i, ] <- apply( ensemble, 2, function(x) entH( table(x)/N ) )
    pop <- 1:N
    die <- sample( pop, bdp * N )
    pop <- setdiff(pop, die)
    birth <- sample( pop, bdp * N )
    pop <- c( pop, birth )
    ensemble <- ensemble[ pop, ]
  }
  ent
}