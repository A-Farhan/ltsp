wf_in_ensemble <- function( N, ng, nm, frac, nuc = c('m', 'a') ){
  
  # probability of choosing the same parent in WF Model,
  Ppar <- round(dpois(x = 0:N, lambda = 1) * N)
  Ppar <- Ppar[Ppar > 0]

  ensemble = sapply( 1:nm, function(x) sample(nuc, N, replace = T, prob = c(frac, 1 - frac) ) )
  ent <- matrix( 0, nrow = ng, ncol = nm)

  # code to execute for each generation
  for (i in 1:ng){
    # shannon entropy at each position in each generation
    ent[i, ] <- apply( ensemble, 2, function(x) entH( table(x) / N ) )
    
    # list to hold indices of parents leaving n offsprings
    ls_par <- list()
    # temporary population
    tpop <- 1:N
    for (j in 1:length(Ppar)) {
      ls_par[[j]] <- rep( sample( tpop, Ppar[j]), j - 1 )
        tpop <- setdiff(tpop, ls_par[[j]])
    }
    ensemble <- ensemble[ unlist( ls_par[-1] ), ]
  }
  ent
}
