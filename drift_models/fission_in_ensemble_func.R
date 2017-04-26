# Fission Model Function for a multi-type system

fission_in_ensemble <- function( frac, bdp, N, nm, ng, nuc = c('m','a') ){
 # Arguments
  # frac = fraction of population with mutant base at the mutated site
  # bdp = fraction of population picked to die or divide
  # N = population size
  # nm = number of mutated sites
  # ng = number of generations
  # nuc = nucleotide base; mutant or ancestral
  
  # randomly generate an ensemble of genomes with mutant bases at each site with the given frequency 'frac' 
  ensemble = sapply( 1:nm, function(x) sample(nuc, N, replace = T, prob = c(frac, 1 - frac) ) )
  # matrix to hold shannon's entropy values of all sites in the ensemble  
  ent <- matrix(0, nrow = ng, ncol = nm)
  # for each generation  
  for( i in 1:ng ){
    # calculate entropy at each site from the ensemble
    ent[i, ] <- apply( ensemble, 2, function(x) entH( table(x)/N ) )
    # give indices to all individuals  
    pop <- 1:N
    # sample indivudals to die  
    die <- sample( pop, bdp * N )
    # remove this set  
    pop <- setdiff(pop, die)
    # sample indivudals to divide  
    birth <- sample( pop, bdp * N )
    # add this to the existing set of indices  
    pop <- c( pop, birth )
    # use this set of indices to generate next ensemble  
    ensemble <- ensemble[ pop, ]
  }
  ent
}
