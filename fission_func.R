# Fission Model Function in bi-allele system
fission <- function(size, n_gen, freq, rate ){

  # birth-death probability  
  bdp = rate
  # no. of generations
  ng = n_gen
  # population size
  N = size
  # number of individuals dying/birthing
  n_bd = bdp * N
  # vector of frequency of mutant in each generation
  x <- numeric(ng)
  # initial frequency of mutant 
  x[1] <- freq
  # next generation onwards    
  for (i in 2:ng ){
    # get freq of previous generation
    f = x[i-1]
    # convert to numbers
    fN = round(f * N)
    # probability of birthing
    prob2 <- dhyper(x = 0:fN, m = fN, n = N - fN, k = n_bd )
    # sample number of birthing individuals
    B = sample( 0:fN, 1, prob = prob2 ) 
    left = fN - B
    # if nothing is left then numbers in next generation 
    # would be double of birthing indivudals in this gen
    if( left == 0 )
      x[i] = B * 2 / N
    # else, dying indivudals are sampled similarly out of the rest
    # & freq in next generation would be double of birthing individuals and those remaining 
    # after death process 
    else{
      prob0 <- dhyper(x = 0:left, m = left, n = N - n_bd - left, k = n_bd )
      D = sample( 0:left, 1, prob = prob0 ) 
      x[i] = ( B * 2 + left - D ) / N
    }
  }
  return(x)
}  



