# function for Wright Fisher Model 
# credits: https://gist.github.com/indapa/881702

wright_fisher <- function(size, gen, freq) {
  
  #Preliminaries
  #number of chromosomes
  N = size
  
  #number of generations
  ng = gen
  f = numeric(ng) 
  
  #freq of minor alleles in the first generation
  f[1] = freq
  
  #Simulation
  #binomial sampling at each generation - determines allele count in next generations
  for (i in 2:ng){
    k = f[i-1]
    P = dbinom(0:N, N, prob = k )
    f[i] = sample( 0:N, 1, prob = P ) / N
  }
  return(f)
}



