# function to calculate shannon's entropy
entH <- function(x){
  # removing '0' frequencies as is considered as 0
  x = x[ x > 0 ]
  # calculating shannons entropy for the rest & rounding off to 3 digits
  y = - sum( x * log2(x) ) 
  return(y)
}
