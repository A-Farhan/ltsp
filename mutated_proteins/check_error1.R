# script to check if the bug in earlier code was because of duplicate row names when reading phmmer output as in case of multiple hits, target names would be same

# set working directory
setwd('~/ltsp/analysis/r1/cons_genes/')
# genes = list.dirs('.', recursive = F)
 genes = 'oppD'
# compare 2 methods
dum = NULL
logvec = NULL

for( g in 1:length(genes) ){
  
  # path of relevant directory for the gene
  gdir = file.path( genes[g], 'phmmer_fwd' )
  
  # list of files
  fls = list.files(gdir, full.names = T, pattern = '.fasta$') 
  
  # initiate lists to hold data from each file using 2 methods
  datals1 = sapply( 1:length(fls), function(x) NULL )
  datals2 = sapply( 1:length(fls), function(x) NULL )

  for(i in 1:length(fls)){
    
    # mthod_v1
    datals1[[i]] = tryCatch( read.table( fls[i], fill = T, stringsAsFactors = F, header = F ), 
                            error = function(e) NULL )
    
    # method_v2
    ld <- readLines( fls[i] )
    # remove all commented lines 
    ld <- ld[grep('^#', ld, invert = T)]
    datals2[[i]] <- ld
  }

  
  #remove any row with NA, extra row because of splitting statements by spaces
  for(d in 1:length(datals1) ){
    df = datals1[[d]]
     if(d==1) tmp2 = df
    if(!is.null(df)){
      df = df[,1:11]
      tmp = !apply( df, 1, anyNA )
      datals1[[d]] = df[ tmp, ]
    }
  }
  
  # count no. of entries in each table for both methods
  n_hits1 = sapply(datals1, nrow)
  n_hits1 = sapply( n_hits1, function(x) if(is.null(x)) return(0) else return(x) ) 
  n_hits2 = sapply(datals2, length)
  
  logvec = c( logvec, all( n_hits1 == n_hits2 ) )
}