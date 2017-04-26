# script to convert mutated proteins profiles to matrix from which conservation at each residue of query protein could be measured

# set working directory
setwd('~/ltsp/analysis/r1/mutprot/')

# require
require(seqinr)

# vector of genes ( same as names of subdirectories )
genes <- list.dirs(full.names = F, recursive = F)

# empty list of genes to hold profiles
prots_msa <- lapply(genes, function(x) list(nam = NULL, seq = NULL) )
# each component is named by the gene
names(prots_msa) <- genes

# loop for each protein
for ( i in 1:length(genes) ) {
  
  gdir = genes[i]

  # read the profile
    # file path
  aln_fl = file.path( gdir, 'aligned.fasta' )
  # check if file exists
  if(file.exists(aln_fl)) data <- readLines(aln_fl) else next
  # extract headers
  headers_indices <- grep('>>', data )

  # get no. of sequences in the profile
  nseq = length(headers_indices)

  # extract microbe name from headers
  headers <- data[headers_indices]  
  os_names <- gsub( '.+OS=| GN=.+', '', headers )
  prots_msa[[i]]$nam <- os_names

  for(j in 1:nseq ){
    # use headers index to extract corresponding sequence from the data
    a = headers_indices[j]+1
    b = ifelse( j == nseq, length(data), headers_indices[j+1]-1 ) 
    prots_msa[[i]]$seq[j] <- paste0( data[a:b], collapse = '')
  }

  aln <- prots_msa[[i]]
  
  # a list of character vector of seq
  charseq <- lapply( aln$seq, s2c )

  # find index of query protein seq
  query_index = grep('Escherichia coli', aln$nam) 

  # find query protein length
  len <- length( charseq[[ query_index ]] )

  # create matrix of characters from aligned seq
  gmat <- matrix( unlist(charseq), ncol = len, byrow = T )

  # to reduce alignment to residues present in query
  gmat <- gmat[ , grep('-', gmat[ query_index, ], fixed = T, invert = T ) ]
  gmat <- tolower(gmat)
  
  # write matrix to file
  write.table( gmat, file.path( gdir, 'profile.mat' ), col.names = F, row.names = F, quote = F)

}
