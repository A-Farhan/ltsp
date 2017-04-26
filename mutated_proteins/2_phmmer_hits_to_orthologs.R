# Script to get orthologs of each mutated protein (with NSM ) in the study from phmmer output 
# using bidirectional best-hit

# set working directory
setwd('~/ltsp/analysis/r1/mutprot/')

# capture run time pt1
ptm <- proc.time()

# required package to work with biological sequences
require( seqinr )

# vector of replicates
rep = c( 1,2,3,6 )
# vector of days 
day = c(0,2,4,6,10,16,22,28)

# get the list of non-synonymous mutations above 5 % 
data <- read.csv('nsm_above5.csv', stringsAsFactors = F)
genes <- unique(data$gene)

# get all reference proteomes in a list
# get file names for all reference proteomes of gamma proteobacteria
refprots = list.files( <path_to_reference_proteomes_fasta_files> )

ls_prot <- unlist( lapply( file.path( 'gambac_refprot', refprots ), function(x) 
    read.fasta(file = x, seqtype = 'AA', as.string = T ) ), recursive = F )


for ( i in 1:length(genes) ) {
  
  gdir = file.path(genes[i])
  dr = file.path(gdir, 'phmmer_fwd')
  files = list.files( dr, pattern = '.fasta$') 
  # get the resultant tables into a list
  lsF <- sapply( files, function(x) readLines( file.path( dr, x ) ), USE.NAMES = T, simplify = F )
  # remove all commented lines 
  lsF <- sapply(lsF, function(X) X[grep('^#', X, invert = T)])
  # remove those OS records in which no ortholog is reported
  lsF <- lsF[sapply(lsF, length) >= 1 ]
  
  dr = file.path(gdir, 'phmmer_rev')
  files = list.files( dr, pattern = '.fasta$') 
  # get the resultant tables into a list
  lsR <- sapply( files, function(x) readLines( file.path( dr, x ) ), USE.NAMES = T, simplify = F )
  # remove all commented lines 
  lsR <- sapply(lsR, function(X) X[grep('^#', X, invert = T)])
  # remove those OS records in which no ortholog is reported
  lsR <- lsR[sapply(lsR, length) >= 1 ]
  
  # only consider OS records present in both forward and reverse output of PHMMER
  common_os = intersect( names(lsF), names(lsR) )
  lsF <- lsF[common_os]
  lsR <- lsR[common_os]
  
  # get a list of evalues for each ortholog in each OS
  f_data <- sapply( lsF, function(E) lapply( strsplit(E, ' '), function(x) x[!x==''] ) )
  evalues_F <- sapply( f_data, function(X) as.numeric( sapply(X, '[', 5) ) )
  
  r_data <- sapply( lsR, function(E) lapply( strsplit(E, ' '), function(x) x[!x==''] ) )
  evalues_R <- sapply( r_data, function(X) as.numeric( sapply(X, '[', 5) ) )
  
  # Bi-directional Best hit
  targets <- NULL
  for (i in 1:length(common_os)){
    best_match_f <- f_data[[i]] [ match( min( evalues_F[[i]] ), evalues_F[[i]] ) ] ;
    min_ef = sapply(best_match_f, '[', 1)
    best_match_r <- r_data[[i]] [ match( min( evalues_R[[i]] ), evalues_R[[i]] ) ];
    min_er = sapply(best_match_r, '[', 3)
    matching_prot = intersect( min_ef, min_er )
    if(length(matching_prot)>=1){
      ortho_prot = matching_prot[1]
      targets[common_os[i]] = ortho_prot
    }
  }
  
  # create directory to hold fasta files of homologs	
  resdir <- file.path( gdir, 'homologs')
  dir.create(resdir, showWarnings = F )
  
  # make fasta files for each homolog protein sequence and save under name of the reference proteome
  sapply( names(targets), function(x) write.fasta(sequences = ls_prot[ targets[x] ],
                                                  names = getAnnot( ls_prot[ targets[x] ] ), 
                                                  file.out = file.path( resdir, x ), as.string = T ) )
  
  # create directory to hold NEEDLE results
  resdir <- file.path( gdir, 'needle')
  dir.create(resdir, showWarnings = F )
  
  # run NEEDLE for global alignment to get percent identity  
  
  for( k in names(targets) ) 
    system( paste0( 'needle -auto Y -outfile', ' ', file.path( resdir, k ), ' ',
                    file.path( gdir, 'proseq.fasta'), ' ', file.path( gdir, 'homologs', k ) ),
            ignore.stdout = F, ignore.stderr = F, intern = F, wait = T)
  
  
  percent_id <- cbind( sapply( names(targets), function(x) {
    tmp <- grep( 'Identity', readLines( file.path( resdir, x ) ), fixed = T, value = T );
    return( sub( "\\%).*", "", sub( ".*\\(", "", tmp ) ) ) 
  } ), targets )  	
  
  # consider only homologs with at least 50 % identity
  percent_id <- percent_id[ as.numeric(percent_id[ , 1]) >= 50, ]
  selhm <- rownames( percent_id )
  
  # get selected homologs fasta sequence in one file
  A = file.path( gdir, 'sel_hm.fasta' )
  file.create(A)
  for( i in 1:length(selhm) ){
    B = file.path( gdir, 'homologs', selhm[i] )
    file.append( A, B)
  }
}

