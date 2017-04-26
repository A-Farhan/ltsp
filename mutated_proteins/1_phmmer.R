# Script to run PHMMER with mutated proteins as both target and query
# Not run

# set working directory
setwd('~/ltsp/analysis/r1/mutprot/')

# required package to get data from NCBI 
require( rentrez )
# required package to work with biological sequences
require( seqinr )

# get the list of non-synonymous mutations above 5 % 
data <- read.csv('nsm_above5.csv')
genes <- unique(data$gene)
# get their GI accession numbers to download protein fasta sequence from NCBI
gi <- unique(data$gi)

# get file names for all reference proteomes of gamma proteobacteria
refprots = list.files( <path_to_proteomes_fasta_files> )

# load all reference proteomes into a list
ls_prot <- unlist( lapply( file.path( 'gambac_refprot', refprots ), function(x) 
    read.fasta(file = x, seqtype = 'AA', as.string = T ) ), recursive = F )

# start the program
for ( i in 1:length(genes) ) {
  
  # make directory for each gene
  gdir = file.path( genes[i] )
  dir.create( gdir, showWarnings = F )
  
  # get protein sequence from NCBI and save to file
  write( entrez_fetch(db = 'protein', id = gi[ i ], rettype = 'fasta'), file.path(gdir, 'proseq.fasta' ) )
  
  # create directory to hold PHMMER Forward results
  resdir = file.path(gdir, 'phmmer_fwd')
  dir.create( resdir, showWarnings = F )
  
  # make each reference proteome the target database for homolog search
  # and save results in tabulated form
  # use e-value inclusion threshold of 1e-10
  # save the run errors to a file
  
  # Run PHMMER-FWD
  for( j in 1:length(refprots) )
    system( paste0( 'phmmer --cpu 32 --incE 1e-10 --tblout', ' ', file.path( resdir, refprots[j] ), ' ',
                    file.path( gdir, 'proseq.fasta'), ' ', file.path( 'gambac_refprot', refprots[j] ), ' ', 
                    '2>> ', file.path( resdir, 'err.log') ),
            ignore.stdout = F, ignore.stderr = F, intern = F, wait = T )

  # create directory to hold PHMMER Reverse results
  resdir = file.path(gdir, 'phmmer_rev')
  dir.create( resdir, showWarnings = F )
  
  # test all the proteins of each proteome to be homolog for the given protein sequence as target 
  # and save results in tabulated form
  # use e-value inclusion threshold of 1e-10
  # save the run errors to a file
  
  # Run PHMMER-REV
  for( k in 1:length(refprots) )
    system( paste0( 'phmmer --cpu 32 --incE 1e-10 --tblout', ' ', file.path( resdir, refprots[k] ), ' ',
                    file.path( 'gambac_refprot', refprots[k] ), ' ', file.path( gdir, 'proseq.fasta'), ' ', 
                    '2>> ', file.path( resdir, 'err.log') ),
            ignore.stdout = F, ignore.stderr = F, intern = F, wait = T)
}

