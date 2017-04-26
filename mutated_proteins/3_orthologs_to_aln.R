# script to build profiles of mutated proteins from orthologs identified 
setwd('~/ltsp/analysis/r1/mutprot/')

genes <- list.dirs(full.names = F, recursive = F)

for ( i in 1:length(genes) ) {
  
  gdir = file.path(genes[i])
  
  # Run MUSCLE with selected homologs passing identity cutoff to do Local alignment
  system( paste0( 'muscle -in', ' ', file.path( gdir, 'sel_hm.fasta' ), ' ', '-out', ' ', 
                file.path( gdir, 'aligned.fasta' ) ),
        ignore.stdout = T, ignore.stderr = T, intern = F, wait = T)
}











