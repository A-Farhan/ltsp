# make boxplots of conservation scores of all residues 
# of proteins with non-synonymous mutations

# set working directory
setwd('~/ltsp/analysis/r1/mutprot/')

# require script
source('~/ltsp/analysis/scripts/func/cons_score_func.R')

# directory of mutated proteins
genes <- list.dirs(full.names = F, recursive = F)

# exclude ychS & ycbR as no ortholgs present
genes <- genes[!genes %in% c('ychS','ycbR') ]

# load all profile matrices into list
profiles <- lapply(genes, function(x)
  read.table( file.path(x, 'profile.mat' ), sep = ' ', header = F, stringsAsFactors = F ) )
# give names to list elements
names(profiles) <- genes
# calculate conservation scores for each profile
prot_cscores <- lapply(profiles, cons_score)

# read data on NSM 
data <- read.csv('nsm_above5.csv', header = T, stringsAsFactors = F )
# extract columns for 'gene' & 'annotation'
data <- data[,c(1,4)]
# extract data for 'genes'
data <- data[ data[,1] %in% genes, ]
# modify annotation column
tmp <- data[,2]
data[,2] <- sapply( strsplit( gsub('[^[:alnum:]|*]', '_', tmp ), '_' ),'[', 1 )
# extract residue position from annotation 
data$pos <- as.numeric( gsub('[^[:digit:]]','', data[,2] ) )
# save original and mutant residue as separate columns
tmp <- gsub('[^[:alpha:]|*]','', data[,2] )
tmp <- do.call( rbind, strsplit(tmp,'') )
data[,4:5] <- t( apply( tmp, 1, tolower) )
# save binary values for mutation being conservative in another column
tmp <- apply( data[ ,4:5] , 1, function(x) sapply(func_class, function(y) sum( x %in% y ) ) )
data$consv <- as.numeric( apply(tmp, 2, function(x) 2 %in% x) )
# split data into list by gene names
nsm <- split( data, data$gene )
# extract mutated residue positions into lists for each gene
nsm_pos <- lapply(nsm, '[', 3 )
# initiate a list of conservation score for above list of positions
nsm_score <- lapply(nsm_pos, function(x) NULL )
# loop for each protein
for(i in 1:length(nsm_pos) ){
  # get mutated positions of the protein from above list
  prot <- nsm_pos[[i]]
  # initiate vector to hold conservation score
  score <- NULL
  # loop for each position
  for(j in 1:nrow(prot)){
    # get position
    pos <- prot[j,1]
    # extract its conservation score 
    tmp <- prot_cscores[[i]][pos] 
    # append it to the score vector
    score <- c(score, tmp)
  }
  # save score vector into the above list for scores
  nsm_score[[i]] <- score  
}
# save the list as another column in data 
data$cscores <- round( unlist(nsm_score), 2 )

# write data to file
# columns: gene, annotation, conservative bit, conservation score
write.csv( data[,c(1,2,6,7)], 'nsm_cscores.csv', row.names = F )

