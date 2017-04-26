# script to plot consrvation scores for residues of protein mutated in ltsp

# set working directory
setwd('~/ltsp/analysis/r1/mutprot/')
# script to perform analysis
source('5_cal_cscores.R')

# read data on consrvation scores of protein residues mutated in ltsp
data <- read.csv('nsm_cscores.csv', header = T, stringsAsFactors = F )
# split conservation scores into list for conservative & non-conservative mutations
tab <- cbind( factor(data[,1]), data[,c(4,3)] )
tab <- split( tab[,1:2], tab[,3] )

# plotting
pdf('boxplot_prot_cscores.pdf', useDingbats = F, width = 9 )
# define colors
cb_ggp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette(cb_ggp2)

par.org <- par()
par( bty = 'n', mar = c(6, 4, 2, 0), mgp = c(4,1,0) )
boxplot( prot_cscores, pch = '.', axes = F, varwidth = T, xlab = 'Genes', ylab =  '', 
         cex.lab = 1.4, border = 'green')
axis(1, at=seq(1, 22, by=1), labels = F )
text( seq(1, 22, by=1), labels = genes, par('usr')[3] - 0.06, srt = 45, pos = 1, xpd = T, cex = 1.3 )
axis(2, at=seq(0, 1, by=0.2), labels = T, cex = 1.1, line = -1 )
mtext(text = 'Conservation score', side = 2, line = 2, cex = 1.4 )
points( tab[[1]], pch = 24, bg = 2 )
points( tab[[2]], pch = 25, bg = 3 )
segments(x0 = 0.5, x1 = 22.5, y0 = 0.9, y1 = 0.9, lty = 2 )
# mtext(text = 'the dashed line marks the average median score', side = 1, line = 5, cex = 1)
x_i <- factor(data$gene)
y_i <- data$cscores
labs = data$annotation

data$gene = factor(data$gene)
tmp_ls <- split(x = data, f = data[,1])
for(i in 1:length(tmp_ls)){
  gene = tmp_ls[[i]]
  labs = gene[,2]
  coords = gene[,c(1,4)]
  x = coords[,1]
  ys = coords[,2]
  ix1 = which(ys==1)
  if(length(ix1)>1)
    for(j in 2:length(ix1)) ys[ix1[j]] = ys[ix1[j-1]]-0.04
  text(x,ys, labels = labs, pos = 1, srt = 45, cex = 0.8, xpd = T, col = 'darkblue' )
}

dev.off()
