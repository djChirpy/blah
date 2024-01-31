
library(ggplot2)
library(reshape2)
library(ggraph)

args = commandArgs(trailingOnly=TRUE)
inFile = args[1]

microPerc <- read.table(inFile, header = T, quote = '', stringsAsFactors = F, sep = '\t')
noDog <- subset(microPerc, select = -c(Chordata_percent,Chordata_reads, unclassified_reads, unclassified_percent, Chordata_microBiomePercentage, unclassified_microBiomePercentage))
keep <- grepl('_microBiomePercentage', colnames(noDog))
noDogKeep <- noDog[keep]
noDogKeep <- noDogKeep[, order(colSums(noDogKeep))]
colnames(noDogKeep) <- gsub('_microBiomePercentage','',colnames(noDogKeep))

display <- c('unclassified', 'Proteobacteria', 'Bacteroidetes', 'Firmicutes', 'Actinobacteria', 'Spirochaetes', 'Apicomplexa', 'Fusobacteria', 'Cyanobacteria', 'Tenericutes', 'Candidatus.Saccharibacteria', 'Planctomycetes', 'Euryarchaeota','NA')
displayKeep <- colnames(noDogKeep) %in% display
noDogKeep <- noDogKeep[displayKeep]
noDogKeep$other <- 100-rowSums(noDogKeep)
noDogKeep$name <- row.names(noDogKeep)
meltNoDog <- melt(noDogKeep)
colorsTest <- c('#6f19b2','#00bcff','#990000')
colorSet <- colorRampPalette(colorsTest)(13)
colorSet <- c(colorSet, '#808080')
colScale <- scale_colour_manual(name = display, values = colorSet)
fillScale <- scale_fill_manual(name = display, values = colorSet)
edgeScale <- scale_edge_color_manual(name = display, values = colorSet)

pdf(paste(inFile, '.pdf'), height = 7, width = 7)
ggplot(meltNoDog, aes(x=name, y = value, fill=variable, width=1)) + geom_bar(stat = 'identity', position = 'stack',size=0) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y=element_blank()) + scale_fill_manual(values = colorSet, name='Phylum')
dev.off()



