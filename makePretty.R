library(ggplot2)
library(reshape2)
library(ggraph)
library(igraph)

args = commandArgs(trailingOnly=TRUE)
inFile = basename(args[1])
name = substr(inFile,1,nchar(inFile)-8)
edgeFile = paste(name, ".edges", sep = "")
vertFile = paste(name, ".vertices", sep = "")

percentages <- read.table(inFile, header = T, quote = '', stringsAsFactors = F, sep = '\t',colClasses = c('character',rep('numeric',3)))
data <- data.frame(
  category = c("dog", "Microbiome", "Unclassified"),
  count = c(percentages$dogPercent[1], percentages$microPercent[1], percentages$unclassified[1])
)
data$ymax = cumsum(data$count)
data$ymin = c(0, head(data$ymax, n=-1))
data$category <- as.factor(data$category)
pdf(paste(name, '_summary.pdf', sep = ''), height = 7, width = 7)
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) + 
  scale_fill_manual(values = c('dog' = '#c0d434',
                               'Microbiome' = '#40bcec',
                               'Unclassified' = '#f04494')) +
  theme_void() + 
  theme(legend.position = 'none') #+
  #colScale
#ggplot(meltNoDog, aes(x=name, y = value, fill=variable, width=1)) + geom_bar(stat = 'identity', position = 'stack',size=0) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y=element_blank()) + scale_fill_manual(values = colorSet, name='Phylum')
dev.off()

edges <- read.table(edgeFile,
                    stringsAsFactors = F,
                    sep ='\t',
                    quote = '',
                    fill = T,
                    header = T)
vertices <- read.table(vertFile,
                       stringsAsFactors = F,
                       sep ='\t',
                       quote = '',
                       fill = T)
colnames(edges) <- c('from','to','score','phylum')
dacol <- c('#693793','#ed4096','#52c6e0','#b9d432','#fce000','#ba9adb','#fbd2d7','#1b75bb','#50c878','#f8de7e','#66144b','#d31940','#262261','#008080','#c49102')
#colExpand <- rep(dacol,100)
subset = subset.data.frame(edges, select = -c(phylum))
subset = na.omit(subset)
cutoff = quantile(subset$score, c(.85))
subset = subset[subset$score > cutoff,]
#AscomycotaScore = subset[subset$to == 'Ascomycota', 3]
dogScore = subset[subset$to == 'Canis',3]
#subset = subset[subset$score != AscomycotaScore,]
subset = subset[subset$score != dogScore,]
nodes = unique(c(unique(subset$from), unique(subset$to)))
nodesdf = data.frame(nodes)
colnames(vertices) <- c('name','score','phylum')
testMerge = merge(nodesdf, vertices, by.y = 'name', by.x = 'nodes')
colnames(testMerge) <- c('name','score','phylum')
testMerge$phylum = ifelse(is.na(testMerge$phylum),'other',testMerge$phylum)
testMerge$score = ifelse(testMerge$score > 1000000, 0, testMerge$score)
testMerge$score = ecdf(testMerge$score)(testMerge$score)
testAgg <- aggregate(testMerge$score, by = list(phylum=testMerge$phylum), FUN=sum)
testAgg$name = ifelse(ecdf(testAgg$x)(testAgg$x) > 1-(14.5/length(testAgg$phylum)), testAgg$phylum, 'other')
testMerge$phylum = ifelse(testMerge$phylum %in% testAgg$name, testMerge$phylum, 'other')
mygraph <- graph_from_data_frame(subset, vertices=testMerge)
pdf(paste(name, '_flower.pdf',sep = ''), height=10, width=14)
ggraph(mygraph, layout = 'kk') +
  geom_node_circle(aes(r = score, fill = as.factor(phylum))) +
  geom_edge_link(aes(alpha = 0.1), show.legend = FALSE) +
  theme_void() +
  scale_fill_manual(values = dacol) +
  labs(fill = "Phylum")
dev.off()
