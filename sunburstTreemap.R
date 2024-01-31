library(sunburstR)
pathToGraph <- read.table('/Volumes/seq_vgb/swofford/pathseq/nellie/NelGood_sal_1.bam.pathseq.txt', header = T, stringsAsFactors = F, sep = '\t', quote = '', comment.char = '')
pathToGraph$taxonomy <- gsub('\\|','-',pathToGraph$taxonomy)
pathToGraph <- pathToGraph[pathToGraph$taxonomy != 'root',]
forPlot <- pathToGraph[c('taxonomy','unambiguous')]
sunburst(forPlot,legend=list(w=200,h=20))

pathToGraph$taxonomy <- gsub('\\|',',',pathToGraph$taxonomy)

#p <- treemap(data, index=c('Level1','Level2','Level3'), vSize="count", type="index", title="", palette="Dark2", border.col=c("black", "grey", "grey"), border.lwds=c(1,0.5,0.1), fontsize.labels=c(0.7, 0.4, 0.3), fontcolor.labels=c("white", "white", "black"), fontface.labels=1, bg.labels=c("transparent"), align.labels=list( c("center", "center"), c("left", "top"), c("right", "bottom")), overlap.labels=0.5#, inflate.labels=T))


library(treemap)            
data <- read.table('/Volumes/seq_vgb/swofford/pathseq/test.test.treemap', header=F, sep=';')
colnames(data) <- c('Level1','Level2','Level3','Level4','Level5','Level6','Level7','count')
data$logCount <- log10(data$count)
p <- treemap(data, index=c('Level2','Level3','Level4','Level5'), vSize="logCount", type="index", title="", palette="Dark2", border.col=c("black", "grey", "grey"), border.lwds=c(1,0.5,0.1), fontsize.labels=c(0.7, 0.4, 0.3), fontcolor.labels=c("white", "white", "black"), fontface.labels=1, bg.labels=c("transparent"), align.labels=list( c("center", "center"), c("left", "top"), c("right", "bottom")), overlap.labels=0, inflate.labels=T)


#pathseqTo_igraph
library(data.tree)
library(treemap)
library(plyr)
library(igraph)
library(networkD3)
library(ggraph)
library(circlepackeR)
library(viridis)

nellie <- read.table('/Volumes/seq_vgb/swofford/pathseq/nellie/NelGood_sal_1.bam.pathseq.txt', header = T, stringsAsFactors = F, sep='\t', quote='', comment.char = '')
species <- nellie[nellie$type == 'species',]
species$pathString <- gsub('\\|','/',species$taxonomy)
species$pathString <- gsub('root/','',species$pathString)
species$pathString <- as.character(species$pathString)

keep <-!grepl('unclassified', species$pathString)
species <- species[keep,]
species <- species[species$unambiguous >10,]
species <- unique(species)


tree <- as.Node(species[,])

#aggregate unambigous reads up tree (optional depending on what we're doing)
tree$Do(function(x) { x$score <- Aggregate(node = x, attribute = "unambiguous", aggFun = sum)}, traversal = "post-order")

#optional pruning to decrease level complexity
Prune(tree, function(x) x$level < 5)

treeDf <-ToDataFrameNetwork(tree, 'name', 'score')
treeigraph <- graph_from_data_frame(treeDf)

edges <- treeDf[c('from','to')]
vertices <- treeDf[c('to', 'score', 'name')]
vertices <- rbind(data.frame(to = 'cellular_organisms', score = 0, name = 'Cellular Organisms'), vertices)
vertices$to <- gsub(".*/","", vertices$to)
edges$from <- gsub(".*/", "", edges$from)
edges$to <- gsub(".*/", "", edges$to)
vertices <- unique(vertices)

mygraph <- graph_from_data_frame( edges, vertices=vertices )

ggraph(mygraph, layout = 'circlepack', circular = TRUE, weight = 'score') +
geom_node_circle(aes(fill = as.factor(depth), color = as.factor(depth) )) +
scale_fill_manual(values=c("0" = "white", "1" = viridis(11)[1], "2" = viridis(11)[2], "3" = viridis(11)[3], "4"=viridis(11)[4], "5" = viridis(11)[5], "6" = viridis(11)[6], "7" = viridis(11)[7], "8" = viridis(11)[8], "9"=viridis(11)[9], "10" = viridis(11)[10], "11"=viridis(11)[11])) +
scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black", "5" = "black", "6" = "black", "7" = "black", "8"="black", "9" = "black", "10" = "black", "11" = "black") ) +
coord_fixed() + theme_void() + theme(legend.position = "FALSE")

ggraph(mygraph, 'circlepack', weight = 'score') + 
  geom_edge_link() + 
  geom_node_point(aes(colour = depth, size = 'score')) +
  coord_fixed()

#flower
ggraph(mygraph, layout = 'kk') + geom_node_circle(aes(r = log(score, base = 10000000), fill = as.factor(name), alpha = 0.5)) + geom_edge_link() + theme_void() + theme(legend.position = 'FALSE')

#interactive circlePack
circlepackeR(tree,'unambiguous')

#static circlePpack unweighted:
ggraph(mygraph, layout = 'circlepack', circular = TRUE) +
  geom_node_circle(aes(fill = as.factor(depth), color = as.factor(depth) )) +
  scale_fill_manual(values=c("0" = "white", "1" = viridis(11)[1], "2" = viridis(11)[2], "3" = viridis(11)[3], "4"=viridis(11)[4], "5" = viridis(11)[5], "6" = viridis(11)[6], "7" = viridis(11)[7], "8" = viridis(11)[8], "9"=viridis(11)[9], "10" = viridis(11)[10], "11"=viridis(11)[11])) +
  scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black", "5" = "black", "6" = "black", "7" = "black", "8"="black", "9" = "black", "10" = "black", "11" = "black") ) +
  coord_fixed() + theme_void() + theme(legend.position = "FALSE")

#static Sunburst
ggraph(mygraph, layout = 'partition', circular = TRUE, weight = 'score') + geom_node_arc_bar(aes(fill = name), size = 0.25) + coord_fixed() + theme_void() + theme(legend.position = 'FALSE')

ggraph(mygraph, 'circlepack', weight = 'score') + geom_edge_link() + geom_node_point(aes(colour = depth, size = 'score')) + coord_fixed()

#dendrogram
ggraph(mygraph, 'dendrogram') + geom_edge_diagonal(aes(color = as.factor(from))) + theme_void() + theme(legend.position = 'FALSE')
#dendrogram circular
ggraph(mygraph, 'dendrogram', circular = TRUE) + geom_edge_diagonal(aes(color = as.factor(from))) + theme_void() + theme(legend.position = 'FALSE') + scale_color_viridis_d()
ggraph(mygraph, 'dendrogram', circular = TRUE) + geom_edge_diagonal(color = viridis(59200)) + theme_void() + theme(legend.position = 'FALSE')

#sankey
topLevel <- data.frame(from = c(NA, 'total', 'total','total'), to = c('total','dog','ambiguous','cellular_organisms'), name = c('total','nellie','ambiguous','cellular_organisms'), score = c(nellieTotal,nellieDog, (nellieTotal - nellieDog - sum(treeDf[treeDf$from == 'cellular_organisms','score'])), sum(treeDf[treeDf$from == 'cellular_organisms','score'])))
treeDf <- rbind(topLevel, treeDf)

links <- data.frame(source = treeDf$from, target = treeDf$to, value = treeDf$score, group = treeDf$to)
nodes <- data.frame(node = c(0:(length(unique(treeDf$to))-1)), name = unique(treeDf$to))
results <- merge(links,nodes, by.x = 'source', by.y = 'name')
results <- merge(results,nodes, by.x = 'target', by.y = 'name')
links <- results[c('node.x','node.y','value', 'group')]
nodes$name <- gsub('.*/','', nodes$name)
colnames(links) <- c('source','target','value', 'group')
sankeyNetwork(Links = links, Nodes = nodes, Source = 'source',Target = 'target', Value = 'value', NodeID = 'name', fontSize = 12, nodeWidth = 30, nodePadding = 20, LinkGroup = 'group')

#donut
forDonut <- nellie[nellie$type == 'genus',]
forDonut <- forDonut[forDonut$unambiguous > 100,]



#comboBarChart
forDonut <- rbind(nellie[nellie$type == 'superkingdom',], nellie[nellie$type == 'kingdom',])
forDonut <- forDonut[c('name','reads')]
sumMicro <- data.frame(Source = 'Microbiome', SequencingReads = sum(forDonut$reads))
topLevelDonut <- topLevel[c('name', 'score')]
topLevelDonut <- topLevelDonut[topLevelDonut$name != 'total',]
topLevelDonut <- topLevelDonut[topLevelDonut$name != 'cellular_organisms',]
colnames(topLevelDonut) <- c('Source','SequencingReads')
forDonut <- rbind(sumMicro,topLevelDonut)
forDonut$Source <- gsub('nellie','Pirl',forDonut$Source)
forDonut$Source <- gsub('ambiguous', 'Unknown', forDonut$Source)
forDonut$SequencingReads <- factor(forDonut$SequencingReads, unique(forDonut$SequencingReads))




ggplot(forDonut, aes(x=2, y=SequencingReads, fill = Source)) + geom_bar(stat = 'identity') + theme_bw() + theme_void() + scale_fill_viridis_d()



colorsTest <- c("#683591", "#0000ff","#00bcff","#0083b2","#0000ff","#264f55")
colorSet <- colorRampPalette(colorsTest)(59200)
p1 <- ggplot(forDonut, aes(x=1, y=SequencingReads, fill = Source)) + geom_bar(stat = 'identity') + theme_bw() + theme_void() + scale_fill_manual(values = c('#683591','#b8d432','#ed4096'))+ theme(legend.position = 'left')
p3 <- ggplot(forDonut, aes(x=1, y=SequencingReads, fill = Source)) + geom_bar(stat = 'identity') + theme_bw() + theme_void() + scale_fill_manual(values = c('#264f55','#264f55','#264f55'))+ theme(legend.position = 'FALSE')
p2 <- ggraph(mygraph, 'dendrogram', circular = TRUE) + geom_edge_diagonal(color = colorSet) + theme_void() + theme(legend.position = 'FALSE') + labs(caption = 'Each Leaf, a Species in Pirl\'s Microbiome') + theme(plot.caption = element_text(hjust=0.5, size=rel(1.2)))

lay <- rbind(c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,3,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
             c(1,1,1,1,1,1,1,1,1,1,1,1,NA,NA,NA,NA,3,NA,NA,NA,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
             c(1,1,1,1,1,1,1,1,1,1,1,1,NA,NA,NA,NA,3,NA,NA,NA,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
             c(1,1,1,1,1,1,1,1,1,1,1,1,NA,NA,NA,NA,3,NA,NA,NA,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
             c(1,1,1,1,1,1,1,1,1,1,1,1,NA,NA,NA,NA,3,NA,NA,NA,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
             c(1,1,1,1,1,1,1,1,1,1,1,1,NA,NA,NA,NA,3,NA,NA,NA,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
             c(1,1,1,1,1,1,1,1,1,1,1,1,NA,NA,NA,NA,3,NA,NA,NA,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
             c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,3,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
grid.arrange(p1,p2,p3, layout_matrix=lay,top = 'Pirl\'s Mouth Hosts a Diverse and Vibrant Community')
p1Grob <- ggplotGrob(p1)
p2Grob <- ggplotGrob(p2)

microPercentage <- forDonut[forDonut$Source == 'Microbiome','SequencingReads']/sum(forDonut$SequencingReads)
arrows <- data.frame(x = 0, xend = 1, x1 = 0.285, x2 = 0.285, x1end = 0.6, y = 0, yend = 1, y1 = 1-(microPercentage/2), y1end = 1, y2 = 1-(microPercentage/2), y2end = .01, x2end = 0.6)
#arrows <- data.frame(x = 0, xend = 1, x1 = 0.28, x2 = 0.28, x1end = 1.1, y = 0, yend = 1, y1 = 1-(microPercentage/2), y1end = .5, y2 = 1-(microPercentage/2), y2end = .5, x2end = 1.1)

ggplot(arrows) + geom_curve(data = arrows, aes(x=x, xend=xend, y=y, yend=yend), color = 'white', alpha = 0) + theme_void() + annotation_custom(p1Grob, xmin = 0, xmax = 0.3) + annotation_custom(p2Grob, xmin = 0.3) + geom_curve(data = arrows, aes(x=x1, xend = x1end, y = y1, yend = y1end), curvature = -0.2, color = "#683591", size = 3, lineend = 'round') + geom_curve(data = arrows, aes(x=x2, xend=x2end, y=y2, yend=y2end), color = "#683591", size = 3, lineend = 'round')

             