library(ggplot2)
library(viridis)
library(reshape2)
library(plyr)
library(dplyr)

confusionMatrix <- read.table('/Volumes/seq_vgb/swofford/R01/convertingToSKLInput/testVCFs/confusionMatrix.out', quote='', header = T, sep ='\t', stringsAsFactors = F)
confusionMatrix$Actual_Label <- gsub('angioSarcoma','AS',confusionMatrix$Actual_Label)
colnames(confusionMatrix) <- gsub('angioSarcoma','AS',colnames(confusionMatrix))
confusionMatrix <- confusionMatrix[, order(colnames(confusionMatrix))]
confusionMatrix <- confusionMatrix[order(confusionMatrix$Actual_Label),]

test <- confusionMatrix[,-1]/rowSums(confusionMatrix[,-1])
test$Actual <- confusionMatrix$Actual_Label
testMelt <- melt(test)
ggplot(testMelt, aes(x=Actual_Label, variable, fill = value, colour = 'black', title('test'))) + geom_tile() + scale_fill_viridis() + geom_text(aes(label = paste(round(100*value),'%',sep=''))) + theme(legend.position = 'FALSE') + ylab('Prediction') + ggtitle('Percentage of Actual Types Predicted in Each Class')
ggplot(testMelt, aes(x=Actual, variable, fill = value, colour = 'black', title('test'))) + geom_tile() + scale_fill_viridis() + geom_text(aes(label = round(100*value))) + theme(legend.position = 'FALSE') + ylab('Prediction') + ggtitle('Percentage of Actual Types Predicted in Each Class') + theme(axis.text.x = element_text(angle = 90))

meltedConfusion <- melt(confusionMatrix)
confusionAggregate <- aggregate(. ~Actual_Label,meltedConfusion,sum)
testConfusionMatrixTotals <- confusionMatrix
testConfusionMatrixTotals$totals <- confusionAggregate$value
test <- confusionMatrix[,-1]/confusionAggregate$value
test$Actual <- confusionMatrix$Actual_Label
testMelt <- melt(test)



confusionAdjust <- meltedConfusion[,-1]/confusionAggregate$value
confusionAdjust$Actual <- confusionMatrix$Actual_Label
confusionAdjustMelt <- melt(confusionAdjust)

ggplot(confusionAdjustMelt, aes(x=Actual, variable, fill = value, colour = 'black', title('test'))) + geom_tile() + scale_fill_viridis() + geom_text(aes(label = paste(round(100*value),'%',sep=''))) + theme(legend.position = 'FALSE') + ylab('Prediction') + ggtitle('Percentage of Actual Types Predicted in Each Class')
ggplot(meltedConfusion, aes(x=Actual_label, variable, fill = value, colour = 'black', title('test'))) + geom_tile() + scale_fill_viridis() + geom_text(aes(label = value)) + theme(legend.position = 'FALSE') + ggtitle('Predicted vs. Actual Sample Labels') +ylab('Predicted') + xlab('Actual')
       

confusionMatrix <- read.table('/Volumes/seq_vgb/swofford/R01/convertingToSKLInput/testVCFs/dogPreds042621.txt', quote='', header = T, sep ='\t', stringsAsFactors = F)
confusionMatrix$Actual_label <- gsub('dogBCL','dogLSA', confusionMatrix$Actual_label)
confusionMatrix = confusionMatrix[confusionMatrix$Actual_label != 'dog_MCT_Carcinoma',]
confusionMatrix = confusionMatrix[confusionMatrix$Actual_label != 'dog_MCT_Benign',]
confusionMatrix = aggregate(. ~ Actual_label, data = confusionMatrix, FUN=sum)
meltedConfusion <- melt(confusionMatrix)

test <- confusionMatrix[,-1]/rowSums(confusionMatrix[,-1])
ggplot(meltedConfusion, aes(x=Actual_label, variable, fill = value, colour = 'black', title('test'))) + geom_tile() + scale_fill_viridis() + geom_text(aes(label = value), color = "white") + theme(legend.position = 'FALSE') + ggtitle('Predicted vs. Actual Sample Labels') +ylab('Predicted') + xlab('Actual')
test$Actual_label = confusionMatrix$Actual_label
meltedConfusion <- melt(test)
ggplot(meltedConfusion, aes(x=Actual_label, variable, fill = value, colour = 'black', title('test'))) + geom_tile() + scale_fill_viridis() + geom_text(aes(label = paste(round(100*value),'%',sep='')), color="white") + theme(legend.position = 'FALSE') + ylab('Prediction') + ggtitle('Percentage of Actual Types Predicted in Each Class')

