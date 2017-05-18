#libraries
library(ggplot2)


setwd("~/Dropbox/McMaster/PhD/RScripts/pca/") #working directory goes here
fpkm <- read.csv("~/Documents/readCounts/fpkm/20161230-fpkm.csv", row.names=1)

drought <- fpkm[,c(13:28,44:45)]

drought.log <- log2(drought+1)


## all genes
pca<- prcomp(drought.log, center=TRUE, scale=TRUE)  # PCA with centering and scaling
pca$rotation  # The loadings are here


## more variable genes
rv <- rowVars(as.matrix(drought.log))
names(rv) <- rownames(drought.log)
select <- order(rv, decreasing = TRUE)[1:2000]
pca.top<- prcomp(drought.log[select,], center=TRUE, scale=TRUE)
pca.top$rotation # very similar....so maybe ok, but different explained variance
##########

plot(pca, type = "l")
summary(pca)
summary(pca.top)

########
# FIGS #
########

exprVals<-data.frame(pca$x)
sampleVals<-data.frame(pca$rotation)


dim(exprVals) 
dim(sampleVals)


samples <- substring(colnames(drought)[1:16],1,3)
samples <- c(samples, rep("Y-F2003", 2))
ecotype <- c(rep("Shandong", 8), rep("Yukon", 10))
cond <- substring(colnames(drought)[1:16],2,3)
cond <- c(cond, rep("field", 2))
label <- c(rep("no", 9), "yes", "no", "yes", "no", "no", "yes", "no", "no", "no" )


## extrating all PCA data for ggplot
coords<-data.frame(X=rep(0, 18), Y=rep(0, 18), sampleVals, Samples = samples, ecotype = ecotype, condition = cond, label=label, name = colnames(drought))

### Figures ###
# PC1 vs PC2 #
# PC1 vs PC2 (all samples)



theme_set(theme_bw()) # make background white
#pc12plot<-ggplot(exprVals, aes(PC1, PC2)) + geom_point(shape=19, alpha=0.2) + geom_segment(data=coords, mapping=aes(x=X, y=Y, xend=X+PC1, yend=Y+PC2, colour=Samples), arrow=arrow(), size=1) 

# Individual plots
pc12plot <- ggplot(coords, aes(x = PC1, y = PC2)) + # data that you want to plot
  geom_point(aes(color = condition, shape = ecotype), size=5) + # adding PCA loadings as a "Scatter plot"
  coord_cartesian(xlim=c(-0.5,0.5), ylim=c(-0.5,0.5)) + # "zooming" in on the plot
  geom_text(data=subset(coords, label == "yes"), 
            aes(PC3,PC4,label=name)) # adding text to a subset of the data

pc34plot <- ggplot(coords, aes(x = PC3, y = PC4)) +
  geom_point(aes(color = condition, shape = ecotype), size=5) + 
  coord_cartesian(xlim=c(-0.5,0.5), ylim=c(-0.5,0.5)) + 
  geom_text(data=subset(coords, label == "yes"),
            aes(PC3,PC4,label=name))

## loop for all PC1 combinations
for (tryThis in 2:10){
  toUse<-noquote(paste("PC", tryThis, sep=""))
  fileName<-paste("drought_figs/PC1",toUse, ".pdf",sep="")
  pdf(fileName)
  print(pc12plot <- ggplot(coords, aes_string("PC1", toUse)) +
          geom_point(aes(color = condition, shape = ecotype), size=5) + 
          coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8)) + 
          geom_text(data=subset(coords, label == "yes"), aes(label=name)))
  dev.off()
}

## PC2
for (tryThis in 3:10){
  toUse<-noquote(paste("PC", tryThis, sep=""))
  fileName<-paste("drought_figs/PC2",toUse, ".pdf",sep="")
  pdf(fileName)
  print(pc12plot <- ggplot(coords, aes_string("PC2", toUse)) +
          geom_point(aes(color = condition, shape = ecotype), size=5) + 
          coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8)) + 
          geom_text(data=subset(coords, label == "yes"), aes(label=name)))
  dev.off()
}

## PC3
for (tryThis in 4:10){
  toUse<-noquote(paste("PC", tryThis, sep=""))
  fileName<-paste("drought_figs/PC3",toUse, ".pdf",sep="")
  pdf(fileName)
  print(pc12plot <- ggplot(coords, aes_string("PC3", toUse)) +
          geom_point(aes(color = condition, shape = ecotype), size=5) + 
          coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8)) + 
          geom_text(data=subset(coords, label == "yes"), aes(label=name)))
  dev.off()
}

for (tryThis in 5:10){
  toUse<-noquote(paste("PC", tryThis, sep=""))
  fileName<-paste("drought_figs/PC4",toUse, ".pdf",sep="")
  pdf(fileName)
  print(pc12plot <- ggplot(coords, aes_string("PC4", toUse)) +
          geom_point(aes(color = condition, shape = ecotype), size=5) + 
          coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8)) + 
          geom_text(data=subset(coords, label == "yes"), aes(label=name)))
  dev.off()
}

for (tryThis in 6:10){
  toUse<-noquote(paste("PC", tryThis, sep=""))
  fileName<-paste("drought_figs/PC5",toUse, ".pdf",sep="")
  pdf(fileName)
  print(pc12plot <- ggplot(coords, aes_string("PC5", toUse)) +
          geom_point(aes(color = condition, shape = ecotype), size=5) + 
          coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.8,0.8)) + 
          geom_text(data=subset(coords, label == "yes"), aes(label=name)))
  dev.off()
}

############################################
# Get list of genes positive PC2, neg PC10 #
# Get list of genes positive PC2, pos PC10 #
############################################

### Getting negative and positive genes for PC2 ###

posPC2genes <- which(exprVals[,2]>0)
posPC2 <- data.frame(exprVals[posPC2genes,c(2,10)])

rownames(posPC2) <- rownames(exprVals[posPC2genes,])

# negative pc10
pos2neg10genes <- which(posPC2$PC10 < 0)
pos2neg10 <- data.frame(posPC2[pos2neg10genes,])


# positive pc10
pos2pos10genes <- which(posPC2$PC10 > 0)
pos2pos10 <- data.frame(posPC2[pos2pos10genes,])

expr_annot <- read.csv("~/Documents/readCounts/fpkm/2016-12-30_fpkm_annot.csv", row.names=1)

neg10genes <- match(rownames(pos2neg10), rownames(expr_annot))
pos10genes <- match(rownames(pos2pos10), rownames(expr_annot))


neg10 <- cbind(pos2neg10, expr_annot[neg10genes,7:89])

pos10 <- cbind(pos2pos10, expr_annot[pos10genes,7:89])

write.csv(neg10, file="drought_genes/pos2neg10.csv")
write.csv(pos10, file="drought_genes/pos2pos10.csv")
