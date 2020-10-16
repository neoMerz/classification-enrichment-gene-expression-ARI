wd <- getwd()
setwd(wd)

# Read files
geneData <- read.csv(file = "geneData.csv", sep = ",", row.names = 1, header = TRUE)

# Normalization - Log transformation and Scaling
#geneData.log <-log(geneData[,1:22277])
#geneData.scaled <- as.data.frame(scale(geneData.log))

# Normalization - Only Scaling
geneData.scaled <- as.data.frame(scale(geneData[,1:22277]))

# Append Group and Batch to newly created dataframe
virus <- factor(geneData[,22278])
batch <- factor(geneData[,22279])
geneData.scaled$virus <- virus
geneData.scaled$batch <- batch

# Create subset dataframes based on batches (Rhino / RSV / infA)
rhino <- droplevels(geneData[geneData$batch == "A",])
RSV <- droplevels(geneData[geneData$batch == "B",])
infA <- droplevels(geneData[geneData$batch == "C",])

#Remove Batch Effect 
library("limma")
geneData.noBatchEffect <- as.data.frame(t(removeBatchEffect(t(geneData[,1:22277]), geneData$batch)))

geneData.noBatchEffect$virus <- virus
geneData.noBatchEffect$batch <- batch

# ---------------------------------------------------------------------------------------------
# PCA AND K-MEANS CLUSTERING FOR ENTIRE DATASET (after removing batch effect)
# ---------------------------------------------------------------------------------------------

# PCA
pca <- prcomp(geneData.noBatchEffect[,1:22277])
summary(pca)
screeplot(pca)

plot(pca$x[,1], pca$x[,2], xlab = "PCA1", ylab = "PCA2", main = "PCA for components 1&2", type = "p", pch=10)

# K-MEANS 
library("useful")
k <- 4 # 3 virus + control group
fullDS_kmeans_result <- kmeans(geneData.scaled.noBatchEffect[,1:22277], k)
table(fullDS_kmeans_result$cluster)
#plot k-means cluster with respective labels
plot(fullDS_kmeans_result, data=geneData.scaled.noBatchEffect[,1:22277]) + geom_text(aes(label=virus), hjust=0, vjust=0)

#Even though removeBatchEffect has done a good job of "cleaning", unsupervised classification is not able to distinguish
#between the three virus and the control group.

# ---------------------------------------------------------------------------------------------
# PCA AND K-MEANS CLUSTERING FOR EACH BATCH
# ---------------------------------------------------------------------------------------------

# PCA RHINO
rhinoPCA <- prcomp(rhino[,1:22277])
summary(rhinoPCA)
screeplot(rhinoPCA)

plot(rhinoPCA$x[,1], rhinoPCA$x[,2], xlab = "PCA1", ylab = "PCA2", main = "Rhino PCA for components 1&2", type = "p", pch=10)

# K-MEANS RHINO
k <- 2
rhino_kmeans_result <- kmeans(rhino[,1:22277], k)
table(rhino_kmeans_result$cluster)
#plot k-means cluster with respective labels
plot(rhino_kmeans_result, data=rhino[,1:22277]) + geom_text(aes(label=rhino$virus), hjust=0, vjust=0)

#Unsupervised classification is not good

# ---------------------------------------------------------------------------------------------

# PCA RSV 
RSVPCA <- prcomp(RSV[,1:22277])
summary(RSVPCA)
screeplot(RSVPCA)

plot(RSVPCA$x[,1], RSVPCA$x[,2], xlab = "PCA1", ylab = "PCA2", main = "RSV PCA for components 1&2", type = "p", pch=10)

#K-MEANS RSV
k <- 2
RSV_kmeans_result <- kmeans(RSV[,1:22277], k)
table(RSV_kmeans_result$cluster)
#plot k-means cluster with respective labels
plot(RSV_kmeans_result, data=RSV) + geom_text(aes(label=RSV$virus), hjust=0, vjust=0)

#Unsupervised classification is not good

# ---------------------------------------------------------------------------------------------

# PCA infA
infAPCA <- prcomp(infA[,1:22277])
summary(infAPCA)
screeplot(infAPCA)

plot(infAPCA$x[,1], infA$x[,2], xlab = "PCA1", ylab = "PCA2", main = "infA PCA for components 1&2", type = "p", pch=10)

#K-MEANS infA
k <- 2
infA_kmeans_result <- kmeans(infA[,1:22277], k)
table(infA_kmeans_result$cluster)
#plot k-means cluster with respective labels
plot(infA_kmeans_result, data=infA) + geom_text(aes(label=infA$virus), hjust=0, vjust=0)

#Unsupervised classification is decent

# ---------------------------------------------------------------------------------------------

