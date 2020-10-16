wd <- getwd()
setwd(wd)

# Read files
geneData <- read.csv(file = "geneData.csv", sep = ",", row.names = 1, header = TRUE)

# Normalization - Log transformation and Scaling
#geneData.log <-log(geneData[,1:22277])
#geneData.scaled <- as.data.frame(scale(geneData.log))

# Normalization - Only Scaling
#geneData.scaled <- as.data.frame(scale(geneData[,1:22277]))

# Append Group and Batch to newly created dataframe
virus <- factor(geneData[,22278])
batch <- factor(geneData[,22279])

# Create subset dataframes based on batches (Rhino / RSV / infA)
rhino <- droplevels(geneData[geneData$batch == "A",])
RSV <- droplevels(geneData[geneData$batch == "B",])
infA <- droplevels(geneData[geneData$batch == "C",])


#Remove Batch Effect 
library("limma")
geneData.noBatchEffect <- as.data.frame(t(removeBatchEffect(t(geneData[,1:22277]), geneData$batch)))

geneData.noBatchEffect$virus <- virus
geneData.noBatchEffect$batch <- batch

#Create vector sick vs healty (without considering what type of virus)
sick.healthy <- rep(0, length(virus))
sick.healthy[virus == "healthy"] <- "healthy"
sick.healthy[virus == "infA" | virus == "RSV" | virus == "rhino"] <- "sick"

#add vector to dataframe
geneData.noBatchEffect$sick.healthy <- sick.healthy


# ---------------------------------------------------------------------------------------------
# TRAIN AND TEST SPLIT (for enitre dataset and also for each batch)
# ---------------------------------------------------------------------------------------------

# ENTIRE DATASET
set.seed(23)
smp_size <- floor(0.75* nrow(geneData.noBatchEffect))
train_ind <- sample(seq_len(nrow(geneData.noBatchEffect)), size = smp_size)
X_train <- geneData.noBatchEffect[train_ind,]
X_test <- geneData.noBatchEffect[-train_ind,]
#Save labels in vectors and remove them from the datasets
geneData.noBatchEffect$virus <- NULL
geneData.noBatchEffect$batch <- NULL
#Train and Test target labels
Y_train <- as.factor(X_train$sick.healthy)
Y_test <- as.factor(X_test$sick.healthy)
#Remove Virus and Batch from Train and Test subsets
X_train$virus <- NULL
X_train$batch <- NULL
X_train$sick.healthy <- NULL
X_test$virus <- NULL
X_test$batch <- NULL
X_test$sick.healthy <- NULL
#Transpose Train and test set so that they can be used with SCUDO
X_train <- as.data.frame(t(X_train))
X_test <- as.data.frame(t(X_test))



# RHINO
set.seed(123)
smp_size <- floor(0.75* nrow(rhino))
train_ind <- sample(seq_len(nrow(rhino)), size = smp_size)
X_rhino_train <- rhino[train_ind,]
X_rhino_test <- rhino[-train_ind,]
#Train and Test target labels
Y_rhino_train <- factor(X_rhino_train$virus)
Y_rhino_test <- factor(X_rhino_test$virus)
#Remove Virus and Batch from Train and Test subsets
X_rhino_train$virus <- NULL
X_rhino_train$batch <- NULL
X_rhino_test$virus <- NULL
X_rhino_test$batch <- NULL
#Transpose Train and test set so that they can be used with SCUDO
X_rhino_train <- as.data.frame(t(X_rhino_train))
X_rhino_test <- as.data.frame(t(X_rhino_test))



#RSV
set.seed(123)
smp_size <- floor(0.75* nrow(RSV))
train_ind <- sample(seq_len(nrow(RSV)), size = smp_size)
X_RSV_train <- RSV[train_ind,]
X_RSV_test <- RSV[-train_ind,]
#Train and Test target labels
Y_RSV_train <- factor(X_RSV_train$virus)
Y_RSV_test <- factor(X_RSV_test$virus)
#Remove Virus and Batch from Train and Test subsets
X_RSV_train$virus <- NULL
X_RSV_train$batch <- NULL
X_RSV_test$virus <- NULL
X_RSV_test$batch <- NULL
#Transpose Train and test set so that they can be used with SCUDO
X_RSV_train <- as.data.frame(t(X_RSV_train))
X_RSV_test <- as.data.frame(t(X_RSV_test))



#INFA
set.seed(123)
smp_size <- floor(0.75* nrow(infA))
train_ind <- sample(seq_len(nrow(infA)), size = smp_size)
X_infA_train <- infA[train_ind,]
X_infA_test <- infA[-train_ind,]
#Train and Test target labels
Y_infA_train <- factor(X_infA_train$virus)
Y_infA_test <- factor(X_infA_test$virus)
#Remove Virus and Batch from Train and Test subsets
X_infA_train$virus <- NULL
X_infA_train$batch <- NULL
X_infA_test$virus <- NULL
X_infA_test$batch <- NULL
#Transpose Train and test set so that they can be used with SCUDO
X_infA_train <- as.data.frame(t(X_infA_train))
X_infA_test <- as.data.frame(t(X_infA_test))


# ---------------------------------------------------------------------------------------------
# SCUDO for entire dataset
# ---------------------------------------------------------------------------------------------

library("rScudo")
library("caret")

#analyze training set
trainRes <- scudoTrain(X_train, Y_train, nTop = 30, nBottom = 30, alpha = 0.1)
trainRes

#inspect signatures
upSignatures(trainRes)[1:5,1:5]
consensusUpSignatures(trainRes)[1:5,]
consensusDownSignatures(trainRes)[1:5,]

#generate and plot map of training sample
trainNet <- scudoNetwork(trainRes, N = 0.25)
scudoPlot(trainNet, vertex.label = NA)

#perform validation on test set
testRes <- scudoTest(trainRes, X_test, Y_test, nTop = 30, nBottom = 30)

testNet <- scudoNetwork(testRes, N = 0.4) 
scudoPlot(testNet, vertex.label = NA)

#identify clusters on map
library("igraph")
testClust <- igraph::cluster_spinglass(testNet, spins = 4)
plot(testClust, testNet, vertex.label = NA)
  
#perform classification
classRes <- scudoClassify(X_train, X_test, N=0.25, nTop = 15, nBottom = 15, trainGroups = Y_train, alpha = 0.1)
confusionMatrix(classRes$predicted, Y_test)

# ---------------------------------------------------------------------------------------------
# SCUDO for RHINO
# ---------------------------------------------------------------------------------------------
  
#analyze training set
trainRes_rhino <- scudoTrain(X_rhino_train, Y_rhino_train, nTop = 30, nBottom = 30, alpha = 0.1)
trainRes_rhino

#inspect signatures
upSignatures(trainRes_rhino)[1:5,1:5]
consensusUpSignatures(trainRes_rhino)[1:5,]
consensusDownSignatures(trainRes_rhino)[1:5,]

#generate and plot map of training sample
trainNet_rhino <- scudoNetwork(trainRes_rhino, N = 0.25)
scudoPlot(trainNet_rhino, vertex.label = NA)

#perform validation on test set
testRes_rhino <- scudoTest(trainRes_rhino, X_rhino_test, Y_rhino_test, nTop = 30, nBottom = 30)
testNet_rhino <- scudoNetwork(testRes_rhino, N = 0.4) 
scudoPlot(testNet_rhino, vertex.label = NA)

#identify clusters on map
library("igraph")
testClust_rhino <- igraph::cluster_spinglass(testNet_rhino, spins = 2)
plot(testClust_rhino, testNet_rhino, vertex.label = NA)

#perform classification
classRes_rhino <- scudoClassify(X_rhino_train, X_rhino_test, N=0.25, nTop = 30, nBottom = 30, trainGroups = Y_rhino_train, alpha = 0.1)
confusionMatrix(classRes_rhino$predicted, Y_rhino_test)

# ---------------------------------------------------------------------------------------------
# SCUDO for RSV
# ---------------------------------------------------------------------------------------------

#analyze training set
trainRes_RSV <- scudoTrain(X_RSV_train, Y_RSV_train, nTop = 30, nBottom = 30, alpha = 0.1)
trainRes_RSV

#inspect signatures
upSignatures(trainRes_RSV)[1:5,1:5]
consensusUpSignatures(trainRes_RSV)[1:5,]
consensusDownSignatures(trainRes_RSV)[1:5,]

#generate and plot map of training sample
trainNet_RSV <- scudoNetwork(trainRes_RSV, N = 0.25)
scudoPlot(trainNet_RSV, vertex.label = NA)

#perform validation on test set
testRes_RSV <- scudoTest(trainRes_RSV, X_RSV_test, Y_RSV_test, nTop = 30, nBottom = 30)
testNet_RSV <- scudoNetwork(testRes_RSV, N = 0.4) 
scudoPlot(testNet_RSV, vertex.label = NA)

#identify clusters on map
library("igraph")
testClust_RSV <- igraph::cluster_spinglass(testNet_RSV, spins = 2)
plot(testClust_RSV, testNet_RSV, vertex.label = NA)

#perform classification
classRes_RSV <- scudoClassify(X_RSV_train, X_RSV_test, N=0.25, nTop = 30, nBottom = 30, trainGroups = Y_RSV_train, alpha = 0.1)
confusionMatrix(classRes_RSV$predicted, Y_RSV_test)

# ---------------------------------------------------------------------------------------------
# SCUDO for INFA
# ---------------------------------------------------------------------------------------------

#analyze training set
trainRes_infA <- scudoTrain(X_infA_train, Y_infA_train, nTop = 30, nBottom = 30, alpha = 0.1)
trainRes_infA

#inspect signatures
upSignatures(trainRes_infA)[1:5,1:5]
consensusUpSignatures(trainRes_infA)[1:5,]
consensusDownSignatures(trainRes_infA)[1:5,]

#generate and plot map of training sample
trainNet_infA <- scudoNetwork(trainRes_infA, N = 0.25)
scudoPlot(trainNet_infA, vertex.label = NA)

#perform validation on test set
testRes_infA <- scudoTest(trainRes_infA, X_infA_test, Y_infA_test, nTop = 80, nBottom = 80)
testNet_infA <- scudoNetwork(testRes_infA, N = 0.4) 
scudoPlot(testNet_infA, vertex.label = NA)

#identify clusters on map
library("igraph")
testClust_infA <- igraph::cluster_spinglass(testNet_infA, spins = 2)
plot(testClust_infA, testNet_infA, vertex.label = NA)

#perform classification
classRes_infA <- scudoClassify(X_infA_train, X_infA_test, N=0.25, nTop = 80, nBottom = 80, trainGroups = Y_infA_train, alpha = 0.1)
confusionMatrix(classRes_infA$predicted, Y_infA_test)
