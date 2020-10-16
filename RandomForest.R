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
geneData.noBatchEffect <- as.data.frame(t(removeBatchEffect(t(geneData.scaled[,1:22277]), geneData$batch)))

geneData.noBatchEffect$virus <- virus
geneData.noBatchEffect$batch <- batch

#Create vector sick vs healty (without considering what type of virus)
sick.healthy <- rep(0, length(virus))
sick.healthy[virus == "healthy"] <- "healthy"
sick.healthy[virus == "infA" | virus == "RSV" | virus == "rhino"] <- "sick"

#add vector to dataframe
geneData.noBatchEffect$sick.healthy <- as.factor(sick.healthy)

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
Y_train <- X_train$sick.healthy
Y_test <- X_test$sick.healthy
#Remove Virus and Batch from Train and Test subsets
X_train$virus <- NULL
X_train$batch <- NULL
X_train$sick.healthy <- NULL
X_test$virus <- NULL
X_test$batch <- NULL
X_test$sick.healthy <- NULL


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

# ---------------------------------------------------------------------------------------------
# RANDOM FOREST FOR ENTIRE DATASET
# ---------------------------------------------------------------------------------------------

library(randomForest)
library(rfUtilities)
library(caret)

geneData_rf <- randomForest(x = X_train, y = Y_train, ntree=1000)

geneData_rf.pred.test <- predict(geneData_rf, X_test)
confusionMatrix(geneData_rf.pred.test, Y_test)
#poor classification, same accuracy of No-Information-Rate

#cross-validation #takes over 15 minutes
#control <- trainControl(method = "cv", number = 5)
#geneData.cv <- train(X_train, Y_train, method = "rf", metric = "Accuracy", trControl = control)
#geneData.cv_pred <- predict(geneData.cv, X_test)
#confusionMatrix(geneData.cv_pred, Y_test)

# ---------------------------------------------------------------------------------------------
# RANDOM FOREST FOR RHINO
# ---------------------------------------------------------------------------------------------

rhino_rf <- randomForest(x = X_rhino_train, y = Y_rhino_train, ntree=1000)

rhino_rf.pred.test <- predict(rhino_rf, X_rhino_test)
confusionMatrix(rhino_rf.pred.test, Y_rhino_test) #70.00%

#cross-validation
control <- trainControl(method = "cv", number = 5)
rhino.cv <- train(X_rhino_train, Y_rhino_train, method = "rf", metric = "Accuracy", trControl = control)
rhino.cv_pred <- predict(rhino.cv, X_rhino_test)
confusionMatrix(rhino.cv_pred, Y_rhino_test) #90.00%

#ploting gene importance and saving most important genes
plot(sort(rhino_rf$importance, decreasing = TRUE)[1:4000])
prob.names.rhino <- rownames(rhino_rf$importance)
top500_rhino <- prob.names.rhino[order(rhino_rf$importance, decreasing = TRUE)][1:500]
#save file
#write.csv(top500_rhino, file = "top500_rhino.txt", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------------------------
# RANDOM FOREST FOR RSV
# ---------------------------------------------------------------------------------------------

RSV_rf <- randomForest(x = X_RSV_train, y = Y_RSV_train, ntree=1000)

RSV_rf.pred.test <- predict(RSV_rf, X_RSV_test)
confusionMatrix(RSV_rf.pred.test, Y_RSV_test) #90.00%

#cross-validation
control <- trainControl(method = "cv", number = 5)
RSV.cv <- train(X_RSV_train, Y_RSV_train, method = "rf", metric = "Accuracy", trControl = control)
RSV.cv_pred <- predict(RSV.cv, X_RSV_test)
confusionMatrix(RSV.cv_pred, Y_RSV_test) #90.00%

#ploting gene importance and saving most important genes
plot(sort(RSV_rf$importance, decreasing = TRUE)[1:4000])
prob.names.RSV <- rownames(RSV_rf$importance)
top500_RSV <- prob.names.RSV[order(RSV_rf$importance, decreasing = TRUE)][1:500]
#save file
#write.csv(top500_RSV, file = "top500_RSV.txt", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------------------------
# RANDOM FOREST FOR INFA
# ---------------------------------------------------------------------------------------------

infA_rf <- randomForest(x = X_infA_train, y = Y_infA_train, ntree=1000)

infA_rf.pred.test <- predict(infA_rf, X_infA_test)
confusionMatrix(infA_rf.pred.test, Y_infA_test) #77.78%

#cross-validation
control <- trainControl(method = "cv", number = 5)
infA.cv <- train(X_infA_train, Y_infA_train, method = "rf", metric = "Accuracy", trControl = control)
infA.cv_pred <- predict(infA.cv, X_infA_test)
confusionMatrix(infA.cv_pred, Y_infA_test) #77.78%

#ploting gene importance and saving most important genes
plot(sort(infA_rf$importance, decreasing = TRUE)[1:4000])
prob.names.infA <- rownames(infA_rf$importance)
top500_infA <- prob.names.infA[order(infA_rf$importance, decreasing = TRUE)][1:500]
#save file
#write.csv(top500_infA, file = "top500_infA.txt", quote = FALSE, row.names = FALSE)
