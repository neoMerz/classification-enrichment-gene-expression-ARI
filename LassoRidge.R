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
# Lasso AND RIDGE REGRESSION FOR ENTIRE DATASET
# ---------------------------------------------------------------------------------------------

library("glmnet")
library("caret")

control <- trainControl(method = "cv", number = 10)
cv.fit.all <- train(X_train, Y_train, method = "glmnet", metric = "Accuracy", trControl = control) #alpha = 0/1 
plot(cv.fit.all)

#predict test dataset
pred.all <- predict(cv.fit.all, X_test)
confusionMatrix(pred.all, Y_test)
#this model has basically the same accuracy of no-Information-Rate-classification (for both Lasso and Ridge)
#for this reason there is no point in cross-validating this model

# ---------------------------------------------------------------------------------------------
# Lasso AND RIDGE REGRESSION FOR RHINO
# ---------------------------------------------------------------------------------------------

library("ROCR")

#cross-validation
control <- trainControl(method = "cv", number = 10)
rhino.cv <- train(X_rhino_train, Y_rhino_train, method = "glmnet", metric = "Accuracy", trControl = control)
rhino.cv
pred.rhino.cv <- predict(rhino.cv, X_rhino_test)
confusionMatrix(pred.rhino.cv, Y_rhino_test) #100%

#plot ROCR curve
plot(performance(prediction(as.numeric(pred.rhino.cv), as.numeric(Y_rhino_test)), 'tpr', 'fpr'))
#compute Area Under the Curve (AUC)
auc.tmp.rhino <- performance(prediction(as.numeric(pred.rhino.cv), as.numeric(Y_rhino_test)), "auc")
auc.rhino <- as.numeric(auc.tmp.rhino@y.values)
auc.rhino

# ---------------------------------------------------------------------------------------------
# Lasso AND RIDGE REGRESSION FOR RSV
# ---------------------------------------------------------------------------------------------

#cross-validation
control <- trainControl(method = "cv", number = 10)
RSV.cv <- train(X_RSV_train, Y_RSV_train, method = "glmnet", metric = "Accuracy", trControl = control)
RSV.cv
pred.RSV.cv <- predict(RSV.cv, X_RSV_test)
confusionMatrix(pred.RSV.cv, Y_RSV_test) #90%

#plot ROCR curve
plot(performance(prediction(as.numeric(pred.RSV.cv), as.numeric(Y_RSV_test)), 'tpr', 'fpr'))
#compute Area Under the Curve (AUC)
auc.tmp.RSV <- performance(prediction(as.numeric(pred.RSV.cv), as.numeric(Y_RSV_test)), "auc")
auc.RSV <- as.numeric(auc.tmp.RSV@y.values)
auc.RSV

# ---------------------------------------------------------------------------------------------
# Lasso AND RIDGE REGRESSION FOR INFA
# ---------------------------------------------------------------------------------------------

#cross-validation
control <- trainControl(method = "cv", number = 10)
infA.cv <- train(X_infA_train, Y_infA_train, method = "glmnet", metric = "Accuracy", trControl = control)
infA.cv
pred.infA.cv <- predict(infA.cv, X_infA_test)
confusionMatrix(pred.infA.cv, Y_infA_test) #77.78%


#plot ROCR curve
plot(performance(prediction(as.numeric(pred.infA.cv), as.numeric(Y_infA_test)), 'tpr', 'fpr'))
#compute Area Under the Curve (AUC)
auc.tmp.infA <- performance(prediction(as.numeric(pred.infA.cv), as.numeric(Y_infA_test)), "auc")
auc.infA <- as.numeric(auc.tmp.infA@y.values)
auc.infA



