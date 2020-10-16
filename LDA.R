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
# FEATURES SELECTION and TRAIN TEST SPLIT FOR ENTIRE DATASET
# ---------------------------------------------------------------------------------------------

library("genefilter")

#feature selection
tt.geneData <- rowttests(t(geneData.noBatchEffect[,1:22277]), factor(geneData.noBatchEffect$sick.healthy))
keepers <- which(tt.geneData$p.value < 0.1)
geneData_fs <- geneData.noBatchEffect[,keepers]
geneData_fs$sick.healthy <- factor(geneData.noBatchEffect$sick.healthy)

#train/test split
set.seed(123)
smp_size <- floor(0.75* nrow(geneData_fs))
train_ind <- sample(seq_len(nrow(geneData_fs)), size = smp_size)
X_train <- geneData_fs[train_ind,]
X_test <- geneData_fs[-train_ind,]

#save train/test labels in vectors
Y_train <- X_train$sick.healthy
Y_test <- X_test$sick.healthy

#remove train/test lables from dataset
X_train$sick.healthy <- NULL
X_test$sick.healthy <- NULL


# ---------------------------------------------------------------------------------------------
# FEATURES SELECTION and TRAIN TEST SPLIT FOR EACH BATCH
# ---------------------------------------------------------------------------------------------

# RHINO

#feature selection
tt.rhino <- rowttests(t(rhino[,1:22277]), factor(rhino$virus))
keepers <- which(tt.rhino$p.value < 0.1)
rhino_fs <- rhino[,keepers]
rhino_fs$virus <- factor(rhino$virus)

#train/test split
set.seed(123)
smp_size <- floor(0.75* nrow(rhino_fs))
train_ind <- sample(seq_len(nrow(rhino_fs)), size = smp_size)
X_rhino_train <- rhino_fs[train_ind,]
X_rhino_test <- rhino_fs[-train_ind,]

#save train/test labels in vectors
Y_rhino_train <- X_rhino_train$virus
Y_rhino_test <- X_rhino_test$virus

#remove train/test lables from dataset
X_rhino_train$virus <- NULL
X_rhino_test$virus <- NULL


# RSV

#feature selection
tt.RSV <- rowttests(t(RSV[,1:22277]), factor(RSV$virus))
keepers <- which(tt.RSV$p.value < 0.1)
RSV_fs <- RSV[,keepers]
RSV_fs$virus <- factor(RSV$virus)

#train/test split
set.seed(123)
smp_size <- floor(0.75* nrow(RSV_fs))
train_ind <- sample(seq_len(nrow(RSV_fs)), size = smp_size)
X_RSV_train <- RSV_fs[train_ind,]
X_RSV_test <- RSV_fs[-train_ind,]

#save train/test labels in vectors
Y_RSV_train <- X_RSV_train$virus
Y_RSV_test <- X_RSV_test$virus

#remove train/test lables from dataset
X_RSV_train$virus <- NULL
X_RSV_test$virus <- NULL


# INFA

#feature selection
tt.infA <- rowttests(t(infA[,1:22277]), factor(infA$virus))
keepers <- which(tt.infA$p.value < 0.1)
infA_fs <- infA[,keepers]
infA_fs$virus <- factor(infA$virus)

#train/test split
set.seed(123)
smp_size <- floor(0.75* nrow(infA_fs))
train_ind <- sample(seq_len(nrow(infA_fs)), size = smp_size)
X_infA_train <- infA_fs[train_ind,]
X_infA_test <- infA_fs[-train_ind,]

#save train/test labels in vectors
Y_infA_train <- X_infA_train$virus
Y_infA_test <- X_infA_test$virus

#remove train/test lables from dataset
X_infA_train$virus <- NULL
X_infA_test$virus <- NULL

# ---------------------------------------------------------------------------------------------
# LDA FOR ENTIRE DATASET
# ---------------------------------------------------------------------------------------------

library(MASS)
library(caret)
geneData_lda <- lda(x = X_train, grouping = Y_train)
geneData_lda.pred.test <- predict(geneData_lda, X_test)
confusionMatrix(geneData_lda.pred.test$class, Y_test) #79.31%
#The model has decent accuracy in classifying sick vs healthy individuals, but has poor performance if we try to predict the type of virus
#good in classifying healthy vs sick, bad at classifying the type of virus

#cross-validation
control <- trainControl(method = "cv", number = 10)
geneData.cv <- train(X_train, Y_train, method = "lda", metric = "Accuracy", trControl = control)
geneData.cv_pred <- predict(geneData.cv, X_test)
confusionMatrix(geneData.cv_pred, Y_test) #79.31%

# ---------------------------------------------------------------------------------------------
# LDA FOR RHINO
# ---------------------------------------------------------------------------------------------

rhino_lda <- lda(x = X_rhino_train, grouping = Y_rhino_train)
plot(rhino_lda)

#predict on train set and plot
mod.rhino.values <- predict(rhino_lda, X_rhino_train)
plot(mod.rhino.values$x[,1], ylab = c("LDA Axis"))
text(mod.rhino.values$x[,1], col = as.numeric(Y_rhino_train))

#predict on test set
rhino_lda.pred.test <- predict(rhino_lda, X_rhino_test)
confusionMatrix(rhino_lda.pred.test$class, Y_rhino_test) #100%

#cross-validation
control <- trainControl(method = "cv", number = 10)
rhino.cv <- train(X_rhino_train, Y_rhino_train, method = "lda", metric = "Accuracy", trControl = control)
rhino.cv_pred <- predict(rhino.cv, X_rhino_test)
confusionMatrix(rhino.cv_pred, Y_rhino_test) #100%

# ---------------------------------------------------------------------------------------------
# LDA FOR RSV
# ---------------------------------------------------------------------------------------------

RSV_lda <- lda(x = X_RSV_train, grouping = Y_RSV_train)
plot(RSV_lda)

#predict on train set and plot
mod.RSV.values <- predict(RSV_lda, X_RSV_train)
plot(mod.RSV.values$x[,1], ylab = c("LDA Axis"))
text(mod.RSV.values$x[,1], col = as.numeric(Y_RSV_train))

#predict on test set
RSV_lda.pred.test <- predict(RSV_lda, X_RSV_test)
confusionMatrix(RSV_lda.pred.test$class, Y_RSV_test) #100%

#cross-validation
control <- trainControl(method = "cv", number = 10)
RSV.cv <- train(X_RSV_train, Y_RSV_train, method = "lda", metric = "Accuracy", trControl = control)
RSV.cv_pred <- predict(RSV.cv, X_RSV_test)
confusionMatrix(RSV.cv_pred, Y_RSV_test) #100%

# ---------------------------------------------------------------------------------------------
# LDA FOR INFA
# ---------------------------------------------------------------------------------------------

infA_lda <- lda(x = X_infA_train, grouping = Y_infA_train)
plot(infA_lda)

#predict on train set and plot
mod.infA.values <- predict(infA_lda, X_infA_train)
plot(mod.infA.values$x[,1], ylab = c("LDA Axis"))
text(mod.infA.values$x[,1], col = as.numeric(Y_infA_train))

#predict on test set
infA_lda.pred.test <- predict(infA_lda, X_infA_test)
confusionMatrix(infA_lda.pred.test$class, Y_infA_test) #77.78%

#cross-validation
control <- trainControl(method = "cv", number = 10)
infA.cv <- train(X_infA_train, Y_infA_train, method = "lda", metric = "Accuracy", trControl = control)
infA.cv_pred <- predict(infA.cv, X_infA_test)
confusionMatrix(infA.cv_pred, Y_infA_test) #77.78%

