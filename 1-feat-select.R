rm(list=ls())
library(data.table)
library(tidyverse)
library(glmnet)
library(caTools)
library(MASS)
library(pROC)
library(boot)
library(qpcR)
library(tableone)
library(dplyr)

pro_model <- fread(file = "/mnt/d/hyh/CRCA/revision/pro_model.csv") # 230*5223
pro_valid <- fread(file = "/mnt/d/hyh/CRCA/revision/pro_valid.csv") # 58*5223

# Random seed
random.seed_int <- sample(c(0:999999), size = 1000, replace = F)

select <- function(n){
  
  n <- random.seed_int[n]
  set.seed(n)
  pro_model <- data.frame(pro_model)
  split <- sample.split(pro_model$Y,SplitRatio = 0.7)
  pro_train <- subset(pro_model,split == TRUE)
  pro_test <- subset(pro_model,split == FALSE)
  
  sigpro <- NA
  Y <- pro_train$Y
  for (i in 12:5223){
    x <- pro_train[,i]
    if (sd(x)!=0){
      po_train_sub <-cbind(x,Y)
      po_train_sub <-data.frame(po_train_sub)
      fit_select <- glm(Y~. , data=po_train_sub)
      a <-coef(summary(fit_select))
      if (a[2,4]<0.05) {
        name <-colnames(pro_train)[i]
        sigpro<-c(sigpro, name)
      }
    }
  }
  sigpro <- sigpro[-1]
  sigpro <- c(sigpro,'PatID' , 'Y')
  
  po_train <- as.data.frame(pro_train)
  po_test <- as.data.frame(pro_test)
  po_valid <- as.data.frame(pro_valid)
  po_train <- po_train[,sigpro]
  Status <- po_train$Y
  nfeat <- (ncol(po_train)-2)
  Features <- as.matrix(po_train[,1:(ncol(po_train)-2)])
  fit1 <- glmnet(Features, Status, nlambda = 100, nfold = 10,
                 type.measure='deviance', standardize = T, alpha = 1)
  cv.fit1 <- cv.glmnet(Features, Status, nfold = 10,
                       family="binomial", type.measure='deviance',
                       standardize = T, alpha = 1)   
  opt.lambda1 <- cv.fit1$lambda.1se
  opt.lambda <- cv.fit1$lambda.min
  # lambda.min
  coefs_fit1 <- as.matrix(coef.glmnet(cv.fit1, s=opt.lambda))
  coefs_fit1 <- data.frame(coefs_fit1)
  coefs_fit1 <- filter(coefs_fit1, s1!=0)
  pro_name1 <- rownames(coefs_fit1)[2:nrow(coefs_fit1)]
  # lambda.1se
  coefs_fit2 <- as.matrix(coef.glmnet(cv.fit1, s=opt.lambda1))
  coefs_fit2 <- data.frame(coefs_fit2)
  coefs_fit2 <- filter(coefs_fit2, s1!=0)
  pro_name2 <- rownames(coefs_fit2)[2:nrow(coefs_fit2)]
  
  pro_name <- qpcR:::cbind.na(pro_name1, pro_name2)
  return(pro_name)
}


# For 1000 times
pro_select_min <- NA
pro_select_1se <- NA
for (i in 1:1000){
  pro_min <- select(i)[,1]
  pro_1se <- select(i)[,2]
  pro_select_min <- qpcR:::cbind.na(pro_select_min, pro_min)
  pro_select_1se <- qpcR:::cbind.na(pro_select_1se, pro_1se)
}

# Calculate frequency
pro_min_fre <- NA
for (i in 1:1000){
  for (j in 1:dim(pro_select_min)[1]){
    a <- pro_select_min[j,i]
    if (!is.na(a)){
      pro_min_fre <- c(pro_min_fre, a)
    }
  }
}
pro_min_fre <- pro_min_fre[-1]
pro_min_fre <- data.frame(pro_min_fre)
pro_min_tab <- data.frame(table(pro_min_fre$pro_min_fre))
write.csv(pro_min_tab, 'Out/pro_min_1000(int).csv')

pro_1se_fre <- NA
for (i in 1:1000){
  for (j in 1:dim(pro_select_1se)[1]){
    a <- pro_select_1se[j,i]
    if (!is.na(a)){
      pro_1se_fre <- c(pro_1se_fre, a)
    }
  }
}
pro_1se_fre <- pro_1se_fre[-1]
pro_1se_fre <- data.frame(pro_1se_fre)
pro_1se_tab <- data.frame(table(pro_1se_fre$pro_1se_fre))
write.csv(pro_1se_tab, 'Out/pro_1se_1000(int).csv')