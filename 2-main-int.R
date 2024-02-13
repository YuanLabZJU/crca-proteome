library(data.table)
library(tidyverse)
library(glmnet)
library(caTools)
library(pROC)
library(boot)
library(randomForest)
library(qpcR)
library(ggplot2)
library(Hmisc)
library(corrplot)
library(tableone)
library(gridExtra)
library(readxl)
library(ggrepel)

pro_model <- fread(file = "/mnt/d/hyh/CRCA/revision/pro_model.csv") # 230*5223
pro_valid <- fread(file = "/mnt/d/hyh/CRCA/revision/pro_valid.csv") # 58*5223
pro_merged <- rbind(pro_model, pro_valid)
patInf_all <- fread('/mnt/d/hyh/CRCA/revision/CRCA_288patientinf_all.csv')
survival <- patInf_all |> 
  dplyr::select(PatID, `Survival_time(Day)`) |> 
  rename("time" = "Survival_time(Day)")

# Dummy model
dummy.fit <- glm(Y~1
                 , data = pro_model, family = 'binomial' )
model_pre <- predict.glm(dummy.fit, newdata=pro_model, type = 'response')
model_pre <- as.numeric(model_pre)
model_roc <- pROC::roc(pro_model$Y, model_pre, direction = "<", ci = T)
plot(model_roc, print.auc = TRUE, print.thres = T, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve for Training Set")
valid_pre <- predict.glm(dummy.fit, newdata=pro_valid, type = 'response')
valid_pre <- as.numeric(valid_pre)
valid_roc <- pROC::roc(pro_valid$Y, valid_pre, direction = "<", ci = T)
plot(valid_roc, print.auc = TRUE, print.thres = T, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve for Validating Set")

# Proteomics model
step.fit1 <- glm(Y~PDP1_HUMAN+ALR_HUMAN+ENOG_HUMAN+FYCO1_HUMAN+STXB1_HUMAN
                 +ARH40_HUMAN+LICH_HUMAN+RIMC1_HUMAN+NPC2_HUMAN+MTMR5_HUMAN
                 , data = pro_model, family = 'binomial' )
step.fit1 <- step(step.fit1, direction="both")

model <- pro_model[,-c(1,11)]
model_pre <- predict.glm(step.fit1, newdata=model, type = 'response')
model_pre <- as.numeric(model_pre)
model_roc <- pROC::roc(pro_model$Y, model_pre, direction = "<", ci = T)
plot(model_roc, print.auc = TRUE, print.thres = T, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve for Training Set")   

valid <- pro_valid[,-c(1,11)]
valid_pre <- predict.glm(step.fit1, newdata=valid, type = 'response')
valid_pre <- as.numeric(valid_pre)
valid_roc <- pROC::roc(pro_valid$Y, valid_pre, direction = "<", ci = T)
plot(valid_roc, print.auc = TRUE, print.thres = T, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve for Validating Set") 

plot(model_roc, print.auc = TRUE, print.thres = F, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve of proteomics model", col = "red", 
     print.auc.x = 0.5, print.auc.y = 0.5)
plot(valid_roc, add = T, col = 'blue', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
legend("bottomright", legend = c("Training Set", "Validation Set"), 
       col=c("red", "blue"), lwd=2)



# Clinical model
step.fit2 <- glm(Y~ AS + MS + LL + PT,
                 data = pro_model, family = "binomial")

model <- pro_model[,-c(1,11)]
model_pre <- predict.glm(step.fit2, newdata=model, type = 'response')
model_pre <- as.numeric(model_pre)
model_roc <- pROC::roc(pro_model$Y, model_pre, direction = "<", ci = T)
plot(model_roc, print.auc = TRUE, print.thres = T, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve for Training Set")   

valid <- pro_valid[,-c(1,11)]
valid_pre <- predict.glm(step.fit2, newdata=valid, type = 'response')
valid_pre <- as.numeric(valid_pre)
valid_roc <- pROC::roc(pro_valid$Y, valid_pre, direction = "<", ci = T)
plot(valid_roc, print.auc = TRUE, print.thres = T, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve for Validating Set") 

plot(model_roc, print.auc = TRUE, print.thres = F, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve of clinical model", col = "red", 
     print.auc.x = 0.5, print.auc.y = 0.5)
plot(valid_roc, add = T, col = 'blue', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
legend("bottomright", legend = c("Training Set", "Validation Set"), 
       col=c("red", "blue"), lwd=2)



# Clinical + Proteomics model
step.fit3 <- glm(Y~ AS + MS + LL + PT
                 +PDP1_HUMAN+ALR_HUMAN+ENOG_HUMAN+FYCO1_HUMAN+STXB1_HUMAN
                 +ARH40_HUMAN+RIMC1_HUMAN+NPC2_HUMAN+MTMR5_HUMAN,
                 data = pro_model, family = "binomial")

model <- pro_model[,-c(1,11)]
model_pre <- predict.glm(step.fit3, newdata=model, type = 'response')
model_pre <- as.numeric(model_pre)
model_roc <- pROC::roc(pro_model$Y, model_pre, direction = "<", ci = T)
plot(model_roc, print.auc = TRUE, print.thres = T, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve for Training Set")   

valid <- pro_valid[,-c(1,11)]
valid_pre <- predict.glm(step.fit3, newdata=valid, type = 'response')
valid_pre <- as.numeric(valid_pre)
valid_roc <- pROC::roc(pro_valid$Y, valid_pre, direction = "<", ci = T)
plot(valid_roc, print.auc = TRUE, print.thres = T, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve for Validating Set") 

plot(model_roc, print.auc = TRUE, print.thres = F, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve of clinical + proteomics model", col = "red", 
     print.auc.x = 0.5, print.auc.y = 0.5)
plot(valid_roc, add = T, col = 'blue', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
legend("bottomright", legend = c("Training Set", "Validation Set"), 
       col=c("red", "blue"), lwd=2)






# Subset(70%) with 1000 times to evaluate the stability of established model
random.seed_sub <- sample(0:999999, size = 1000, replace = F)

sub_auc <- NA
for (i in 1:1000){
  n <- random.seed_sub[i]
  set.seed(n)
  split <- sample.split(pro_model$Y, SplitRatio = 0.7)
  pro_sub <- pro_model[split == T]
  
  sub <- pro_sub[,-c(1,11)]
  sub_pre <- predict.glm(step.fit1, newdata=sub, type = 'response')
  sub_pre <- as.numeric(sub_pre)
  sub_roc <- pROC::roc(pro_sub$Y, sub_pre, direction = "<")
  
  sub_auc <- c(sub_auc, sub_roc$auc)
}
sub_auc <- sub_auc[-1]
## Violin plot
Subset <- rep("Subset (70%)", 1000)
sub_auc_frame <- data.frame(sub_auc, Subset)
sub_violin <- ggplot(sub_auc_frame, aes(x = Subset, y = sub_auc)) +
  geom_violin(trim = F, fill = "#87CEFA")+
  geom_boxplot(width=0.1)+
  labs(y = "AUC")+
  theme_classic()
## Histogram
sub_hist <- ggplot(sub_auc_frame, aes(sub_auc))+
  geom_histogram(fill = "#87CEFA", color = "black", bins = 30)+
  geom_density(size = 1, color = "#FF00FF")+
  labs(x = "AUC", y = "Frequency", 
       title = "Frequency histogram for 1000 AUC from 1000 subsets")+
  theme_classic()

# Resample with 1000 times to evaluate the stability of model
random.seed_re <- rep(0, 230)
for (i in 1:1000){
  re <- sample(1:230, size = 230, replace = T)
  random.seed_re <- cbind(random.seed_re, re)
}
random.seed_re <- random.seed_re[,-1]

re_auc <- NA
for (i in 1:1000){
  pro_re <- pro_model[random.seed_re[,i]]
  re <- pro_re[,-c(1,11)]
  re_pre <- predict.glm(step.fit3, newdata=re, type = 'response')
  re_pre <- as.numeric(re_pre)
  re_roc <- pROC::roc(pro_re$Y, re_pre, direction = "<")
  
  re_auc <- c(re_auc, re_roc$auc)
}
re_auc <- re_auc[-1]
## Violin plot
Resample <- rep("Resample", 1000)
re_auc_frame <- data.frame(re_auc, Resample)
re_violin <- ggplot(re_auc_frame, aes(x = Resample, y = re_auc)) +
  geom_violin(trim = F, fill = "#87CEFA")+
  geom_boxplot(width=0.1)+
  labs(y = "AUC")+
  theme_classic()
## Histogram
re_hist_train <- ggplot(re_auc_frame, aes(re_auc))+
  geom_histogram(fill = "#87CEFA", color = "black", bins = 30)+
  geom_density(size = 1, color = "#FF00FF")+
  labs(x = "AUC", y = "Frequency",
       title = "Training set (resampled 1000 times)")+
  theme_classic()  
# Validation
re_auc <- NA
for (i in 1:1000){
  pro_re <- pro_valid[random.seed_re[,i]]
  re <- pro_re[,-c(1,11)]
  re_pre <- predict.glm(step.fit3, newdata=re, type = 'response')
  re_pre <- as.numeric(re_pre)
  re_roc <- pROC::roc(pro_re$Y, re_pre, direction = "<")
  
  re_auc <- c(re_auc, re_roc$auc)
}
re_auc <- re_auc[-1]  
## Violin plot
Resample <- rep("Resample", 1000)
re_auc_frame <- data.frame(re_auc, Resample)
re_violin <- ggplot(re_auc_frame, aes(x = Resample, y = re_auc)) +
  geom_violin(trim = F, fill = "#87CEFA")+
  geom_boxplot(width=0.1)+
  labs(y = "AUC")+
  theme_classic()
## Histogram
re_hist_valid <- ggplot(re_auc_frame, aes(re_auc))+
  geom_histogram(fill = "#87CEFA", color = "black", bins = 30)+
  geom_density(size = 1, color = "#FF00FF")+
  labs(x = "AUC", y = "Frequency",
       title = "Validation set (resampled 1000 times)")+
  theme_classic()  




# Correlation plot
pro_cor <- pro_model[,c("PDP1_HUMAN", "ALR_HUMAN", "ENOG_HUMAN", 
                        "FYCO1_HUMAN", "STXB1_HUMAN", "ARH40_HUMAN", 
                        "LICH_HUMAN", "RIMC1_HUMAN", 
                        "NPC2_HUMAN", "MTMR5_HUMAN")]
cor_ma <- rcorr(as.matrix(pro_cor), type = "pearson")
corrplot(cor_ma$r, type = "upper")

# Correlation with selected proteins (univariate significant)
pro_model <- data.frame(pro_model)
sigpro <- NA
Y <- pro_model$Y
for (i in 12:5223){
  x <- pro_model[,i]
  if (sd(x)!=0){
    po_model_sub <- cbind(x,Y)
    po_model_sub <- data.frame(po_model_sub)
    fit_select <- glm(Y~. , data=po_model_sub)
    a <- coef(summary(fit_select))
    if (a[2,4]<0.05) {
      name <- colnames(pro_model)[i]
      sigpro<- c(sigpro, name)
    }
  }
}
sigpro <- sigpro[-1]

pro_corr <- pro_model[,sigpro]
select_pro <- c("PDP1_HUMAN", "ALR_HUMAN", "ENOG_HUMAN", 
                "FYCO1_HUMAN", "STXB1_HUMAN", "ARH40_HUMAN", 
                "RIMC1_HUMAN", "NPC2_HUMAN", "MTMR5_HUMAN")
corr_table <- rep(NA, 464)
for (i in 1:9){
  select_pro1 <- select_pro[i]
  name <- colnames(pro_corr)
  corr_name <- NA
  corr_value <- NA
  for (j in 1:464){
    corr_name <- c(corr_name, name[j])
    corr_value <- c(corr_value, cor(pro_corr[,j], pro_corr[,select_pro1]))
  }
  corr_name <- corr_name[-1]
  corr_value <- corr_value[-1]
  corr_table <- cbind(corr_table, corr_name, corr_value)
}
corr_table <- corr_table[,-1]
corr_table <- data.frame(corr_table)

## Sort
corr_sort <- corr_table[,1:2]
corr_sort <- corr_sort[order(abs(as.numeric(corr_sort$corr_value)),
                             decreasing = T),]
corr_table[,1:2] <- corr_sort
corr_sort <- corr_table[,3:4]
corr_sort <- corr_sort[order(abs(as.numeric(corr_sort$corr_value)), 
                             decreasing = T),]
corr_table[,3:4] <- corr_sort
corr_sort <- corr_table[,5:6]
corr_sort <- corr_sort[order(abs(as.numeric(corr_sort$corr_value)),
                             decreasing = T),]
corr_table[,5:6] <- corr_sort
corr_sort <- corr_table[,7:8]
corr_sort <- corr_sort[order(abs(as.numeric(corr_sort$corr_value)), 
                             decreasing = T),]
corr_table[,7:8] <- corr_sort
corr_sort <- corr_table[,9:10]
corr_sort <- corr_sort[order(abs(as.numeric(corr_sort$corr_value)), 
                             decreasing = T),]
corr_table[,9:10] <- corr_sort
corr_sort <- corr_table[,11:12]
corr_sort <- corr_sort[order(abs(as.numeric(corr_sort$corr_value)),
                             decreasing = T),]
corr_table[,11:12] <- corr_sort
corr_sort <- corr_table[,13:14]
corr_sort <- corr_sort[order(abs(as.numeric(corr_sort$corr_value)), 
                             decreasing = T),]
corr_table[,13:14] <- corr_sort
corr_sort <- corr_table[,15:16]
corr_sort <- corr_sort[order(abs(as.numeric(corr_sort$corr_value)),
                             decreasing = T),]
corr_table[,15:16] <- corr_sort
corr_sort <- corr_table[,17:18]
corr_sort <- corr_sort[order(abs(as.numeric(corr_sort$corr_value)),
                             decreasing = T),]
corr_table[,17:18] <- corr_sort
colnames(corr_table) <- c("Protein name", "Correlation value",
                          "Protein name", "Correlation value",
                          "Protein name", "Correlation value",
                          "Protein name", "Correlation value",
                          "Protein name", "Correlation value",
                          "Protein name", "Correlation value",
                          "Protein name", "Correlation value",
                          "Protein name", "Correlation value",
                          "Protein name", "Correlation value")



# Table One
pro_merged$group <- 1
pro_merged$group[231:288] <- 2
tableone <- CreateTableOne(vars = c('A', 'S', 'LL', 'PT', 'TS', 'LS', 'AS', 'MS', 'Y'),
                           factorVars = c('S', 'LL', 'PT', 'TS', 'LS', 'AS', 'MS', 'Y'),
                           strata = 'group', data = pro_merged, addOverall = T)
print(tableone)
chi <- as.matrix(pro_merged[,c("TS", "group")])
chi[,1] <- as.factor(chi[,1])
chi <- chisq.test(chi[,1], chi[,2], correct = F)



# Kaplan-Meier plot
# Positive likelihood ratio > 5 → Low risk group
# Negative likelihood ratio < 0.2 → High risk group
# Others → Moderate risk group
library(cutpointr)
library(survival)
library(survminer)

model <- pro_model[,-c(1,11)]
model_pre <- predict.glm(step.fit3, newdata=model, type = 'response')
model_pre <- as.numeric(model_pre)
valid <- pro_valid[,-c(1,11)]
valid_pre <- predict.glm(step.fit3, newdata=valid, type = 'response')
valid_pre <- as.numeric(valid_pre)

model_cutoff <- cbind(pro_model[,c(1,11)], model_pre)
model_cutoff$high <- ifelse(model_cutoff$model_pre >= 0.752, '1', '0')
ma1 <- table(model_cutoff$high, model_cutoff$Y)
plr(tp = ma1[2,2], fp = ma1[2,1], tn = ma1[1,1], fn = ma1[1,2]) # Positive LR
# cut off for low risk group: 0.752 for combined model

## Training set
model_km <- left_join(pro_model[,c(1,11)], survival, by = 'PatID')
model_km <- cbind(model_km, model_pre)
model_km$time[model_km$Y==1] <- 1825
model_km$group[model_km$model_pre >= 0.752] <- 1 # Low risk - high prob of survival
model_km$group[model_km$model_pre < 0.752] <- 2 # High risk - low prob of survival
model_km$status[model_km$Y==1] <- 0
model_km$status[model_km$Y==0] <- 1
surv.fit1 <- survfit(Surv(time, status) ~ group, data = model_km)
survplot1 <- ggsurvplot(surv.fit1, data = model_km, pval = T, pval.coord =c(0, 0.05),
           conf.int = T, palette = c("blue", "red"),conf.int.style = "step",
           legend.labs = c("Low risk group", "High risk group"),
           title = "Kaplan-Meier survival curve for training set", 
           xlab = "Time (Days)",
           ylab = "Cumulative survival (percentage)",
           break.x.by = 300, xlim = c(0,1900),
           risk.table = "abs_pct",tables.height = 0.2,
           tables.theme = theme_cleantable(), ggtheme = theme_bw())

## Validation set
valid_km <- left_join(pro_valid[,c(1,11)], survival, by = 'PatID')
valid_km <- cbind(valid_km, valid_pre)
valid_km$time[valid_km$Y==1] <- 1825
valid_km$group[valid_km$valid_pre >= 0.752] <- 1 # Low risk - high prob of survival
valid_km$group[valid_km$valid_pre < 0.752] <- 2 # High risk - low prob of survival
valid_km$status[valid_km$Y==1] <- 0
valid_km$status[valid_km$Y==0] <- 1
surv.fit2 <- survfit(Surv(time, status) ~ group, data = valid_km)
survplot2 <- ggsurvplot(surv.fit2, data = valid_km, pval = T, pval.coord =c(0, 0.05),
                        conf.int = T, palette = c("blue", "red"),conf.int.style = "step",
           legend.labs = c("Low risk group", "High risk group"),
           title = "Kaplan-Meier survival curve for validation set", 
           xlab = "Time (Days)",
           ylab = "Cumulative survival (percentage)",
           break.x.by = 200, xlim = c(0,1900),
           risk.table = "abs_pct",tables.height = 0.2,
           tables.theme = theme_cleantable(), ggtheme = theme_bw())
pdf(file = "/mnt/d/hyh/CRCA/revision/KM-risk-stratified.pdf", width=10, heigh=8)
survplot1
survplot2
dev.off()

# Basline characteristics by risk stratification
## train
patInf_all_risk <- cbind(patInf_all, 
                         Risk_Strata = c(model_km$group, valid_km$group),
                         Dataset = c(rep("Train", 230), rep("Val", 58))
                         )
table(model_km$group)
table(valid_km$group)
vars <- names(patInf_all_risk)[2:13]
cat_vars <- c('S', 'LL', 'PT', 'TS', 'LS',
              'M', 'AS', "MS", "Y", "Adjuvant chemotherapy")
tableone_rs_train <- CreateTableOne(vars = vars, factorVars = cat_vars, 
                                 strata = c("Risk_Strata"), 
                                 data = patInf_all_risk |> filter(Dataset=="Train"), 
                                 addOverall = T)
tableone_rs_train <- print(tableone_rs_train, contDigits = 1, showAllLevels = T)
tableone_rs_train <- data.frame(rownames(tableone_rs_train), tableone_rs_train)
write.csv(tableone_rs_train, "/mnt/d/hyh/CRCA/revision/RiskStratifiedTableOne-Train.csv")

tableone_rs_valid <- CreateTableOne(vars = vars, factorVars = cat_vars, 
                                    strata = c("Risk_Strata"), 
                                    data = patInf_all_risk |> filter(Dataset=="Val"), 
                                    addOverall = T)
tableone_rs_valid <- print(tableone_rs_valid, contDigits = 1, showAllLevels  = T)
tableone_rs_valid <- data.frame(rownames(tableone_rs_valid), tableone_rs_valid)
write.csv(tableone_rs_valid, "/mnt/d/hyh/CRCA/revision/RiskStratifiedTableOne-Val.csv")

