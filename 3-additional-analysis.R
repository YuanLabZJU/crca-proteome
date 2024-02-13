library("MLmetrics")

F1_Score(pro_model$Y,model_pre)

# F1 score: 
# Clinical model:  train - 0.7964; valid - 0.8372
# Protomics model: train - 0.8459; valid - 0.8000
# Clinical+proteomics model: train - 0.9038; valid - 0.8837

# Dummy model
dummy.fit <- glm(Y~1, data = pro_model, family = 'binomial' )
model_pre <- predict.glm(dummy.fit, newdata=model, type = 'response')
model_pre <- as.numeric(model_pre)
model_roc <- pROC::roc(pro_model$Y, model_pre, direction = "<", ci = T)
plot(model_roc, print.auc = TRUE, print.thres = T, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve for Training Set")
valid_pre <- predict.glm(dummy.fit, newdata=valid, type = 'response')
valid_pre <- as.numeric(valid_pre)
valid_roc <- pROC::roc(pro_valid$Y, valid_pre, direction = "<", ci = T)
plot(valid_roc, print.auc = TRUE, print.thres = T, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve for Validating Set")

cp.data <- cbind(Y=pro_model$Y, Pred=model_pre) |> data.frame()
cp <- cutpointr(cp.data, Pred, Y, direction = "<=", 
                method = maximize_metric, metric = sum_sens_spec)

# Proteomics model
step.fit1 <- glm(Y~PDP1_HUMAN+ALR_HUMAN+ENOG_HUMAN+FYCO1_HUMAN+STXB1_HUMAN
                 +ARH40_HUMAN+LICH_HUMAN+RIMC1_HUMAN+NPC2_HUMAN+MTMR5_HUMAN
                 , data = pro_model, family = 'binomial' )
step.fit1 <- step(step.fit1, direction="both")

model <- pro_model[,-c(1,11)]
model_pre <- predict.glm(step.fit1, newdata=model, type = 'response')
model_pre <- as.numeric(model_pre)
cp.data <- cbind(Y=pro_model$Y, Pred=model_pre) |> data.frame()
cp <- cutpointr(cp.data, Pred, Y, 
                method = maximize_metric, metric = sum_sens_spec)
cp_sum <- summary(cp)
F1_score(cp_sum$confusion_matrix[[1]]$tp, 
         cp_sum$confusion_matrix[[1]]$fp,
         cp_sum$confusion_matrix[[1]]$tn,
         cp_sum$confusion_matrix[[1]]$fn)
valid <- pro_valid[,-c(1,11)]
valid_pre <- predict.glm(step.fit1, newdata=valid, type = 'response')
valid_pre <- as.numeric(valid_pre)
cp.data <- cbind(Y=pro_valid$Y, Pred=valid_pre) |> data.frame()
cp <- cutpointr(cp.data, Pred, Y, 
                method = maximize_metric, metric = sum_sens_spec)
cp_sum <- summary(cp)
F1_score(cp_sum$confusion_matrix[[1]]$tp, 
         cp_sum$confusion_matrix[[1]]$fp,
         cp_sum$confusion_matrix[[1]]$tn,
         cp_sum$confusion_matrix[[1]]$fn)

# Clinical model
step.fit2 <- glm(Y~ AS + MS + LL + PT,
                 data = pro_model, family = "binomial")

model <- pro_model[,-c(1,11)]
model_pre <- predict.glm(step.fit2, newdata=model, type = 'response')
model_pre <- as.numeric(model_pre)
cp.data <- cbind(Y=pro_model$Y, Pred=model_pre) |> data.frame()
cp <- cutpointr(cp.data, Pred, Y, 
                method = maximize_metric, metric = sum_sens_spec)
cp_sum <- summary(cp)
F1_score(cp_sum$confusion_matrix[[1]]$tp, 
         cp_sum$confusion_matrix[[1]]$fp,
         cp_sum$confusion_matrix[[1]]$tn,
         cp_sum$confusion_matrix[[1]]$fn)
valid <- pro_valid[,-c(1,11)]
valid_pre <- predict.glm(step.fit2, newdata=valid, type = 'response')
valid_pre <- as.numeric(valid_pre)
cp.data <- cbind(Y=pro_valid$Y, Pred=valid_pre) |> data.frame()
cp <- cutpointr(cp.data, Pred, Y, 
                method = maximize_metric, metric = sum_sens_spec)
cp_sum <- summary(cp)
F1_score(cp_sum$confusion_matrix[[1]]$tp, 
         cp_sum$confusion_matrix[[1]]$fp,
         cp_sum$confusion_matrix[[1]]$tn,
         cp_sum$confusion_matrix[[1]]$fn)


# Clinical + Proteomics model
step.fit3 <- glm(Y~ AS + MS + LL + PT
                 +PDP1_HUMAN+ALR_HUMAN+ENOG_HUMAN+FYCO1_HUMAN+STXB1_HUMAN
                 +ARH40_HUMAN+CE051_HUMAN+NPC2_HUMAN+MTMR5_HUMAN,
                 data = pro_model, family = "binomial")

model <- pro_model[,-c(1,11)]
model_pre <- predict.glm(step.fit3, newdata=model, type = 'response')
model_pre <- as.numeric(model_pre)
cp.data <- cbind(Y=pro_model$Y, Pred=model_pre) |> data.frame()
cp <- cutpointr(cp.data, Pred, Y, 
                method = maximize_metric, metric = sum_sens_spec)
cp_sum <- summary(cp)
F1_score(cp_sum$confusion_matrix[[1]]$tp, 
         cp_sum$confusion_matrix[[1]]$fp,
         cp_sum$confusion_matrix[[1]]$tn,
         cp_sum$confusion_matrix[[1]]$fn)
valid <- pro_valid[,-c(1,11)]
valid_pre <- predict.glm(step.fit3, newdata=valid, type = 'response')
valid_pre <- as.numeric(valid_pre)
cp.data <- cbind(Y=pro_valid$Y, Pred=valid_pre) |> data.frame()
cp <- cutpointr(cp.data, Pred, Y, 
                method = maximize_metric, metric = sum_sens_spec)
cp_sum <- summary(cp)
F1_score(cp_sum$confusion_matrix[[1]]$tp, 
         cp_sum$confusion_matrix[[1]]$fp,
         cp_sum$confusion_matrix[[1]]$tn,
         cp_sum$confusion_matrix[[1]]$fn)

## Boxplot with comparison test
library(gghubr)
library(rstatix)
pro_model$YY[pro_model$Y==1] <- 'Yes'
pro_model$YY[pro_model$Y==0] <- 'No'
pro_valid$YY[pro_valid$Y==1] <- 'Yes'
pro_valid$YY[pro_valid$Y==0] <- 'No'
pro_plot <- pro_merged |> 
  dplyr::select(PatID, PDP1_HUMAN, ALR_HUMAN, ENOG_HUMAN,
                FYCO1_HUMAN, ARH40_HUMAN, STXB1_HUMAN, 
                CE051_HUMAN, NPC2_HUMAN, MTMR5_HUMAN, Y) |> 
  mutate(Cohort = ifelse(str_detect(PatID, "P"), "Training", "Validation"),
         Survivor = ifelse(Y==0, "No", "Yes")) |> 
  data.frame()
prots <- c("PDP1_HUMAN", "ALR_HUMAN", "ENOG_HUMAN",
           "FYCO1_HUMAN", "ARH40_HUMAN", "STXB1_HUMAN", 
           "RIMC1_HUMAN", "NPC2_HUMAN", "MTMR5_HUMAN")

pdf("boxplots.pdf", width=5, height=5)
for (prot_name in prots){
  pro_plot$prot_curr <- pro_plot[,prot_name]
  stat.test <- pro_plot %>% 
    group_by(Cohort) %>%
    t_test(prot_curr ~ Survivor) %>% 
    add_significance("p")
  bxp <- ggboxplot(
    pro_plot,x="Cohort",y="prot_curr",
    color= "Survivor", xlab = prot_name, ylab = "Expression (z-score)",
    palette=c("#87CEFA","violet"))
  stat.test <- stat.test %>% 
    add_xy_position(x='Cohort',dodge = 1)
  print(bxp + stat_pvalue_manual(stat.test,  label = "P={p}{p.signif}", 
                           tip.length = 0, hide.ns = FALSE))
}
dev.off()
