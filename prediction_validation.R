library(ranger)

#Imputing evaluation datasets
for (a in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[a])
  species_data_temp["eval_imputed"] <- list(rfImpute(formula1,
                                                     species_data_temp$eval[,c(2:25,35)]))
  assign(myspecies[a], species_data_temp)
  print(myspecies[a])
}

#Prediciting probabilities
for (b in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[b])
  ##Prediction with glm_dib
  species_data_temp$eval$glm_dib_pred_prob <- predict(species_data_temp$glm_dib,
                                                      newdata = species_data_temp$eval,
                                                      type = "response")
  ##Prediction with glm_clim
  species_data_temp$eval$glm_clim_pred_prob <- predict(species_data_temp$glm_clim,
                                                       newdata = species_data_temp$eval,
                                                       type = "response")
  ##Prediction with glm_dib_clim
  species_data_temp$eval$glm_dib_clim_pred_prob <- predict(species_data_temp$glm_dib_clim,
                                                       newdata = species_data_temp$eval,
                                                       type = "response")
  ##Prediction with ct_dib
  pred <- predict(species_data_temp$pruned_tree_dib,
                  newdata = species_data_temp$eval,
                  type = "prob")
  species_data_temp$eval$ct_dib_pred_prob <- pred[,2]
  ##Prediction with ct_clim
  pred <- predict(species_data_temp$pruned_tree_clim,
                  newdata = species_data_temp$eval,
                  type = "prob")
  species_data_temp$eval$ct_clim_pred_prob <- pred[,2]
  ##Prediction with ct_dib_clim
  pred <- predict(species_data_temp$pruned_tree_dib_clim,
                  newdata = species_data_temp$eval,
                  type = "prob")
  species_data_temp$eval$ct_dib_clim_pred_prob <- pred[,2]
  ##Prediction with rf_dib
  pred <- predict(species_data_temp$rforest_dib,
                  data = species_data_temp$eval_imputed,
                  type = "response")
  species_data_temp$eval$rf_dib_pred_prob <- pred[["predictions"]]
  ##Prediction with rf_clim
  pred <- predict(species_data_temp$rforest_clim,
                  data = species_data_temp$eval_imputed,
                  type = "response")
  species_data_temp$eval$rf_clim_pred_prob <- pred[["predictions"]]
  ##Prediction with rf_dib_clim
  pred <- predict(species_data_temp$rforest_dib_clim,
                  data = species_data_temp$eval_imputed,
                  type = "response")
  species_data_temp$eval$rf_dib_clim_pred_prob <- pred[["predictions"]]
  assign(myspecies[b], species_data_temp)
  print(myspecies[b])
}
rm(pred)

##Validation via AUC
library(pROC)
for (c in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[c])
  species_data_temp["impact_metrics"] <- list(data.frame(model = c("glm_dib", "glm_clim","glm_dib_clim",
                                                                   "ct_dib", "ct_clim", "ct_dib_clim",
                                                                   "rf_dib", "rf_clim", "rf_dib_clim"),
                                                         AUC = c(1:9),
                                                         TSS_threshold = c(1:9),
                                                         TSS = c(1:9),
                                                         F1_threshold = c(1:9),
                                                         F1 = c(1:9)))
  ###auc_glm_dib
  roc <- roc(species_data_temp$eval$occurrenceStatus, 
             species_data_temp$eval$glm_dib_pred_prob)
  species_data_temp$impact_metrics[1,2] <- auc(roc)
  ###auc_glm_clim
  roc <- roc(species_data_temp$eval$occurrenceStatus, 
             species_data_temp$eval$glm_clim_pred_prob)
  species_data_temp$impact_metrics[2,2] <- auc(roc)
  ###auc_glm_dib_clim
  roc <- roc(species_data_temp$eval$occurrenceStatus, 
             species_data_temp$eval$glm_dib_clim_pred_prob)
  species_data_temp$impact_metrics[3,2] <- auc(roc)
  ###auc_ct_dib
  roc <- roc(species_data_temp$eval$occurrenceStatus, 
             species_data_temp$eval$ct_dib_pred_prob)
  species_data_temp$impact_metrics[4,2] <- auc(roc)
  ###auc_ct_clim
  roc <- roc(species_data_temp$eval$occurrenceStatus, 
             species_data_temp$eval$ct_clim_pred_prob)
  species_data_temp$impact_metrics[5,2] <- auc(roc)
  ###auc_ct_dib_clim
  roc <- roc(species_data_temp$eval$occurrenceStatus, 
             species_data_temp$eval$ct_dib_clim_pred_prob)
  species_data_temp$impact_metrics[6,2] <- auc(roc)
  ###auc_rf_dib
  roc <- roc(species_data_temp$eval$occurrenceStatus, 
             species_data_temp$eval$rf_dib_pred_prob)
  species_data_temp$impact_metrics[7,2] <- auc(roc)
  ###auc_rf_clim
  roc <- roc(species_data_temp$eval$occurrenceStatus, 
             species_data_temp$eval$rf_clim_pred_prob)
  species_data_temp$impact_metrics[8,2] <- auc(roc)
  ###auc_rf_dib_clim 
  roc <- roc(species_data_temp$eval$occurrenceStatus, 
             species_data_temp$eval$rf_dib_clim_pred_prob)
  species_data_temp$impact_metrics[9,2] <- auc(roc)
  assign(myspecies[c], species_data_temp)
  print(myspecies[c])
}

#thresholds via maximizing F1/TSS-Score
for (d in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[d])
  ##F1/TSS_glm_dib
  ###Defining vectors for tpr/fpr/fnr/tnr
  tpr <- c()
  fpr  <- c()
  fnr <- c()
  tnr <- c()
  for (dd in seq(0, 1, 0.005)) {
    ###Factorizing predictions and Occurrences 
    species_data_temp$eval$pred_class_glm_dib <- factor(ifelse(species_data_temp$eval$glm_dib_pred_prob > dd, 1, 0),
                                                        levels = c(1, 0))
    species_data_temp$eval$occurrenceStatus <- factor(species_data_temp$eval$occurrenceStatus, 
                                                      levels = c(1, 0))
    ###Creating confusion matrix
    pr <- species_data_temp$eval$pred_class_glm_dib
    ac <- species_data_temp$eval$occurrenceStatus
    kreuztabelle_t <- table(ac, pr)
    ###Calculating tpr/fpr/fnr/tnr
    tpr <- c(tpr, kreuztabelle_t[1, 1] / sum(kreuztabelle_t[1,]))
    fpr <- c(fpr, kreuztabelle_t[2, 1] / sum(kreuztabelle_t[2,]))
    tnr <- c(tnr, kreuztabelle_t[2, 2] / sum(kreuztabelle_t[2,]))
    fnr <- c(fnr, kreuztabelle_t[1, 2] / sum(kreuztabelle_t[1,]))
  }
  ###Creating dataframe with possible thresholds (p) and corresponding TSS-scores
  TSS_df <- data.frame(p = seq(0, 1, 0.005), TSS = tpr + tnr - 1)
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[1, 4] <- TSS_df[which.max(TSS_df$TSS), "TSS"]
  species_data_temp$impact_metrics[1, 3] <- TSS_df[which.max(TSS_df$TSS), "p"]
  p <- TSS_df[which.max(TSS_df$TSS), "p"]
  
  ###Creating dataframe with possible thresholds (p) and corresponding F1-scores
  f1_df <- data.frame(p = seq(0, 1, 0.005), f1 = tpr / (tpr + 0.5 * (fpr + fnr)))
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[1, 6] <- f1_df[which.max(f1_df$f1), "f1"]
  species_data_temp$impact_metrics[1, 5] <- f1_df[which.max(f1_df$f1), "p"]
  ###Adjusting factorisation for TSS score
  species_data_temp$eval$pred_class_glm_dib <- factor(ifelse(species_data_temp$eval$glm_dib_pred_prob > p, 1, 0),
                                                      levels = c(1, 0))
  ##F1/TSS_glm_clim
  ###Defining vectors for tpr/fpr/fnr/tnr
  tpr <- c()
  fpr  <- c()
  fnr <- c()
  tnr <- c()
  for (de in seq(0, 1, 0.005)) {
    ###Factorizing predictions and Occurrences 
    species_data_temp$eval$pred_class_glm_clim <- factor(ifelse(species_data_temp$eval$glm_clim_pred_prob > de, 1, 0),
                                                        levels = c(1, 0))
    species_data_temp$eval$occurrenceStatus <- factor(species_data_temp$eval$occurrenceStatus, 
                                                      levels = c(1, 0))
    ###Creating confusion matrix
    pr <- species_data_temp$eval$pred_class_glm_clim
    ac <- species_data_temp$eval$occurrenceStatus
    kreuztabelle_t <- table(ac, pr)
    ###Calculating tpr/fpr/fnr/tnr
    tpr <- c(tpr, kreuztabelle_t[1, 1] / sum(kreuztabelle_t[1,]))
    fpr <- c(fpr, kreuztabelle_t[2, 1] / sum(kreuztabelle_t[2,]))
    tnr <- c(tnr, kreuztabelle_t[2, 2] / sum(kreuztabelle_t[2,]))
    fnr <- c(fnr, kreuztabelle_t[1, 2] / sum(kreuztabelle_t[1,]))
  }
  ###Creating dataframe with possible thresholds (p) and corresponding f1-scores
  TSS_df <- data.frame(p = seq(0, 1, 0.005), TSS = tpr + tnr - 1)
  ###Extracting p and TSS for max(TSS)
  species_data_temp$impact_metrics[2, 4] <- TSS_df[which.max(TSS_df$TSS), "TSS"]
  species_data_temp$impact_metrics[2, 3] <- TSS_df[which.max(TSS_df$TSS), "p"]
  p <- TSS_df[which.max(TSS_df$TSS), "p"]
  
  ###Creating dataframe with possible thresholds (p) and corresponding F1-scores
  f1_df <- data.frame(p = seq(0, 1, 0.005), f1 = tpr / (tpr + 0.5 * (fpr + fnr)))
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[2, 6] <- f1_df[which.max(f1_df$f1), "f1"]
  species_data_temp$impact_metrics[2, 5] <- f1_df[which.max(f1_df$f1), "p"]
  ###Adjusting factorisation for TSS score
  species_data_temp$eval$pred_class_glm_clim <- factor(ifelse(species_data_temp$eval$glm_clim_pred_prob > p, 1, 0),
                                                      levels = c(1, 0))
  ##F1/TSS_glm_dib_clim
  ###Defining vectors for tpr/fpr/fnr/tnr
  tpr <- c()
  fpr  <- c()
  fnr <- c()
  tnr <- c()
  for (df in seq(0, 1, 0.005)) {
    ###Factorizing predictions and Occurrences 
    species_data_temp$eval$pred_class_glm_dib_clim <- factor(ifelse(species_data_temp$eval$glm_dib_clim_pred_prob > df, 1, 0),
                                                        levels = c(1, 0))
    species_data_temp$eval$occurrenceStatus <- factor(species_data_temp$eval$occurrenceStatus, 
                                                      levels = c(1, 0))
    ###Creating confusion matrix
    pr <- species_data_temp$eval$pred_class_glm_dib_clim
    ac <- species_data_temp$eval$occurrenceStatus
    kreuztabelle_t <- table(ac, pr)
    ###Calculating tpr/fpr/fnr/tnr
    tpr <- c(tpr, kreuztabelle_t[1, 1] / sum(kreuztabelle_t[1,]))
    fpr <- c(fpr, kreuztabelle_t[2, 1] / sum(kreuztabelle_t[2,]))
    tnr <- c(tnr, kreuztabelle_t[2, 2] / sum(kreuztabelle_t[2,]))
    fnr <- c(fnr, kreuztabelle_t[1, 2] / sum(kreuztabelle_t[1,]))
  }
  ###Creating dataframe with possible thresholds (p) and corresponding TSS-scores
  TSS_df <- data.frame(p = seq(0, 1, 0.005), TSS = tpr + tnr - 1)
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[3, 4] <- TSS_df[which.max(TSS_df$TSS), "TSS"]
  species_data_temp$impact_metrics[3, 3] <- TSS_df[which.max(TSS_df$TSS), "p"]
  p <- TSS_df[which.max(TSS_df$TSS), "p"]
  
  ###Creating dataframe with possible thresholds (p) and corresponding F1-scores
  f1_df <- data.frame(p = seq(0, 1, 0.005), f1 = tpr / (tpr + 0.5 * (fpr + fnr)))
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[3, 6] <- f1_df[which.max(f1_df$f1), "f1"]
  species_data_temp$impact_metrics[3, 5] <- f1_df[which.max(f1_df$f1), "p"]
  ###Adjusting factorisation for TSS score
  species_data_temp$eval$pred_class_glm_dib_clim <- factor(ifelse(species_data_temp$eval$glm_dib_clim_pred_prob > p, 1, 0),
                                                      levels = c(1, 0))
  ##TSS_ct_dib
  ###Defining vectors for tpr/fpr/fnr/tnr
  tpr <- c()
  fpr  <- c()
  fnr <- c()
  tnr <- c()
  for (dg in seq(0, 1, 0.005)) {
    ###Factorizing predictions and Occurrences 
    species_data_temp$eval$pred_class_ct_dib <- factor(ifelse(species_data_temp$eval$ct_dib_pred_prob > dg, 1, 0),
                                                             levels = c(1, 0))
    species_data_temp$eval$occurrenceStatus <- factor(species_data_temp$eval$occurrenceStatus, 
                                                      levels = c(1, 0))
    ###Creating confusion matrix
    pr <- species_data_temp$eval$pred_class_ct_dib
    ac <- species_data_temp$eval$occurrenceStatus
    kreuztabelle_t <- table(ac, pr)
    ###Calculating tpr/fpr/fnr/tnr
    tpr <- c(tpr, kreuztabelle_t[1, 1] / sum(kreuztabelle_t[1,]))
    fpr <- c(fpr, kreuztabelle_t[2, 1] / sum(kreuztabelle_t[2,]))
    tnr <- c(tnr, kreuztabelle_t[2, 2] / sum(kreuztabelle_t[2,]))
    fnr <- c(fnr, kreuztabelle_t[1, 2] / sum(kreuztabelle_t[1,]))
  }
  ###Creating dataframe with possible thresholds (p) and corresponding TSS-scores
  TSS_df <- data.frame(p = seq(0, 1, 0.005), TSS = tpr + tnr - 1)
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[4, 4] <- TSS_df[which.max(TSS_df$TSS), "TSS"]
  species_data_temp$impact_metrics[4, 3] <- TSS_df[which.max(TSS_df$TSS), "p"]
  p <- TSS_df[which.max(TSS_df$TSS), "p"]
  
  ###Creating dataframe with possible thresholds (p) and corresponding F1-scores
  f1_df <- data.frame(p = seq(0, 1, 0.005), f1 = tpr / (tpr + 0.5 * (fpr + fnr)))
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[4, 6] <- f1_df[which.max(f1_df$f1), "f1"]
  species_data_temp$impact_metrics[4, 5] <- f1_df[which.max(f1_df$f1), "p"]
  ###Adjusting factorisation for TSS score
  species_data_temp$eval$pred_class_ct_dib <- factor(ifelse(species_data_temp$eval$ct_dib_pred_prob > p, 1, 0),
                                                           levels = c(1, 0))
  ##TSS_ct_clim
  ###Defining vectors for tpr/fpr/fnr/tnr
  tpr <- c()
  fpr  <- c()
  fnr <- c()
  tnr <- c()
  for (dh in seq(0, 1, 0.005)) {
    ###Factorizing predictions and Occurrences 
    species_data_temp$eval$pred_class_ct_clim <- factor(ifelse(species_data_temp$eval$ct_clim_pred_prob > dh, 1, 0),
                                                       levels = c(1, 0))
    species_data_temp$eval$occurrenceStatus <- factor(species_data_temp$eval$occurrenceStatus, 
                                                      levels = c(1, 0))
    ###Creating confusion matrix
    pr <- species_data_temp$eval$pred_class_ct_clim
    ac <- species_data_temp$eval$occurrenceStatus
    kreuztabelle_t <- table(ac, pr)
    ###Calculating tpr/fpr/fnr/tnr
    tpr <- c(tpr, kreuztabelle_t[1, 1] / sum(kreuztabelle_t[1,]))
    fpr <- c(fpr, kreuztabelle_t[2, 1] / sum(kreuztabelle_t[2,]))
    tnr <- c(tnr, kreuztabelle_t[2, 2] / sum(kreuztabelle_t[2,]))
    fnr <- c(fnr, kreuztabelle_t[1, 2] / sum(kreuztabelle_t[1,]))
  }
  ###Creating dataframe with possible thresholds (p) and corresponding TSS-scores
  TSS_df <- data.frame(p = seq(0, 1, 0.005), TSS = tpr + tnr - 1)
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[5, 4] <- TSS_df[which.max(TSS_df$TSS), "TSS"]
  species_data_temp$impact_metrics[5, 3] <- TSS_df[which.max(TSS_df$TSS), "p"]
  p <- TSS_df[which.max(TSS_df$TSS), "p"]
  
  ###Creating dataframe with possible thresholds (p) and corresponding F1-scores
  f1_df <- data.frame(p = seq(0, 1, 0.005), f1 = tpr / (tpr + 0.5 * (fpr + fnr)))
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[5, 6] <- f1_df[which.max(f1_df$f1), "f1"]
  species_data_temp$impact_metrics[5, 5] <- f1_df[which.max(f1_df$f1), "p"]
  ###Adjusting factorisation for TSS score
  species_data_temp$eval$pred_class_ct_clim <- factor(ifelse(species_data_temp$eval$ct_clim_pred_prob > p, 1, 0),
                                                     levels = c(1, 0))
  ##TSS_ct_dib_clim
  ###Defining vectors for tpr/fpr/fnr/tnr
  tpr <- c()
  fpr  <- c()
  fnr <- c()
  tnr <- c()
  for (di in seq(0, 1, 0.005)) {
    ###Factorizing predictions and Occurrences 
    species_data_temp$eval$pred_class_ct_dib_clim <- factor(ifelse(species_data_temp$eval$ct_dib_clim_pred_prob > di, 1, 0),
                                                       levels = c(1, 0))
    species_data_temp$eval$occurrenceStatus <- factor(species_data_temp$eval$occurrenceStatus, 
                                                      levels = c(1, 0))
    ###Creating confusion matrix
    pr <- species_data_temp$eval$pred_class_ct_dib_clim
    ac <- species_data_temp$eval$occurrenceStatus
    kreuztabelle_t <- table(ac, pr)
    ###Calculating tpr/fpr/fnr/tnr
    tpr <- c(tpr, kreuztabelle_t[1, 1] / sum(kreuztabelle_t[1,]))
    fpr <- c(fpr, kreuztabelle_t[2, 1] / sum(kreuztabelle_t[2,]))
    tnr <- c(tnr, kreuztabelle_t[2, 2] / sum(kreuztabelle_t[2,]))
    fnr <- c(fnr, kreuztabelle_t[1, 2] / sum(kreuztabelle_t[1,]))
  }
  ###Creating dataframe with possible thresholds (p) and corresponding TSS-scores
  TSS_df <- data.frame(p = seq(0, 1, 0.005), TSS = tpr + tnr - 1)
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[6, 4] <- TSS_df[which.max(TSS_df$TSS), "TSS"]
  species_data_temp$impact_metrics[6, 3] <- TSS_df[which.max(TSS_df$TSS), "p"]
  p <- TSS_df[which.max(TSS_df$TSS), "p"]
  
  ###Creating dataframe with possible thresholds (p) and corresponding F1-scores
  f1_df <- data.frame(p = seq(0, 1, 0.005), f1 = tpr / (tpr + 0.5 * (fpr + fnr)))
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[6, 6] <- f1_df[which.max(f1_df$f1), "f1"]
  species_data_temp$impact_metrics[6, 5] <- f1_df[which.max(f1_df$f1), "p"]
  ###Adjusting factorisation for TSS score
  species_data_temp$eval$pred_class_ct_dib_clim <- factor(ifelse(species_data_temp$eval$ct_dib_clim_pred_prob > p, 1, 0),
                                                     levels = c(1, 0))
  ##TSS_rf_dib
  ###Defining vectors for tpr/fpr/fnr/tnr
  tpr <- c()
  fpr  <- c()
  fnr <- c()
  tnr <- c()
  for (dj in seq(0, 1, 0.005)) {
    ###Factorizing predictions and Occurrences 
    species_data_temp$eval$pred_class_rf_dib <- factor(ifelse(species_data_temp$eval$rf_dib_pred_prob > dj, 1, 0),
                                                       levels = c(1, 0))
    species_data_temp$eval$occurrenceStatus <- factor(species_data_temp$eval$occurrenceStatus, 
                                                      levels = c(1, 0))
    ###Creating confusion matrix
    pr <- species_data_temp$eval$pred_class_rf_dib
    ac <- species_data_temp$eval$occurrenceStatus
    kreuztabelle_t <- table(ac, pr)
    ###Calculating tpr/fpr/fnr/tnr
    tpr <- c(tpr, kreuztabelle_t[1, 1] / sum(kreuztabelle_t[1,]))
    fpr <- c(fpr, kreuztabelle_t[2, 1] / sum(kreuztabelle_t[2,]))
    tnr <- c(tnr, kreuztabelle_t[2, 2] / sum(kreuztabelle_t[2,]))
    fnr <- c(fnr, kreuztabelle_t[1, 2] / sum(kreuztabelle_t[1,]))
  }
  ###Creating dataframe with possible thresholds (p) and corresponding TSS-scores
  TSS_df <- data.frame(p = seq(0, 1, 0.005), TSS = tpr + tnr - 1)
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[7, 4] <- TSS_df[which.max(TSS_df$TSS), "TSS"]
  species_data_temp$impact_metrics[7, 3] <- TSS_df[which.max(TSS_df$TSS), "p"]
  p <- TSS_df[which.max(TSS_df$TSS), "p"]
  
  ###Creating dataframe with possible thresholds (p) and corresponding F1-scores
  f1_df <- data.frame(p = seq(0, 1, 0.005), f1 = tpr / (tpr + 0.5 * (fpr + fnr)))
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[7, 6] <- f1_df[which.max(f1_df$f1), "f1"]
  species_data_temp$impact_metrics[7, 5] <- f1_df[which.max(f1_df$f1), "p"]
  ###Adjusting factorisation for TSS score
  species_data_temp$eval$pred_class_rf_dib <- factor(ifelse(species_data_temp$eval$rf_dib_pred_prob > p, 1, 0),
                                                     levels = c(1, 0))
  ##TSS_rf_clim
  ###Defining vectors for tpr/fpr/fnr/tnr
  tpr <- c()
  fpr  <- c()
  fnr <- c()
  tnr <- c()
  for (dk in seq(0, 1, 0.005)) {
    ###Factorizing predictions and Occurrences 
    species_data_temp$eval$pred_class_rf_clim <- factor(ifelse(species_data_temp$eval$rf_clim_pred_prob > dk, 1, 0),
                                                       levels = c(1, 0))
    species_data_temp$eval$occurrenceStatus <- factor(species_data_temp$eval$occurrenceStatus, 
                                                      levels = c(1, 0))
    ###Creating confusion matrix
    pr <- species_data_temp$eval$pred_class_rf_clim
    ac <- species_data_temp$eval$occurrenceStatus
    kreuztabelle_t <- table(ac, pr)
    ###Calculating tpr/fpr/fnr/tnr
    tpr <- c(tpr, kreuztabelle_t[1, 1] / sum(kreuztabelle_t[1,]))
    fpr <- c(fpr, kreuztabelle_t[2, 1] / sum(kreuztabelle_t[2,]))
    tnr <- c(tnr, kreuztabelle_t[2, 2] / sum(kreuztabelle_t[2,]))
    fnr <- c(fnr, kreuztabelle_t[1, 2] / sum(kreuztabelle_t[1,]))
  }
  ###Creating dataframe with possible thresholds (p) and corresponding TSS-scores
  TSS_df <- data.frame(p = seq(0, 1, 0.005), TSS = tpr + tnr - 1)
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[8, 4] <- TSS_df[which.max(TSS_df$TSS), "TSS"]
  species_data_temp$impact_metrics[8, 3] <- TSS_df[which.max(TSS_df$TSS), "p"]
  p <- TSS_df[which.max(TSS_df$TSS), "p"]
  
  ###Creating dataframe with possible thresholds (p) and corresponding F1-scores
  f1_df <- data.frame(p = seq(0, 1, 0.005), f1 = tpr / (tpr + 0.5 * (fpr + fnr)))
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[8, 6] <- f1_df[which.max(f1_df$f1), "f1"]
  species_data_temp$impact_metrics[8, 5] <- f1_df[which.max(f1_df$f1), "p"]
  ###Adjusting factorisation for TSS score
  species_data_temp$eval$pred_class_rf_clim <- factor(ifelse(species_data_temp$eval$rf_clim_pred_prob > p, 1, 0),
                                                     levels = c(1, 0))
  ##TSS_rf_dib_clim
  ###Defining vectors for tpr/fpr/fnr/tnr
  tpr <- c()
  fpr  <- c()
  fnr <- c()
  tnr <- c()
  for (dl in seq(0, 1, 0.005)) {
    ###Factorizing predictions and Occurrences 
    species_data_temp$eval$pred_class_rf_dib_clim <- factor(ifelse(species_data_temp$eval$rf_dib_clim_pred_prob > dl, 1, 0),
                                                       levels = c(1, 0))
    species_data_temp$eval$occurrenceStatus <- factor(species_data_temp$eval$occurrenceStatus, 
                                                      levels = c(1, 0))
    ###Creating confusion matrix
    pr <- species_data_temp$eval$pred_class_rf_dib_clim
    ac <- species_data_temp$eval$occurrenceStatus
    kreuztabelle_t <- table(ac, pr)
    ###Calculating tpr/fpr/fnr/tnr
    tpr <- c(tpr, kreuztabelle_t[1, 1] / sum(kreuztabelle_t[1,]))
    fpr <- c(fpr, kreuztabelle_t[2, 1] / sum(kreuztabelle_t[2,]))
    tnr <- c(tnr, kreuztabelle_t[2, 2] / sum(kreuztabelle_t[2,]))
    fnr <- c(fnr, kreuztabelle_t[1, 2] / sum(kreuztabelle_t[1,]))
  }
  ###Creating dataframe with possible thresholds (p) and corresponding TSS-scores
  TSS_df <- data.frame(p = seq(0, 1, 0.005), TSS = tpr + tnr - 1)
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[9, 4] <- TSS_df[which.max(TSS_df$TSS), "TSS"]
  species_data_temp$impact_metrics[9, 3] <- TSS_df[which.max(TSS_df$TSS), "p"]
  p <- TSS_df[which.max(TSS_df$TSS), "p"]
  
  ###Creating dataframe with possible thresholds (p) and corresponding F1-scores
  f1_df <- data.frame(p = seq(0, 1, 0.005), f1 = tpr / (tpr + 0.5 * (fpr + fnr)))
  ###Extracting p and f1 for max(f1)
  species_data_temp$impact_metrics[9, 6] <- f1_df[which.max(f1_df$f1), "f1"]
  species_data_temp$impact_metrics[9, 5] <- f1_df[which.max(f1_df$f1), "p"]
  ###Adjusting factorisation for TSS score
  species_data_temp$eval$pred_class_rf_dib_clim <- factor(ifelse(species_data_temp$eval$rf_dib_clim_pred_prob > p, 1, 0),
                                                     levels = c(1, 0))
  assign(myspecies[d], species_data_temp)
  print(myspecies[d])
}
rm(d, dd, de, df, dg, dh, di, dj, dk, dl,  p, f1, f1_df, tpr, fpr, fnr, tnr, kreuztabelle_t, pr, ac)


#Creating dataframes with all AUC and TSS Scores
TSS_overview <- data.frame(species = c(1:21),
                           glm_dib = c(1:21),
                           glm_clim = c(1:21),
                           glm_dib_clim = c(1:21),
                           ct_dib = c(1:21),
                           ct_clim = c(1:21),
                           ct_dib_clim = c(1:21),
                           rf_dib = c(1:21),
                           rf_clim = c(1:21),
                           rf_dib_clim = c(1:21))
AUC_overview <- data.frame(species = c(1:21),
                           glm_dib = c(1:21),
                           glm_clim = c(1:21),
                           glm_dib_clim = c(1:21),
                           ct_dib = c(1:21),
                           ct_clim = c(1:21),
                           ct_dib_clim = c(1:21),
                           rf_dib = c(1:21),
                           rf_clim = c(1:21),
                           rf_dib_clim = c(1:21))

for (e in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[e])
  TSS_overview[e, 1] <- myspecies[e]
  TSS_overview[e, 2:10] <- species_data_temp$impact_metrics$TSS
  AUC_overview[e, 1] <- myspecies[e]
  AUC_overview[e, 2:10] <- species_data_temp$impact_metrics$AUC
  print(myspecies[e])
}

setwd("C:/Users/luca4/OneDrive/UniUnterlagen/TUM/7.Semester/BA/Results")

write.csv(file = "TSS.csv", TSS_overview)
write.csv(file = "AUC.csv", AUC_overview)

#Prediciting maps for RForest models
library(ranger)
library(dplyr)
library(terra)

#Creating prediction maps
for (f in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[f])
  
  disturbance_df <- terra::as.data.frame(disturbance, xy = TRUE, cells = TRUE)
  species_data_temp$disturbance_pred <- na.omit(disturbance_df)
  species_data_temp$disturbance_pred$pred_dib <- predict(species_data_temp$rforest_dib,
                                                     species_data_temp$disturbance_pred,
                                                     type = "response")[["predictions"]]
  
  bioclim_df <- terra::as.data.frame(bioclim, xy = TRUE, cells = TRUE)
  species_data_temp$bioclim_pred <- na.omit(bioclim_df)
  species_data_temp$bioclim_pred$pred_clim <- predict(species_data_temp$rforest_clim, 
                                                 species_data_temp$bioclim_pred,
                                                 type = "response")[["predictions"]]
  
  dib_clim_fc_df <- terra::as.data.frame(dib_clim_fc, xy = TRUE, cells = TRUE)
  species_data_temp$dib_clim_pred <- na.omit(dib_clim_fc_df)
  species_data_temp$dib_clim_pred$pred_dib_clim <- predict(species_data_temp$rforest_dib_clim, 
                                                  species_data_temp$dib_clim_pred,
                                                  type = "response")[["predictions"]]
  
  preds <- species_data_temp$bioclim_pred[c("cell", "pred_clim")]
  species_data_temp$disturbance_pred <- left_join(species_data_temp$disturbance_pred,
                                                            preds,
                                                            by = c("cell"))
  
  preds <- species_data_temp$dib_clim_pred[c("cell", "pred_dib_clim")]
  species_data_temp$disturbance_pred <- left_join(species_data_temp$disturbance_pred,
                                                                preds,
                                                                by = c("cell"))
  
  #write.csv(file = paste0("prediction_maps/species_", f, ".csv"), species_data_temp$disturbance_pred)
  rm(preds)
  assign(myspecies[f], species_data_temp)
  print(myspecies[f])
}

#Creating prediction maps
for (g in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[g])
  ##create raster for dib
  preds <- species_data_temp$disturbance_pred[c("x", "y", "pred_dib")]
  raster <- rast(preds, type = "xyz") 
  writeRaster(raster, paste0("prediction_maps/pred_dib/species", g, ".tif"))
  ##create raster for clim
  preds <- species_data_temp$disturbance_pred[c("x", "y", "pred_clim")]
  raster <- rast(preds, type = "xyz") 
  writeRaster(raster, paste0("prediction_maps/pred_clim/species", g, ".tif"))
  ##create raster for dib_clim
  preds <- species_data_temp$disturbance_pred[c("x", "y", "pred_dib_clim")]
  raster <- rast(preds, type = "xyz") 
  writeRaster(raster, paste0("prediction_maps/pred_dib_clim/species", g, ".tif"))
  ##create diff-raster
  species_data_temp$disturbance_pred$delta <- (abs(species_data_temp$disturbance_pred$pred_dib_clim) - abs(species_data_temp$disturbance_pred$pred_clim))
  preds <- species_data_temp$disturbance_pred[c("x", "y", "delta")]
  raster <- rast(preds, type = "xyz") 
  writeRaster(raster, paste0("prediction_maps/deltas/species", g, ".tif"))
  print(myspecies[g])
}

