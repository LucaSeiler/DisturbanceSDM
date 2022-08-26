library(sp)
library(terra)
library(dplyr)
library(mapview)
library(ggplot2)
library(tictoc)
set.seed(420)

#Generating pseudoabsences
for (a in 1:length(myspecies)) {
  ##Random-Sampling of pseudoabsences
  species_data_temp <- get(myspecies[a])
  pseudoabsences <- spatSample(disturbance$dib1, 
                               size = nrow(species_data_temp),
                               method = "random",
                               as.df = TRUE,
                               xy = TRUE,
                               values = FALSE,
                               na.rm = TRUE)
  pseudoabsences <- as.data.frame(pseudoabsences)
  ##Adding climate and disturbance parameters
  pseudoabsences$ID <- c(1:nrow(pseudoabsences))
  ###Vectorizing Data
  pseudoabsences_vec <- vect(pseudoabsences, 
                                geom=c("x", "y"), 
                                crs="+init=epsg:4326")
  ###Extracting Predictor-data for Pseudoabsence point
  pseudoabsences_dib_clim_fc <- terra::extract(dib_clim_fc, 
                                               pseudoabsences_vec)
  ###Joining Predictor data and occurrences
  pseudoabsences <- full_join(pseudoabsences_dib_clim_fc, 
                              pseudoabsences, 
                              by="ID")
  ##Dummy-Coding
  species_data_temp$occurrenceStatus <- 1
  pseudoabsences$occurrenceStatus <- 0
  ##Joining pseudoabsences and presences
  species_data_temp <- bind_rows(species_data_temp, 
                                 pseudoabsences)
  assign(myspecies[a], species_data_temp)
  print(myspecies[a])
}

#Plotting ecological niches
ggplot(data = get(myspecies[3])) +
  geom_boxplot(aes(x = occurrenceStatus, 
                   y = dib4, 
                   col = factor(occurrenceStatus)))


#Splitting Evaluation-Datasets
for (b in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[b])
  training <- sample(nrow(species_data_temp), 
                     0.8*nrow(species_data_temp))
  train <- species_data_temp[training, ]
  eval <- species_data_temp[-training, ]
  ##Creating list of datasets
  list <- list(train, eval)
  nms <- c("train", "eval")
  names(list) = nms
  assign(myspecies[b], list)
  print(myspecies[b])
}
rm(b, list, train, nms, eval, training, species_data_temp)

#Fitting
##GLM
##Testing correlation
kol <- cor(`Anastrangalia sanguinolenta`$train[,2:24])

##Fitting glms with disturbance and climate parameters
for (c in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[c])
  fit <- glm(occurrenceStatus ~ (bio1 + bio7 + bio9 + bio12 + bio15 + bio17) * 
               (dib1 + dib2 + dib3 + dib4) *
               forestcover,
             data = species_data_temp$train,
             family = "binomial")
  ###Creating lists of datasets and glms
  species_data_temp["glm_dib_clim"] <- list(fit)
  assign(myspecies[c], species_data_temp)
  print(myspecies[c])
}

##Fitting glms with only climate parameters
for (d in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[d])
  fit <- glm(occurrenceStatus ~ (bio1 + bio7 + bio9 + bio12 + bio15 + bio17),
             data = species_data_temp$train,
             family = "binomial")
  ###Creating lists of datasets and glms
  species_data_temp["glm_clim"] <- list(fit)
  assign(myspecies[d], species_data_temp)
  print(myspecies[d])
}

##Fitting glms with only disturbance parameters
for (d in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[d])
  fit <- glm(occurrenceStatus ~ (dib1 * dib2 * dib3 * dib4),
             data = species_data_temp$train,
             family = "binomial")
  ###Creating lists of datasets and glms
  species_data_temp["glm_dib"] <- list(fit)
  assign(myspecies[d], species_data_temp)
  print(myspecies[d])
}

##Classification Tree
library(rpart) 
library(tictoc)

###Building trees with disturbance and climate parameters
####Defining the formula 
formula1 <- as.formula(paste("occurrenceStatus ~ ", 
                            paste(colnames(species_data_temp$train[c(2:25)]),
                                  collapse = " + ")))

for (e in 1:length(myspecies)) {
  ####Building initial tree
  species_data_temp <- get(myspecies[e])
  tic("growing")
  tree <- rpart(formula1,
                data = species_data_temp$train,
                method = "class",
                #min number of observations per split = 20
                #min number of obersvations per leaf
                control = rpart.control(cp = .0001))  
  toc()
  ####Pruning tree via optimizing the complexity parameter by minimizing xerror (cross-validatet)
  best <- tree$cptable[which.min(tree$cptable[,"xerror"]),"CP"]
  pruned_tree <- prune(tree, cp=best)
  ####Adding to list
  species_data_temp["pruned_tree_dib_clim"] <- list(pruned_tree)
  assign(myspecies[e], species_data_temp)
  print(myspecies[e])
}

###Building trees with only climate parameters
####Defining the formula 
formula2 <- as.formula(paste("occurrenceStatus ~ ", 
                             paste(colnames(species_data_temp$train[c(6:24)]),
                                   collapse = " + ")))

for (f in 1:length(myspecies)) {
  ####Building initial tree
  species_data_temp <- get(myspecies[f])
  tic("growing")
  tree <- rpart(formula2,
                data = species_data_temp$train,
                method = "class",
                control = rpart.control(cp = .0001))  
  toc()
  ####Pruning tree via xerror
  best <- tree$cptable[which.min(tree$cptable[,"xerror"]),"CP"]
  pruned_tree <- prune(tree, cp=best)
  ####Adding to list
  species_data_temp["pruned_tree_clim"] <- list(pruned_tree)
  assign(myspecies[f], species_data_temp)
  print(myspecies[f])
}

###Building trees with only disturbance parameters
####Defining the formula 
formula3 <- as.formula(paste("occurrenceStatus ~ ", 
                             paste(colnames(species_data_temp$train[c(2:5)]),
                                   collapse = " + ")))
for (f in 1:length(myspecies)) {
  ####Building initial tree
  species_data_temp <- get(myspecies[f])
  tic("growing")
  tree <- rpart(formula3,
                data = species_data_temp$train,
                method = "class",
                control = rpart.control(cp = .0001))  
  toc()
  ####Pruning tree via xerror
  best <- tree$cptable[which.min(tree$cptable[,"xerror"]),"CP"]
  pruned_tree <- prune(tree, cp=best)
  ####Adding to list
  species_data_temp["pruned_tree_dib"] <- list(pruned_tree)
  assign(myspecies[f], species_data_temp)
  print(myspecies[f])
}


##RandomForest
library(ranger)
library(randomForest)

###Growing forest with disturbance and climate parameters
####Defining the formula 
formula1 <- as.formula(paste("occurrenceStatus ~ ", 
                             paste(colnames(species_data_temp$train[c(2:25)]),
                                   collapse = " + ")))

####Cleaning out NAs in zraining dataset
for (g in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[g])
  dib_NA_rows_all <- which(is.na(species_data_temp$train$dib1) | 
                             is.na(species_data_temp$train$dib2) | 
                             is.na(species_data_temp$train$dib3) | 
                             is.na(species_data_temp$train$dib4))
  length(dib_NA_rows_all)
  if (length(dib_NA_rows_all) > 0) {
    species_data_temp["train_rf"] <- list(species_data_temp$train[-dib_NA_rows_all, ])
  }
  assign(myspecies[g] ,species_data_temp)
  print(myspecies[g])
}


####Growing the forest
for (h in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[h])
  rforest <- ranger(formula1,
                    data = species_data_temp$train_rf,
                    num.trees = 500,
                    mtry = floor(sqrt(24)),
                    respect.unordered.factors = "order",
                    importance = "impurity",
                    verbose = TRUE,
                    seed = 420)
  ####Adding to dataset
  species_data_temp["rforest_dib_clim"] <- list(rforest)
  assign(myspecies[h], species_data_temp)
  print(myspecies[h])
}


###Growing forest with only climate parameters
####Defining the formula
formula2 <- as.formula(paste("occurrenceStatus ~ ", 
                             paste(colnames(species_data_temp$train[c(6:24)]),
                                   collapse = " + ")))
####Growing the forest
for (i in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[i])
  rforest <- ranger(formula2,
                    data = species_data_temp$train_rf,
                    num.trees = 500,
                    mtry = floor(sqrt(24)),
                    respect.unordered.factors = "order",
                    importance = "impurity",
                    verbose = TRUE,
                    seed = 420)
  ####Adding to dataset
  species_data_temp["rforest_clim"] <- list(rforest)
  assign(myspecies[i], species_data_temp)
  print(myspecies[i])
}

###Growing forest with only disturbance parameters
####Defining the formula
formula3 <- as.formula(paste("occurrenceStatus ~ ", 
                             paste(colnames(species_data_temp$train[c(2:5)]),
                                   collapse = " + ")))
####Growing the forest
for (i in 1:length(myspecies)) {
  species_data_temp <- get(myspecies[i])
  rforest <- ranger(formula3,
                    data = species_data_temp$train_rf,
                    num.trees = 500,
                    mtry = floor(sqrt(24)),
                    respect.unordered.factors = "order",
                    importance = "impurity",
                    verbose = TRUE,
                    seed = 420)
  ####Adding to dataset
  species_data_temp["rforest_dib"] <- list(rforest)
  assign(myspecies[i], species_data_temp)
  print(myspecies[i])
}
