## ---------------------------
##
## Script name: backwardsElimination.R
##
## Author: Mr. Ashleigh C. Myall
##
## Date Created: 01-11-2020
##
## Copyright (c) Ashleigh C. Myall 2020
## Email: a.myall19@imperial.ac.uk
##
## ---------------------------
##
## Notes: Main script to run feature selection (backwards elimination) over gene 
## expression dataLoad in data. Run backwards elimination.
##
## ---------------------------

################################################################################

### Load Libs

################################################################################

library(parallel)
library(ranger)
library(caret)
library(varSelRF)


################################################################################

###
###                               Functions
###

################################################################################

# ------------------------------------------------------------------------------


#' Parallel Backwards Elimination function
#'
#' @description Parallel implemenation of backwards elimination with the ranger 
#' random forest function. Take an input of data list, runs backwards elimination 
#' over the train, and test split, then reports optimal model performance on val data.
#'
#' @param datapath A character of input file location (RDS list of train, test, val data).
#' @param name A character for output files
#' @param trees Number of trees for random forest
#' @param dropFrac Percentage of variables to drop at each step in backwards elinimation
#' @param runs number of repeitions of backwards elinimation to run
#' @param featureSet optional arguement to run over a subset of features
#' @param cores number of cores to run over in parallel
#' @return List of run outputs, optimal model performance, gene name mapping to model
## ---------------------------
##
## Script name: backwardsElimination.R
##
## Author: Mr. Ashleigh C. Myall
##
## Date Created: 01-11-2020
##
## Copyright (c) Ashleigh C. Myall 2020
## Email: a.myall19@imperial.ac.uk
##
## ---------------------------
##
## Notes: Main script to run feature selection (backwards elimination) over gene 
## expression dataLoad in data. Run backwards elimination.
##
## ---------------------------

################################################################################

### Load Libs

################################################################################

library(parallel)
library(ranger)
library(caret)
library(varSelRF)
library(dplyr)
library(ggplot2)


################################################################################

###
###                               Functions
###

################################################################################

# ------------------------------------------------------------------------------


#' Parallel Backwards Elimination function
#'
#' @description Parallel implemenation of backwards elimination with the ranger 
#' random forest function. Take an input of data list, runs backwards elimination 
#' over the train, and test split, then reports optimal model performance on val data.
#'
#' @param datapath A character of input file location (RDS list of train, test, val data).
#' @param name A character for output files
#' @param trees Number of trees for random forest
#' @param dropFrac Percentage of variables to drop at each step in backwards elinimation
#' @param runs number of repeitions of backwards elinimation to run
#' @param featureSet optional arguement to run over a subset of features
#' @param cores number of cores to run over in parallel
#' @return List of run outputs, optimal model performance, gene name mapping to model

doParBackElim <- function(data_path, name, trees = 5000, dropFrac = 0.001, runs = 100,
                        featureSet = NULL,cores = 40){
  
  data <- readRDS(data_path)
  
  train_dat <- data$train
  test_data <- data$test
  val_data <- data$val
  
  
  #Subset if freatureset not NULL
  if(!is.null(featureSet)){
    features <- readRDS(featureSet)
    train_dat <- train_dat[, which(colnames(train_dat) %in% c(features,"sample.type")), with = F]
    test_data <- test_data[, which(colnames(test_data) %in% c(features,"sample.type")), with = F]
    val_data <- val_data[, which(colnames(val_data) %in% c(features,"sample.type")), with = F]
  }
  
  
  col.lim <- (ncol(train_dat) - 1)
  
  # Format names to not interfere with random forest run
  name.mapping <- data.frame("orig" = colnames(train_dat)[1:col.lim], "new" = paste("Gene",1:col.lim,sep = ""))
  colnames(train_dat)[1:col.lim] <- paste("Gene",1:col.lim,sep = "")
  colnames(test_data)[1:col.lim] <- paste("Gene",1:col.lim,sep = "")
  colnames(val_data)[1:col.lim] <- paste("Gene",1:col.lim,sep = "")
  #saveRDS(name.mapping,file = paste("./Output/",name,".name.mapping.RDS",sep = ""))
  
  
  # Make list the length of intention
  myls <- vector("list", length = runs)
  names(myls) <- 1:runs
  
  
  print(paste("Doing: ",name,sep = ""))
  
  # Too count progress
  for (i in 1:length(myls)) {
    myls[[i]]$count <- i
  }
  
  save.list <- mclapply(myls, function(x){
    
    # Run varselRF
    backPs <- varSelRF2(xdata = train_dat[,1:col.lim],Class = train_dat$sample.type,
                        ntree = trees,
                        vars.drop.frac = dropFrac,
                        whole.rang = T,
                        verbose = T,
                        recompute.var.imp = T)
    
    # Calc class weights - to adjust random forest objective function.
    data.s <-rbind(data$train,data$test,data$val)
    cls.tab <- table(data.s$sample.type)
    cls.tab <- (1-(cls.tab/sum(cls.tab)))
    cls.weight <- as.numeric(cls.tab)
    names(cls.weight) <- names(cls.tab)
    
    # Run ranger (fast random forest) over optimised feature set to create final model.
    rfd <- ranger(data = subset(train_dat, select = c(backPs$selected.vars,"sample.type")),
                  dependent.variable.name = 'sample.type',
                  num.trees= trees,
                  class.weights = cls.weight)
    
    
    # Make prediction
    pred.test <- predict(rfd, subset(test_data, select = c(backPs$selected.vars,"sample.type")))
    pred.val <- predict(rfd, subset(val_data, select = c(backPs$selected.vars,"sample.type")))
    
    confus.test <- confusionMatrix(pred.test$predictions,test_data$sample.type,mode='everything')
    confus.val <- confusionMatrix(pred.val$predictions,val_data$sample.type,mode='everything')
    
    test.accuracy <- confus.test[["overall"]][1]
    val.accuracy <- confus.val[["overall"]][1]
    
    
    return.list <- list(
      "VarSelRF" = backPs,
      "Ranger" = rfd,
      "conf.test" = confus.test,
      "conf.val" = confus.val,
      "test.accuracy" = test.accuracy,
      "val.accuracy" = val.accuracy
    )
    
    
    return(return.list)
    
  },mc.cores = cores)
  
  saveRDS(save.list, file = paste("./Output/",name,"_backwards_elim.RDS",sep = ""))
  save.list$name.mapping = name.mapping
  return(save.list)
  
}


# ------------------------------------------------------------------------------

## Varself RF modified function to work with Ranger

varSelRF2 <- function (xdata, Class, c.sd = 1, mtryFactor = 1, ntree = 5000, 
                       ntreeIterat = 2000, vars.drop.num = NULL, vars.drop.frac = 0.2, 
                       whole.range = TRUE, recompute.var.imp = FALSE, verbose = FALSE, 
                       returnFirstForest = TRUE, fitted.rf = NULL, keep.forest = FALSE, rfpackage="ranger", invert.importance=F) 
{
  if (!is.factor(Class)) 
    stop("Class should be a factor")
  if ((is.null(vars.drop.num) & is.null(vars.drop.frac)) | 
      (!is.null(vars.drop.num) & !is.null(vars.drop.frac))) 
    stop("One (and only one) of vars.drop.frac and vars.drop.num must be NULL and the other set")
  max.num.steps <- dim(xdata)[2]
  num.subjects <- dim(xdata)[1]
  if (is.null(colnames(xdata))) 
    colnames(xdata) <- paste("v", 1:dim(xdata)[2], sep = "")
  n.vars <- vars <- OOB.rf <- OOB.sd <- rep(NA, max.num.steps)
  oobError <- function(rf) {
    if(class(rf) == "randomForest"){
      ooo <- rf$confusion[, -dim(rf$confusion)[2]]
    }else if(class(rf) == "ranger"){
      ooo <- rf$confusion
    }
    s.ooo <- sum(ooo)
    diag(ooo) <- 0
    sum(ooo)/s.ooo
  }
  data.df <- data.frame(xdata,class=Class)
  if (!is.null(fitted.rf)) {
    if (ncol(fitted.rf$importance) < 2) 
      stop("The fitted rf was not fitted with importance = TRUE")
    n.ntree <- fitted.rf$ntree
    mtry <- fitted.rf$mtry
    n.mtryFactor <- mtry/sqrt(ncol(xdata))
    if ((n.ntree != ntree) | (n.mtryFactor != mtryFactor)) 
      warning("Using as ntree and mtry the parameters obtained from fitted.rf", 
              immediate. = TRUE)
    ntree <- n.ntree
    mtryFactor <- n.mtryFactor
    rm(n.ntree, n.mtryFactor)
    rf <- fitted.rf
  } else {
    mtry <- floor(sqrt(ncol(xdata)) * mtryFactor)
    if(rfpackage == "ranger"){
      rf <- ranger(data = data.df, dependent.variable.name="class", num.trees = ntree, 
                   mtry = mtry, importance = "impurity", keep.inbag=T)
    }else{
      rf <- randomForest(x = xdata, y = Class, ntree = ntree, mtry = mtry, importance = TRUE, keep.forest = keep.forest)
    }
  }
  if (returnFirstForest){ 
    FirstForest <- rf
  }else{
    FirstForest <- NULL
  }
  m.iterated.ob.error <- m.initial.ob.error <- oobError(rf)
  sd.iterated.ob.error <- sd.initial.ob.error <- sqrt(m.iterated.ob.error * 
                                                        (1 - m.iterated.ob.error) * (1/num.subjects))
  if (verbose) {
    print(paste("Initial OOB error: mean = ", round(m.initial.ob.error, 
                                                    4), "; sd = ", round(sd.initial.ob.error, 4), sep = ""))
  }
  importances <- ranger::importance(rf, type = 1, scale = FALSE)
  if(invert.importance)
    importances <- importances * -1
  selected.vars <- order(importances, decreasing = TRUE)
  ordered.importances <- importances[selected.vars]
  initialImportances <- importances
  initialOrderedImportances <- ordered.importances
  j <- 1
  n.vars[j] <- dim(xdata)[2]
  vars[j] <- paste(colnames(xdata), collapse = " + ")
  OOB.rf[j] <- m.iterated.ob.error
  OOB.sd[j] <- sd.iterated.ob.error
  var.simplify <- TRUE
  while (var.simplify) {
    if (verbose) {
      print("gc inside loop of varSelRF")
      print(gc())
    }
    else {
      gc()
    }
    last.rf <- rf
    last.vars <- selected.vars
    previous.m.error <- m.iterated.ob.error
    previous.sd.error <- sd.iterated.ob.error
    if (length(selected.vars) <= 2) {
      var.simplify <- FALSE
      break
    }
    if (recompute.var.imp & (j > 1)) {
      importances <- ranger::importance(rf, type = 1, scale = FALSE)
      if(invert.importance)
        importances <- importances * -1
      tmp.order <- order(importances, decreasing = TRUE)
      selected.vars <- selected.vars[tmp.order]
      ordered.importances <- importances[tmp.order]
    }
    num.vars <- length(selected.vars)
    if (is.null(vars.drop.num)) 
      vars.drop <- round(num.vars * vars.drop.frac)
    else vars.drop <- vars.drop.num
    if (num.vars >= (vars.drop + 2)) {
      if (vars.drop == 0) {
        vars.drop <- 1
        if ((num.vars - vars.drop) < 1) 
          stop("vars.drop = 0 and num.vars -vars.drop < 1!")
      }
      selected.vars <- selected.vars[1:(num.vars - vars.drop)]
      ordered.importances <- ordered.importances[1:(num.vars - 
                                                      vars.drop)]
    } else {
      selected.vars <- selected.vars[1:2]
      ordered.importances <- ordered.importances[1:2]
    }
    if ((length(selected.vars) < 2) | (any(selected.vars < 
                                           1))) {
      var.simplify <- FALSE
      break
    }
    mtry <- floor(sqrt(length(selected.vars)) * mtryFactor)
    if (mtry > length(selected.vars)) 
      mtry <- length(selected.vars)
    if (recompute.var.imp){
      if(rfpackage == "ranger"){
        rf <- ranger(data = data.df[, c(colnames(xdata)[selected.vars],"class")], dependent.variable.name = "class", 
                     importance = "impurity", num.trees = ntree, mtry = mtry, 
                     keep.inbag = T)
      }else{
        rf <- randomForest(x = xdata, y = Class, ntree = ntree, mtry = mtry, importance = TRUE, keep.forest = keep.forest)
      }
    }else {
      if(rfpackage == "ranger"){
        rf <- ranger(data = data.df[, c(colnames(xdata)[selected.vars],"class")], dependent.variable.name="class", 
                     importance = "none", num.trees = ntreeIterat, mtry = mtry, 
                     keep.inbag = T)
      }else{
        rf <- randomForest(x = xdata, y = Class, ntree = ntree, mtry = mtry, importance = F, keep.forest = keep.forest)
      }
    }
    m.iterated.ob.error <- oobError(rf)
    sd.iterated.ob.error <- sqrt(m.iterated.ob.error * (1 - 
                                                          m.iterated.ob.error) * (1/num.subjects))
    if (verbose) {
      print(paste("..... iteration ", j, "; OOB error: mean = ", 
                  round(m.iterated.ob.error, 4), "; sd = ", round(sd.iterated.ob.error, 
                                                                  4), "; num. vars = ", length(selected.vars), 
                  sep = ""))
    }
    j <- j + 1
    n.vars[j] <- length(selected.vars)
    vars[j] <- paste(colnames(xdata)[selected.vars], collapse = " + ")
    OOB.rf[j] <- m.iterated.ob.error
    OOB.sd[j] <- sd.iterated.ob.error
    if (!whole.range & ((m.iterated.ob.error > (m.initial.ob.error + 
                                                c.sd * sd.initial.ob.error)) | (m.iterated.ob.error > 
                                                                                (previous.m.error + c.sd * previous.sd.error)))) 
      var.simplify <- FALSE
  }
  if (!whole.range) {
    if (!is.null(colnames(xdata))) 
      selected.vars <- sort(colnames(xdata)[last.vars])
    else selected.vars <- last.vars
    out <- list(selec.history = data.frame(Number.Variables = n.vars, 
                                           Vars.in.Forest = vars, OOB = OOB.rf, sd.OOB = OOB.sd)[1:j, 
                                                                                                 ], rf.model = last.rf, selected.vars = selected.vars, 
                selected.model = paste(selected.vars, collapse = " + "), 
                best.model.nvars = length(selected.vars), initialImportances = initialImportances, 
                initialOrderedImportances = initialOrderedImportances, 
                ntree = ntree, ntreeIterat = ntreeIterat, mtryFactor = mtryFactor, 
                firstForest = FirstForest)
    class(out) <- "varSelRF"
    return(out)
  } else {
    n.vars <- n.vars[1:j]
    vars <- vars[1:j]
    OOB.rf <- OOB.rf[1:j]
    OOB.sd <- OOB.sd[1:j]
    min.oob.ci <- min(OOB.rf) + c.sd * OOB.sd[which.min(OOB.rf)]
    best.pos <- which(OOB.rf <= min.oob.ci)[which.min(n.vars[which(OOB.rf <= 
                                                                     min.oob.ci)])]
    selected.vars <- sort(unlist(strsplit(vars[best.pos], 
                                          " + ", fixed = TRUE)))
    out <- list(selec.history = data.frame(Number.Variables = n.vars, 
                                           Vars.in.Forest = vars, OOB = OOB.rf, sd.OOB = OOB.sd), 
                rf.model = NA, selected.vars = selected.vars, selected.model = paste(selected.vars, 
                                                                                     collapse = " + "), best.model.nvars = n.vars[best.pos], 
                initialImportances = initialImportances, initialOrderedImportances = initialOrderedImportances, 
                ntree = ntree, ntreeIterat = ntreeIterat, mtryFactor = mtryFactor, 
                firstForest = FirstForest)
    class(out) <- "varSelRF"
    return(out)
  }
}

# ------------------------------------------------------------------------------

## Plotting for backwards elimination

plot_acc_size = function(res){
  plot.df = bind_rows(lapply(seq_along(res)[-length(res)], function(x,res){
    r = res[[x]]
    return(data.frame("size" = length(r$VarSelRF$selected.vars),
                      "accuracy" = r$test.accuracy))
    
  },res))
  
  p = ggplot(plot.df,aes(x=size,y=accuracy))+
    geom_point(size=1.5,color="red",alpha=0.5)+
    theme_bw()+
    geom_smooth(method=lm)+
    labs(title="Backwards elimination optimal model size by test accuracy",
         x="Model size", y = "Accuracy")
  
  return(p)
}

plot_runs_oob = function(res){
  plot.df = bind_rows(lapply(seq_along(res)[-length(res)], function(x,res){
    r = res[[x]]
    df = r[["VarSelRF"]][["selec.history"]][,c(1,3)]
    #df$run_ind = rep(x,nrow(df))
    return(df)
  },res))
  
  plot.df = plot.df %>%
    group_by(Number.Variables) %>%
    summarize(
      m = mean(OOB, na.rm=TRUE),
      sd = sd(OOB)
    ) 
  
  p = ggplot(plot.df,aes(x=Number.Variables,y=m))+
    geom_line()+
    geom_ribbon(aes(ymin=m-2*sd, ymax=m+2*sd), alpha=0.4)+
    scale_x_reverse()+
    theme_bw()+
    #geom_smooth(method=lm)+
    labs(title="Backwards elimination run model sizes by OOB",
         x="Run model size", y = "OOB")
  return(p)
}

plot_freq_analysis = function(res){
  all.genes = unlist(lapply(seq_along(res)[-length(res)], function(x,res){
    return(res[[x]]$VarSelRF$selected.vars)
  },res))
  
  freq.gene = as.data.frame(table(all.genes))
  colnames(freq.gene) = c("Gene","f")
  
  freq.gene$Gene = factor(freq.gene$Gene,levels = freq.gene$Gene[order(freq.gene$f)])
  
  p = ggplot(freq.gene,aes(x=Gene,y=f))+
    geom_bar(stat = "identity",width = 0.5) +
    theme(axis.text.x = element_text(angle = 90)) + coord_flip()+
    theme_bw()+
    labs(title="Gene frequency amongst optimal models",
         x="Gene", y = "Freq")
  
  return(p)
}


################################################################################

##########
##########                      Main Work Flow 
##########

################################################################################

# ------------------------------------------------------------------------------

##
##  0. Run Backwards Elination
##

# Set cores for parallel
cores = detectCores() - 1

# Run backwards elimination in parallel. Output res is a list of the runs including:
# optimal model performance, optimal model feature set, and the run statastics.

res = doParBackElim(data_path = "./example data/example_bc.RDS",
            name = "example_dat",
            trees = 100,
            dropFrac = 0.1,
            runs = 30,
            cores = 1) # <--Do runs in parallel


# ------------------------------------------------------------------------------

##
##  1. Visualise results
##


## 1.1 Model size within backwards elimination runs by OOB
p1 = plot_runs_oob(res)

## 1.2 test accruacy by opitaml model size 
p2 = plot_acc_size(res)

## 1.3 Gene frequency analysis
p3 = plot_freq_analysis(res)

## plot all together
ggarrange(p3,
          ggarrange(p1, p2, nrow = 2, labels = c("B", "C")),
          nrow = 1, 
          ncol = 2,
          labels = "A"   
) 



doParBackElim <- function(data_path, name, trees = 5000, dropFrac = 0.001, runs = 100,
                        featureSet = NULL,cores = 40){
  
  data <- readRDS(data_path)
  
  train_dat <- data$train
  test_data <- data$test
  val_data <- data$val
  
  
  #Subset if freatureset not NULL
  if(!is.null(featureSet)){
    features <- readRDS(featureSet)
    train_dat <- train_dat[, which(colnames(train_dat) %in% c(features,"sample.type")), with = F]
    test_data <- test_data[, which(colnames(test_data) %in% c(features,"sample.type")), with = F]
    val_data <- val_data[, which(colnames(val_data) %in% c(features,"sample.type")), with = F]
  }
  
  
  col.lim <- (ncol(train_dat) - 1)
  
  # Format names to not interfere with random forest run
  name.mapping <- data.frame("orig" = colnames(train_dat)[1:col.lim], "new" = paste("Gene",1:col.lim,sep = ""))
  colnames(train_dat)[1:col.lim] <- paste("Gene",1:col.lim,sep = "")
  colnames(test_data)[1:col.lim] <- paste("Gene",1:col.lim,sep = "")
  colnames(val_data)[1:col.lim] <- paste("Gene",1:col.lim,sep = "")
  saveRDS(name.mapping,file = paste("./Output/",name,".name.mapping.RDS",sep = ""))
  
  
  # Make list the length of intention
  myls <- vector("list", length = runs)
  names(myls) <- 1:runs
  
  
  print(paste("Doing: ",name,sep = ""))
  
  # Too count progress
  for (i in 1:length(myls)) {
    myls[[i]]$count <- i
  }
  
  save.list <- mclapply(myls, function(x){
    
    # Run varselRF
    backPs <- varSelRF2(xdata = train_dat[,1:col.lim],Class = train_dat$sample.type,
                        ntree = trees,
                        vars.drop.frac = dropFrac,
                        whole.rang = T,
                        verbose = T,
                        recompute.var.imp = T)
    
    # Calc class weights - to adjust random forest objective function.
    data.s <-rbind(data$train,data$test,data$val)
    cls.tab <- table(data.s$sample.type)
    cls.tab <- (1-(cls.tab/sum(cls.tab)))
    cls.weight <- as.numeric(cls.tab)
    names(cls.weight) <- names(cls.tab)
    
    # Run ranger (fast random forest) over optimised feature set to create final model.
    rfd <- ranger(data = subset(train_dat, select = c(backPs$selected.vars,"sample.type")),
                  dependent.variable.name = 'sample.type',
                  num.trees= trees,
                  class.weights = cls.weight)
    
    
    # Make prediction
    pred.test <- predict(rfd, subset(test_data, select = c(backPs$selected.vars,"sample.type")))
    pred.val <- predict(rfd, subset(val_data, select = c(backPs$selected.vars,"sample.type")))
    
    confus.test <- confusionMatrix(pred.test$predictions,test_data$sample.type,mode='everything')
    confus.val <- confusionMatrix(pred.val$predictions,val_data$sample.type,mode='everything')
    
    test.accuracy <- confus.test[["overall"]][1]
    val.accuracy <- confus.val[["overall"]][1]
    
    
    return.list <- list(
      "VarSelRF" = backPs,
      "Ranger" = rfd,
      "conf.test" = confus.test,
      "conf.val" = confus.val,
      "test.accuracy" = test.accuracy,
      "val.accuracy" = val.accuracy
    )
    
    
    return(return.list)
    
  },mc.cores = cores)
  
  saveRDS(save.list, file = paste("./Output/",name,".save.list.RDS",sep = ""))
  save.list$name.mapping = name.mapping
  return(save.list)
  
}


# ------------------------------------------------------------------------------

## Varself RF modified function to work with Ranger

varSelRF2 <- function (xdata, Class, c.sd = 1, mtryFactor = 1, ntree = 5000, 
                       ntreeIterat = 2000, vars.drop.num = NULL, vars.drop.frac = 0.2, 
                       whole.range = TRUE, recompute.var.imp = FALSE, verbose = FALSE, 
                       returnFirstForest = TRUE, fitted.rf = NULL, keep.forest = FALSE, rfpackage="ranger", invert.importance=F) 
{
  if (!is.factor(Class)) 
    stop("Class should be a factor")
  if ((is.null(vars.drop.num) & is.null(vars.drop.frac)) | 
      (!is.null(vars.drop.num) & !is.null(vars.drop.frac))) 
    stop("One (and only one) of vars.drop.frac and vars.drop.num must be NULL and the other set")
  max.num.steps <- dim(xdata)[2]
  num.subjects <- dim(xdata)[1]
  if (is.null(colnames(xdata))) 
    colnames(xdata) <- paste("v", 1:dim(xdata)[2], sep = "")
  n.vars <- vars <- OOB.rf <- OOB.sd <- rep(NA, max.num.steps)
  oobError <- function(rf) {
    if(class(rf) == "randomForest"){
      ooo <- rf$confusion[, -dim(rf$confusion)[2]]
    }else if(class(rf) == "ranger"){
      ooo <- rf$confusion
    }
    s.ooo <- sum(ooo)
    diag(ooo) <- 0
    sum(ooo)/s.ooo
  }
  data.df <- data.frame(xdata,class=Class)
  if (!is.null(fitted.rf)) {
    if (ncol(fitted.rf$importance) < 2) 
      stop("The fitted rf was not fitted with importance = TRUE")
    n.ntree <- fitted.rf$ntree
    mtry <- fitted.rf$mtry
    n.mtryFactor <- mtry/sqrt(ncol(xdata))
    if ((n.ntree != ntree) | (n.mtryFactor != mtryFactor)) 
      warning("Using as ntree and mtry the parameters obtained from fitted.rf", 
              immediate. = TRUE)
    ntree <- n.ntree
    mtryFactor <- n.mtryFactor
    rm(n.ntree, n.mtryFactor)
    rf <- fitted.rf
  } else {
    mtry <- floor(sqrt(ncol(xdata)) * mtryFactor)
    if(rfpackage == "ranger"){
      rf <- ranger(data = data.df, dependent.variable.name="class", num.trees = ntree, 
                   mtry = mtry, importance = "impurity", keep.inbag=T)
    }else{
      rf <- randomForest(x = xdata, y = Class, ntree = ntree, mtry = mtry, importance = TRUE, keep.forest = keep.forest)
    }
  }
  if (returnFirstForest){ 
    FirstForest <- rf
  }else{
    FirstForest <- NULL
  }
  m.iterated.ob.error <- m.initial.ob.error <- oobError(rf)
  sd.iterated.ob.error <- sd.initial.ob.error <- sqrt(m.iterated.ob.error * 
                                                        (1 - m.iterated.ob.error) * (1/num.subjects))
  if (verbose) {
    print(paste("Initial OOB error: mean = ", round(m.initial.ob.error, 
                                                    4), "; sd = ", round(sd.initial.ob.error, 4), sep = ""))
  }
  importances <- ranger::importance(rf, type = 1, scale = FALSE)
  if(invert.importance)
    importances <- importances * -1
  selected.vars <- order(importances, decreasing = TRUE)
  ordered.importances <- importances[selected.vars]
  initialImportances <- importances
  initialOrderedImportances <- ordered.importances
  j <- 1
  n.vars[j] <- dim(xdata)[2]
  vars[j] <- paste(colnames(xdata), collapse = " + ")
  OOB.rf[j] <- m.iterated.ob.error
  OOB.sd[j] <- sd.iterated.ob.error
  var.simplify <- TRUE
  while (var.simplify) {
    if (verbose) {
      print("gc inside loop of varSelRF")
      print(gc())
    }
    else {
      gc()
    }
    last.rf <- rf
    last.vars <- selected.vars
    previous.m.error <- m.iterated.ob.error
    previous.sd.error <- sd.iterated.ob.error
    if (length(selected.vars) <= 2) {
      var.simplify <- FALSE
      break
    }
    if (recompute.var.imp & (j > 1)) {
      importances <- ranger::importance(rf, type = 1, scale = FALSE)
      if(invert.importance)
        importances <- importances * -1
      tmp.order <- order(importances, decreasing = TRUE)
      selected.vars <- selected.vars[tmp.order]
      ordered.importances <- importances[tmp.order]
    }
    num.vars <- length(selected.vars)
    if (is.null(vars.drop.num)) 
      vars.drop <- round(num.vars * vars.drop.frac)
    else vars.drop <- vars.drop.num
    if (num.vars >= (vars.drop + 2)) {
      if (vars.drop == 0) {
        vars.drop <- 1
        if ((num.vars - vars.drop) < 1) 
          stop("vars.drop = 0 and num.vars -vars.drop < 1!")
      }
      selected.vars <- selected.vars[1:(num.vars - vars.drop)]
      ordered.importances <- ordered.importances[1:(num.vars - 
                                                      vars.drop)]
    } else {
      selected.vars <- selected.vars[1:2]
      ordered.importances <- ordered.importances[1:2]
    }
    if ((length(selected.vars) < 2) | (any(selected.vars < 
                                           1))) {
      var.simplify <- FALSE
      break
    }
    mtry <- floor(sqrt(length(selected.vars)) * mtryFactor)
    if (mtry > length(selected.vars)) 
      mtry <- length(selected.vars)
    if (recompute.var.imp){
      if(rfpackage == "ranger"){
        rf <- ranger(data = data.df[, c(colnames(xdata)[selected.vars],"class")], dependent.variable.name = "class", 
                     importance = "impurity", num.trees = ntree, mtry = mtry, 
                     keep.inbag = T)
      }else{
        rf <- randomForest(x = xdata, y = Class, ntree = ntree, mtry = mtry, importance = TRUE, keep.forest = keep.forest)
      }
    }else {
      if(rfpackage == "ranger"){
        rf <- ranger(data = data.df[, c(colnames(xdata)[selected.vars],"class")], dependent.variable.name="class", 
                     importance = "none", num.trees = ntreeIterat, mtry = mtry, 
                     keep.inbag = T)
      }else{
        rf <- randomForest(x = xdata, y = Class, ntree = ntree, mtry = mtry, importance = F, keep.forest = keep.forest)
      }
    }
    m.iterated.ob.error <- oobError(rf)
    sd.iterated.ob.error <- sqrt(m.iterated.ob.error * (1 - 
                                                          m.iterated.ob.error) * (1/num.subjects))
    if (verbose) {
      print(paste("..... iteration ", j, "; OOB error: mean = ", 
                  round(m.iterated.ob.error, 4), "; sd = ", round(sd.iterated.ob.error, 
                                                                  4), "; num. vars = ", length(selected.vars), 
                  sep = ""))
    }
    j <- j + 1
    n.vars[j] <- length(selected.vars)
    vars[j] <- paste(colnames(xdata)[selected.vars], collapse = " + ")
    OOB.rf[j] <- m.iterated.ob.error
    OOB.sd[j] <- sd.iterated.ob.error
    if (!whole.range & ((m.iterated.ob.error > (m.initial.ob.error + 
                                                c.sd * sd.initial.ob.error)) | (m.iterated.ob.error > 
                                                                                (previous.m.error + c.sd * previous.sd.error)))) 
      var.simplify <- FALSE
  }
  if (!whole.range) {
    if (!is.null(colnames(xdata))) 
      selected.vars <- sort(colnames(xdata)[last.vars])
    else selected.vars <- last.vars
    out <- list(selec.history = data.frame(Number.Variables = n.vars, 
                                           Vars.in.Forest = vars, OOB = OOB.rf, sd.OOB = OOB.sd)[1:j, 
                                                                                                 ], rf.model = last.rf, selected.vars = selected.vars, 
                selected.model = paste(selected.vars, collapse = " + "), 
                best.model.nvars = length(selected.vars), initialImportances = initialImportances, 
                initialOrderedImportances = initialOrderedImportances, 
                ntree = ntree, ntreeIterat = ntreeIterat, mtryFactor = mtryFactor, 
                firstForest = FirstForest)
    class(out) <- "varSelRF"
    return(out)
  } else {
    n.vars <- n.vars[1:j]
    vars <- vars[1:j]
    OOB.rf <- OOB.rf[1:j]
    OOB.sd <- OOB.sd[1:j]
    min.oob.ci <- min(OOB.rf) + c.sd * OOB.sd[which.min(OOB.rf)]
    best.pos <- which(OOB.rf <= min.oob.ci)[which.min(n.vars[which(OOB.rf <= 
                                                                     min.oob.ci)])]
    selected.vars <- sort(unlist(strsplit(vars[best.pos], 
                                          " + ", fixed = TRUE)))
    out <- list(selec.history = data.frame(Number.Variables = n.vars, 
                                           Vars.in.Forest = vars, OOB = OOB.rf, sd.OOB = OOB.sd), 
                rf.model = NA, selected.vars = selected.vars, selected.model = paste(selected.vars, 
                                                                                     collapse = " + "), best.model.nvars = n.vars[best.pos], 
                initialImportances = initialImportances, initialOrderedImportances = initialOrderedImportances, 
                ntree = ntree, ntreeIterat = ntreeIterat, mtryFactor = mtryFactor, 
                firstForest = FirstForest)
    class(out) <- "varSelRF"
    return(out)
  }
}

# ------------------------------------------------------------------------------

## Function for random forest class (Galgo)

rF.class <- function (chr, parent, tr, te, result) 
{
  require(randomForest)
  avg <- parent$data$avg
  train <- parent$data$data[tr, as.numeric(chr),drop=F]
  test <- parent$data$data[te, as.numeric(chr),drop=F]
  train.cl <- parent$classes[tr]
  test.cl <- parent$classes[te]
  trfd <- randomForest(train,train.cl,test,test.cl,ntree=parent$data$ntree,keep.forest=F,type="prob",mtry=parent$data$mtry,classwt=parent$data$cls.weight)
  #res <<- list(chr=chr,tr=tr,te=te,res=result,train=train,test=test,train.cl=train.cl,test.cl=test.cl)
  
  if(result == 1){	
    if(avg){
      classes <- levels(parent$classes)
      vals <- sapply(sapply(classes,function(x) which(test.cl == x)),function(y) sum(trfd$test$predicted[y] == test.cl[y])/length(y))
      xx <- mean(vals)
    }else{
      xx <- sum(trfd$test$predicted == test.cl)/length(test.cl)
    }
    if (is.numeric(xx)) 
      return(xx)
    else 0.01
  }else{
    return(trfd$test$predicted)
  }
}



################################################################################

##########
##########                      Main Work Flow 
##########

################################################################################

# ------------------------------------------------------------------------------

##
##  1. Backwards Elination
##

# Set cores for parallel
cores = detectCores() - 1

# Run backwards elimination in parallel. Output res is a list of run outputs, 
# optimal model performance, gene name mapping to model.
res = doParBackElim(data_path = "./example data/example_bc.RDS",
            name = "example_dat",
            trees = 100,
            dropFrac = 0.1,
            runs = 1,
            cores = 1) # <--Do runs in parallel


