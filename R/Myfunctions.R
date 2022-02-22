library(tidyverse)
library(glmnet)
library(dplyr)
library(randomForest)
library(e1071)
#load("ukb_phenotypes_rolledup.rda")

# Z: covariate of interest; A: observed response; X: Confounding features
# FDR: desirable FDR level.
# d.interaction: 'T' for using dICRT; F for using 'd0CRT'
# k: dimension of the control variables in distilled data, only useful when "d.interaction = T"

# model = c('Binomial_lasso', 'RF'), for linear model, logistics model (binary Y) and random forest (dICRT (RF)).
#
# M = 50000: Number of resampling if needed
# MC_free: whether to use resampling-free d0CRT or dICRT

dCRT<-function(A, Z, X, mean_Z, 
               Gen_Z = NULL, 
               model = 'Binomial_lasso', k = NULL, M = 5000, RF.num.trees = c(100, 30),
               MC_free = F){
  
  ##### logistic/RF #####
  if (model == 'Binomial_lasso'){
    model <- 'binomial'
  }
  ##### d0CRT/dICRT #####
  if (MC_free == T){
    type = 'cov'
  }else{
    type = 'beta'
  }
  
  p<-ncol(X)
  n<-length(A)
  
  ##### Distill Z #####
  Z_res_ob <- Z - mean_Z

  
  ##### Distill A #####
  if(model== 'binomial'){
    ##### dICRT #####
    if(k>=1){
      cv_lasso <- cv.glmnet(X, A, alpha = 1, family = model, 
                            lambda = NULL, dfmax = as.integer(p / 2))
      lamb <- cv_lasso$lambda.min
      opt_model <- glmnet(X, A, alpha = 1, lambda = lamb, 
                          family = model, dfmax = as.integer(p / 2))
      
      beta_leave <- opt_model$beta
      beta_sort <- sort(abs(as.vector(beta_leave)), decreasing = T, index.return = T)
      index_use <- beta_sort$ix[1:k]
      X_use <- X[,index_use]
      offsets <- as.vector(opt_model$a0 + X %*% opt_model$beta)
      eps_res <- A - 1 / (1 + exp(- offsets))
      
      ##### resample-free/not #####
      if(type== 'cov'){}
      
      if(type== 'beta'){
        weight_inter <- 1 / sqrt(k)
        W_inter <- cbind(1, X_use) * as.vector(eps_res)
        W_inter_X <- cbind(1, X_use) * as.vector(Z_res_ob)
        XTX <- t(W_inter_X) %*% W_inter_X
        XTX_inv <- solve(XTX)
        Z_dI <- diag(c(1, rep(weight_inter, k))) %*% XTX_inv %*% t(W_inter) %*% Z_res_ob
        #Z_dI <- diag(c(1, rep(weight_inter, k))) %*% XTX_inv %*% t(W_inter_X) %*% eps_res
        imp_obe <- sum(Z_dI^2)
      }
    } else{
    ##### d0CRT #####
    cv_lasso_null <- cv.glmnet(X, A, alpha = 1, family = model, 
                               lambda = NULL, dfmax = as.integer(p / 2))
    lamb_null <- cv_lasso_null$lambda.min
    model_res_null <- glmnet(X, A, alpha = 1, lambda = lamb_null,
                             family = model, dfmax = as.integer(p / 2))
    eps_res <- A - 1 / (1 + exp(- predict(model_res_null, X)))
    
    ##### resample free/not #####
    if(type== 'cov'){}
    if(type== 'beta'){
      imp_obe <- abs(mean(Z_res_ob * eps_res)) / mean(Z_res_ob^2)
      imp_obe_raw <- mean(Z_res_ob * eps_res) / mean(Z_res_ob^2)
    }
    }
  } 
  
  if(model== 'RF'){
    n_tree <- RF.num.trees[1]
    n_tree_imp <- RF.num.trees[2]
    ##### dICRT #####
    if (k >= 1){
      rf.fit <- randomForest(x = X, y = as.factor(A), ntree = n_tree)
      offsets <- rf.fit$predicted
      imp_fit <- as.vector(rf.fit$importance)
      imp_sort <- sort(imp_fit, decreasing = T, index.return = T)
      index_use <- imp_sort$ix[1:k]
      X_use <- X[,index_use]
      
      # Observed importance
      
      rf.indx.fit <- randomForest(x = as.matrix(cbind(X - mean_Z, offsets, X_use)), 
                                  as.factor(A), ntree = n_tree_imp)
      imp_obe <- rf.indx.fit$importance[1]
      imp_obe_raw <- imp_obe
      
    }else{
      rf.fit <- randomForest(x = X - mean_Z, y = as.factor(A), ntree = n_tree)
      eps_res <- A - as.numeric(rf.fit$predicted)
      imp_obe <- abs(mean(Z_res_ob * eps_res))
      imp_obe_raw <- mean(Z_res_ob * eps_res)
      
      # Observed importance
      if (MC_free == F | is.null(Gen_Z)){

      }
    }
  }
  ############### Estimation p-values ##################
  
  ##### d0CRT #####
  if (k == 0){
    if (MC_free == T){}else{
      Z_resample <- Gen_Z(num = M,mean_Z=mean_Z)
      Z_res_resample <- Z_resample - mean_Z
      
      if (type == 'cov'){
        t_lst <- t(Z_res_resample) %*% eps_res / n
        t_lst_raw <- c(imp_obe_raw, t_lst)
        t_lst <- c(imp_obe, abs(t_lst))
      }
      if (type == 'beta'){
        t_lst <- unlist(lapply(c(1:M), function(j){
          mean(Z_res_resample[,j] * eps_res) / mean((Z_res_resample[,j])^2)
        })) 
        t_lst_raw <- c(imp_obe_raw, t_lst)
        t_lst <- c(imp_obe, abs(t_lst))
      }
      pvl <- mean(ifelse(t_lst >= imp_obe, 1, 0))
    }
    
    
  }else{
    ##### dICRT #####
    if (model == 'RF'){
      #### RF #####
      Z_vec <- c(1)
      t_lst <- c()
      for (m in 1:M) {
        if (is.null(Gen_Z)){
          delta_gen <- rnorm(n, 0, sqrt(sigma2_x))
          Z_sample <- mean_X + delta_gen
        }else{
          Z_sample <- Gen_Z(mean_Z=mean_Z, num = 1)
        }
        rf.sam.fit <- randomForest(x = cbind(as.numeric(Z_sample) - mean_Z, offsets, X_use), 
                                   as.factor(A), ntree = n_tree_imp)
        imp_resample <- rf.sam.fit$importance[1]
        Z <- ifelse(imp_resample >= imp_obe, 1, 0)
        Z_vec <- c(Z_vec, Z)
        t_lst <- c(t_lst, imp_resample)
        t_lst_raw <- t_lst
      }
      pvl <- mean(Z_vec)
    }
    
    if (model == 'binomial'){
      #### Logistic #####
      if (type == 'cov'){}
      
      if (type == 'beta'){
        Z_resample <- Gen_Z(mean_Z=mean_Z, num = M)
        Z_res_sample <- Z_resample - mean_Z
      }
      
      W_inter <- cbind(1, X_use) * as.vector(eps_res)
      weight_inter <- 1 / sqrt(k)
      
      t_lst <- unlist(lapply(c(1:M), function(j){
        #W_inter_X <- cbind(1, X_use) * as.vector(Z_res_sample[,j])
        #XTX <- t(W_inter_X) %*% W_inter_X
        #XTX_inv <- solve(XTX)
        Z_dI <- diag(c(1, rep(weight_inter, k))) %*% XTX_inv %*% t(W_inter) %*% Z_res_sample[,j]
        #Z_dI <- diag(c(1, rep(weight_inter, k))) %*% XTX_inv %*% t(W_inter_X) %*% eps_res
        sum(Z_dI^2)
      })) 
      t_lst_raw <- t_lst
      pvl <- mean(c(1, ifelse(t_lst >= imp_obe, 1, 0)))  
    }
    
    
  }
  
  
  return(list(pvl = pvl, imp_obe = imp_obe, t_lst_raw = t_lst_raw))
}

####### data generation functions ######
#### X: data set we use to generate A and Z ####
#### num: number of obs we use in X; n obs will be randomly chosen from X,default is 502510 ####
#### d: number of Z variates we want to generate from X ####
#### model: 'logistic' , 'nonlinear' ,the model to generate data ####



fit_cond_Z<-function(X, Z, model='logistic', n_tree=100){
  n<-length(Z)
  p<-ncol(X)
  if(model=='logistic'){
    model='binomial'
    cv_lasso <- cv.glmnet(X, Z, alpha = 1, family = model, 
                          dfmax = as.integer(p / 2))
    lamb <- cv_lasso$lambda.min
    opt_model <- glmnet(X, Z, alpha = 1, lambda = lamb, 
                        family = model, dfmax = as.integer(p / 2))
    #### fitted Z ####
    mean_Z<-1 / (1 + exp(- predict(opt_model, X)))
  }
  if(model=='RF'){
    rf.fit <- randomForest(x = X, y = as.factor(Z), ntree = n_tree)
    #### fitted Z ####
    mean_Z <- as.numeric(predict(rf.fit,type = "prob")[,2])
  }
  return(mean_Z=as.vector(mean_Z))
}

example_Gen_Z<-function(mean_Z, num){
  n<-length(mean_Z)
  resample_Z<-t(sapply(mean_Z,FUN = function(c){rbinom(num,1,prob = c)}))
  return(resample_Z)
}


dCRT_multiple<-function(A, Z, X, 
                        model = 'Binomial_lasso', k = NULL, M = 5000, RF.num.trees = c(100, 30),
                        MC_free = F, FDR=0.1){
  d<-ncol(Z)
  pvl<- sapply(c(1:d), FUN = function(j){
    dCRT(A=A, X=X, Z=Z[,j], mean_Z=fit_cond_Z(X=X,Z=Z[,j]),
         model = model, k = k, M = M, MC_free = MC_free, Gen_Z = example_Gen_Z, RF.num.trees = RF.num.trees)}) %>%
    p.adjust(method = 'BH')
  selection_set_CRT <- which(pvl <= FDR)
  return(list(pvl, selection_set_CRT))
}


  




