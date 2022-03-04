#Split 1:input into K sets
#randomly split the n samples into K folds
#The function is use for cross-fitting
#' @import caret
Split<-function(K,input){
  I1 = rep(0,input)
  Newids = sample(1:input,input,replace = F)
  m = input/K
  for(i in 1:K){
    I1[Newids[(((i-1)*m):(i*m-1))+1]] = i
  }
  return(I1)
}

L<-function(R,C,givenEstimator,Ctest){
  #R is a vector; C is a data.frame
  #R can be 0/1 variable or continuous variable
  #Fit R ~ C by a machine learning model specified by givenEstimator
  #Output: Rtest by inputing Ctest to the model trained by (R,C)
  #If R is (0,1), then the output should be probability

  Data = cbind(R,C)
  colnames(Data)[1] = 'R'

  #Using 5-fold to tune parameter for the machine learning model
  tt = length(table(R))
  #If R is 0/1 variable, transform it into a factor variable
  if(tt==2){
    Data$R<-gsub('0','N',Data$R)
    Data$R<-gsub('1','P',Data$R)
    Data$R = as.factor(Data$R)
  }

  fitControl <- trainControl(## 3-fold CV
    method = "cv",
    number = 2,
    verboseIter=FALSE)

  n <- nrow(Data)
  p <- ncol(Data)

  #Fit the model
  if (givenEstimator == 'glmnet'){
    glmnetGrid <- expand.grid(alpha = c(0.5, 1), lambda = 0.25 * exp(0.1 * c(1:30)) * sqrt(log(p) / n))

    Fit <- train(R ~ ., data = Data,
                 method = givenEstimator,
                 trControl = fitControl,
                 verbose = FALSE, tuneGrid = glmnetGrid)

  }

  if (givenEstimator == 'gbm'){
    gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 5),
                            n.trees = c(200),
                            shrinkage = 0.02,
                            n.minobsinnode = c(5, 10))
    Fit <- train(R ~ ., data = Data,
                 method = givenEstimator,
                 trControl = fitControl,
                 verbose = FALSE, tuneGrid = gbmGrid)
  }

  if(tt==2){
    Rtest = predict(Fit,Ctest,type='prob')[,2]
  }else{
    Rtest = predict(Fit,Ctest)
  }

  return(Rtest)
}

Lr0<-function(Y,A,X,K,givenEstimator,Xtest){
  n = length(Y)
  I2 = Split(K,n)
  Mp = rep(0,n)
  ap = rep(0,n)

  aBar = rep(0,nrow(Xtest))

  for(j in 1:K){
    idNj = (I2!=j)
    idj = (I2==j)
    #Fit hat_M^{-k,-j} = L(Y,(A,X); {-k,-j})
    Mp[idj] = L(Y[idNj],cbind(A,X)[idNj,],givenEstimator,cbind(A,X)[idj,])
    #Fit hat_a^{-k,-j} = L(A,X; {-k,-j})
    atemp = L(A[idNj],X[idNj,],givenEstimator,rbind(X[idj,],Xtest))
    ap[idj] = atemp[1:nrow(X[idj,])]
    aBar = aBar + atemp[-(1:nrow(X[idj,]))]
  }
  #Equation (3.8): aBar:= hat_a^{-k}(X^{k}) = 1/K sum_{j=1}^K hat_a^{-k,-j}(X^{k})
  aBar = aBar/K
  if(sum(is.infinite(Mp))>0){
    print('Infinite in Mp')
  }

  Mp[Mp<=1e-8] = 1e-8
  Mp[Mp>=1-1e-8] = 1-1e-8

  Ares2 = A - ap

  #Wi = logit(hat_M^{-k,-j}(A_i,X_i))
  Wp = logit(Mp)
  #Equation (3.7): solve beta^{-k}
  betaNk = sum(Ares2*Wp)/sum(Ares2^2)

  #t^{-k} = L(W,X;{-k}); tNk = t^{-k}(X^{k})
  tNk = L(Wp,X,givenEstimator,Xtest)

  #Defined in Equation (3.8)
  rNk = tNk - betaNk*aBar
  return(rNk)
}

#Fit the model logit(Pr(Y=1|A,X)) = beta0*A + r_0(X)
#Return r^{-k}(X^k), m^{-k}(X^k), k =1,...,K
#K: the number of folder, usually is set as 5
#givenEstimator: the machine learning to be used, now have four choices {'gbm', 'svmLinear2', 'rf', 'nnet'}
DML<-function(Y,A,X,K,givenEstimator){
  n = length(Y)

  mXp = rep(0,n); rXp = rep(0,n)
  I1 = Split(K,n)
  for(k in 1:K){
    idNk0 = (I1!=k)&(Y==0)
    idNk = (I1!=k)
    idk = (I1==k)

    #Fit hat_m^{-k} = L(A,X;{-k} cap {Y==0})
    #Then obtain hat_m^{-k}(X^{-k})
    mXp[idk] = L(A[idNk0],X[idNk0,],givenEstimator,X[idk,])

    #Estimate hat_r^{-k} by Y^{-k}, A^{-k}, X^{-k}
    #Then obtain hat_r^{-k}(X^{k})
    rXp[idk] = Lr0(Y[idNk],A[idNk],X[idNk,],K,givenEstimator,X[idk,])
  }
  return(list('mXp'=mXp,'rXp'=rXp))
}

#Solve the estimation function (3.6) in the paper
Estimate<-function(Y,A,dml){
  resA = A-dml$mXp
  C = sum(resA*(1-Y)*exp(dml$rXp))

  g<-function(beta){
    sum(Y*exp(-beta*A)*resA) - C
  }
  lo = -10; up = 10
  beta0 = uniroot(g,c(lo,up))$root
  beta0
}

#Use the (Gaussian Multiplier) Boostrap method to obtain Confidence Interval
# Boostrap<-function(Y,A,dml,B=1000){
#   resA = A-dml$mXp
#   Betas = rep(0,B)
#   for(b in 1:B){
#     e = rnorm(length(Y),mean=1,sd=1)
#     C = sum(e*resA*(1-Y)*exp(dml$rXp))
#
#     g<-function(beta){
#       sum(e*Y*exp(-beta*A)*resA) - C
#     }
#     lo = -10; up = 10
#     beta0 = tryCatch(uniroot(g,c(lo,up))$root,error = function(e){0})
#     Betas[b] = beta0
#   }
#   c(quantile(Betas[Betas!=0],c(0.025,0.975)),mean(Betas[Betas!=0]),sd(Betas[Betas!=0]))
# }

#Use the (Gaussian Multiplier) Boostrap method to obtain Confidence Interval
Boostrap<-function(Y,A,dml,B=1000){
  beta_obs = Estimate(Y,A,dml)
  resA = A-dml$mXp
  Betas = rep(0,B)
  for(b in 1:B){
    e = rnorm(length(Y),mean=1,sd=1)
    C = sum(e*resA*(1-Y)*exp(dml$rXp))

    g<-function(beta){
      sum(e*Y*exp(-beta*A)*resA) - C
    }
    lo = -10; up = 10
    beta0 = tryCatch(uniroot(g,c(lo,up))$root,error = function(e){0})
    #beta0 = uniroot(g,c(lo,up))$root
    Betas[b] = beta0
  }
  c(quantile(Betas[Betas!=0],c(0.025,0.975)),
    `mean`=mean(Betas[Betas!=0]),
    `sd`=sd(Betas[Betas!=0]),
    `pval`=mean(abs(beta_obs) < abs(Betas[Betas!=0])))
}

logit<-function(x){
  x = log(x/(1-x))
  if(sum(is.na(x))>0){
    print('Na in logit!')
  }
  x
}
