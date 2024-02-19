el_classifier <- function(beta_tr, type_tr,
                          beta_val, type_val,
                          alpha = 0.0,
                          lambda = "(-4,4)",
                          slope = "default",
                          nfolds = 10){
  
  require(Hmisc)
  require(rms)
  require(pROC)
  require(glmnet)
  
  lambda <- if (lambda == "(-4,4)"){
    exp(seq(-4,4,by=0.1))
  } else {
    if (lambda == "unrestricted"){
      NULL
    } else {
      x <- strsplit(lambda, split = "[(|,|)+]")
      lambda.min <- x[[1]][2]
      lambda.max <- x[[1]][3]
      exp(seq(lambda.min,lambda.max,by=0.1))
    }
  }
  
  
  # create training and validation datasets
  dat_tr <- list(x = beta_tr,
                 y = type_tr)
  
  dat_val <- list(x = beta_val,
                  y = type_val)
  
  # perform k-fold cross validation
  NFOLD <- nfolds
  foldid <- rep(1:NFOLD,ceiling(ncol(dat_tr$x)/NFOLD))[1:ncol(dat_tr$x)]
  
  fit.cv <- cv.glmnet(x = t(dat_tr$x), 
                      y = as.factor(as.vector(dat_tr$y)),
                      lambda = lambda,
                      type.measure="auc",
                      alpha=alpha,
                      family="binomial",
                      nfold=NFOLD,
                      foldid=foldid)
  
  
  # compute training auc
  predictor <- predict(fit.cv,
                       newx=t(dat_tr$x),
                       s="lambda.min")
  
  tr_predictor <- predict(fit.cv,
                          newx=t(dat_tr$x),
                          s="lambda.min",
                          type='response')
  
  tr_roc <- roc(dat_tr$y, as.numeric(predictor),quiet = TRUE)
  
  
  # compute validation auc
  predictor <- predict(fit.cv,
                       newx=t(dat_val$x),
                       s="lambda.min")
  
  val_predictor <- predict(fit.cv,
                           newx=t(dat_val$x),
                           s="lambda.min",
                           type='response')
  
  val_roc <- roc(dat_val$y, as.numeric(predictor), quiet = TRUE)
  
  # calibration slope
  x <- if(slope == "default"){
    ifelse(type_val=="Control",1,0)
  } else {
    ifelse(type_val=="Control",0,1)
  }
  r <- val.prob(val_predictor[,1], x, pl=FALSE)
  
  # # -------- ROC curve ------- #
  # y <- rep(0,length(dat_val$y))
  # y[dat_val$y=='Endometrial'] <- 1
  # 
  # x <- as.numeric(predictor)
  # 
  # y <- y[order(x, decreasing = TRUE)]
  # TPR=cumsum(y)/sum(y)
  # FPR=cumsum(!y)/sum(!y)
  # 
  # plot(c(0,FPR),c(0,TPR),type='l',
  #      xlim=c(0,1),ylim=c(0,1),
  #      xlab='False positives',
  #      ylab='True postivies',
  #      main='PC1 ROC curve')
  # points(c(-1,2),c(-1,2),lty='dashed',type='l')
  # grid()
  # 
  
  # print results
  #cat('Training AUC = ',tr_roc$auc,'\n','Validation AUC = ',val_roc$auc,'\n',sep='')
  
  # Warning message if lambda.min too close to minimum/maximum
  if (fit.cv$lambda.min < (min(lambda)+min(lambda)*0.1) | fit.cv$lambda.min > (max(lambda) - max(lambda)*0.1)) warning("Optimum lambda is close to limits of lambda (range may be too small)")
  
  # return results
  return(list(tr_roc=tr_roc,
              val_roc=val_roc,
              fit.cv=fit.cv,
              oob_auc=mean(fit.cv$cvm),
              tr_predictor=tr_predictor,
              val_predictor=val_predictor,
              r=r))
  
}
