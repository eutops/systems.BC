glmnet_gaussian <- function(beta_tr, y_tr,
                            beta_val, y_val,
                            alpha = 0.0,
                            lambda = "(-4,4)"){
  
  
  ## Edit 30/10/2020:
  # Adding in a customisable lambda
  ####  Default lambda = exp(-4,4) (back-compatibility)
  ####  "unrestricted" = default
  #### other: define own sequence by (-x, y)
  
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
                 y = y_tr)
  
  dat_val <- list(x = beta_val,
                  y = y_val)
  
  # perform k-fold cross validation
  NFOLD <- 10
  foldid <- rep(1:NFOLD,ceiling(ncol(dat_tr$x)/NFOLD))[1:ncol(dat_tr$x)]
  
  fit.cv <- cv.glmnet(x = t(dat_tr$x), 
                      y = as.numeric(dat_tr$y),
                      lambda = lambda,
                      type.measure="mse",
                      alpha=alpha,
                      family="gaussian",
                      nfold=NFOLD,
                      foldid=foldid)
  
  
  # compute training auc
  tr_predictor <- predict(fit.cv,
                          newx=t(dat_tr$x),
                          s="lambda.min")
  
  
  tr_mse <- mean((tr_predictor-dat_tr$y)^2)
  
  # compute validation auc
  val_predictor <- predict(fit.cv,
                           newx=t(dat_val$x),
                           s="lambda.min")
  
  val_mse <- mean((val_predictor-dat_val$y)^2)
  
  ### Warning message if lambda.min too close to minimum/maximum
  if (!is.null(lambda)) {
    if (fit.cv$lambda.min < (min(lambda) + min(lambda) * 0.1) |
        fit.cv$lambda.min > (max(lambda) - max(lambda) * 0.1)) {
      warning("Optimum lambda is close to limits of lambda (range may be too small)")
    }
  }
  
  # return results
  return(list(tr_mse=tr_mse,
              val_mse=val_mse,
              fit.cv=fit.cv,
              oob_mse=min(fit.cv$cvm),
              tr_predictor=tr_predictor,
              val_predictor=val_predictor))
  
}