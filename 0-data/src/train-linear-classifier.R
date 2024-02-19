# Description: train linear elastic net classifers as a function of the number of
# input variables.
# 
# Author: James E. Barrett
# Contact: regmjeb@ucl.ac.uk
# Date: 21 October 2019
# Edit: 10 November 2020: customisable lambda fuction
# Edit: 24 March 2020: adding in calibration intercept/slope for evaluation

train_linear_classifier <- function(beta_tr,
                                    type_tr,
                                    beta_val,
                                    type_val,
                                    alpha,
                                    lambda = "(-4,4)",
                                    n_seq,
                                    slope = "default",
                                    nfolds=10){
   
   if(max(n_seq) > nrow(beta_tr)){
      stop('n larger than beta number of rows')
   }
   
   val_auc <- rep(NA, length(n_seq))  # validation AUC
   tr_auc <- rep(NA, length(n_seq))   # training AUC
   oob_auc <- rep(NA, length(n_seq))  # out-of-bag validation AUC estimates
   lambda <- lambda # Customisable lambda
   cal_intercept <- rep(NA, length(n_seq)) # calibration intercept
   cal_slope <- rep(NA, length(n_seq)) # calibration slope
   val_predictor <- rep(NA, length(n_seq))
   
   counter <- 1
   for (n in n_seq){
      cat('\nRunning n =', n,'\n')
      res <- el_classifier(beta_tr[1:n,], type_tr,
                           beta_val[1:n,], type_val,
                           alpha=alpha,
                           lambda = lambda,
                           slope = slope,
                           nfolds=nfolds)
      val_auc[counter] <- as.numeric(res$val_roc$auc)
      tr_auc[counter] <- as.numeric(res$tr_roc$auc)
      oob_auc[counter] <- as.numeric(res$oob_auc)
      cal_intercept[counter] <- as.numeric(res$r["Intercept"])
      cal_slope[counter] <- as.numeric(res$r["Slope"])
      # val_predictor[counter] <- as.numeric(res$val_predictor)
      counter <- counter + 1
      cat('\nAUC =', round(res$val_roc$auc, digits=3),'\nIntercept = ', round(res$r["Intercept"], 3), '\nSlope = ', round(res$r["Slope"],3), '\n')
   }
   
   return(list(n_seq=n_seq,
               tr_auc=tr_auc,
               val_auc=val_auc,
               val_predictor=val_predictor,
               oob_auc=oob_auc,
               cal_intercept=cal_intercept,
               cal_slope=cal_slope))
}



