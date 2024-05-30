adjust_neutro <- function(data,
                          indices = c("index_general_age",
                                      "horvath",
                                      "hannum",
                                      "index_epithelial_age",
                                      "index_immune_age",
                                      "PhenoAge",
                                      "DNAmTL",
                                      "pcgtage"),
                          tissue,
                          tissue_adj = FALSE, appendix="_adj"){
  
  # initiate columns
  for (i in indices){
    data[[paste0(i, appendix)]] <- numeric(nrow(data))
  }
  
  if(tissue_adj == F){
    for (i in indices){
      tmp <- data[data$type=="Control",]
      
      if("age" %in% colnames(data)){
      
      fit <- lm(tmp[[i]] ~ age + hepidish_Neutro, data = tmp)
      
      } else {
        fit <- lm(tmp[[i]] ~ hepidish_Neutro, data = tmp)
        
      }
      
      data[,paste0(i, appendix)] <- data[[i]] - as.numeric(predict(fit, newdata = data))
    }
  } else {
    for (t in unique(tissue)){
      for (i in indices){
        tmp <- data[data$type=="Control" & tissue == t,]
        fit <- lm(tmp[[i]] ~ hepidish_Neutro, data = tmp)
        data[tissue == t,paste0(i, appendix)] <- data[tissue == t,][[i]] - as.numeric(predict(fit, newdata = data[tissue == t,]))
      }
    }
  }
  
  return(data)
}