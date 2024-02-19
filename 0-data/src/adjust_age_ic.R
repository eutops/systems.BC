# Adjust values by age and ic
# Indices can be freely defined

adjust_age_ic <- function(data,
                          indices = c("index_general_age",
                                      "horvath",
                                      "hannum",
                                      "index_epithelial_age",
                                      "index_immune_age",
                                      "PhenoAge",
                                      "DNAmTL",
                                      "pcgtage"),
                          tissue,
                          tissue_adj = FALSE){
  
  # initiate columns
  for (i in indices){
    data[[paste0(i, "_adj")]] <- numeric(nrow(data))
  }
  
  if(tissue_adj == F){
    for (i in indices){
      tmp <- data[data$type %in% c("Control", "Normal", "normal"),]
      fit <- lm(tmp[[i]] ~ age + ic, data = tmp)
      data[,paste0(i, "_adj")] <- data[[i]] - as.numeric(predict(fit, newdata = data))
    }
  } else {
    for (t in unique(tissue)){
      for (i in indices){
        tmp <- data[data$type=="Control" & tissue == t,]
        fit <- lm(tmp[[i]] ~ age + ic, data = tmp)
        data[tissue == t,paste0(i, "_adj")] <- data[tissue == t,][[i]] - as.numeric(predict(fit, newdata = data[tissue == t,]))
      }
    }
  }
  
  return(data)
}