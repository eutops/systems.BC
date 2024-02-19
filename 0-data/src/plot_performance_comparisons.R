combine_sets <- function(resid, x, type = "default"){
  if (type == "default"){
    m <- matrix(nrow = length(resid$n_seq), ncol = 5)
    colnames(m) <- c("n", "AUC_val", "AUC_oob", "AUC_tr", "Type")
    
    m[,1] <- resid$n_seq
    m[,2] <- resid$val_auc
    m[,3] <- resid$oob_auc
    m[,4] <- resid$tr_auc
    m[,5] <- rep(x, length(resid$n_seq))
    
    return(m)
  } else {
    m <- matrix(nrow = length(resid$n_seq), ncol = 7)
    colnames(m) <- c("n", "AUC_val", "AUC_oob", "AUC_tr", "Type", "Slope", "intercept")
    
    m[,1] <- resid$n_seq
    m[,2] <- resid$val_auc
    m[,3] <- resid$oob_auc
    m[,4] <- resid$tr_auc
    m[,5] <- rep(x, length(resid$n_seq))
    m[,6] <- resid$cal_slope
    m[,7] <- resid$cal_intercept
    return (m)
  }
}

plot_performance_compare <- function(data, col, sets = 6, legendPos="top"){
  colour.pal.d8 <- c("#EA7580","#F2949A","#F6B3A1","#D8C99E","#3AB9AC","#109DB7","#0C6EA5","#172869")
  colour.pal.c <- colorRampPalette(colour.pal.d8)
  
  data |> 
    ggplot() +
    geom_line(aes(x = n,
                  y = col,
                  colour = Type)) +
    ylim(c(0.2, 1)) +
    xlab("Number of input CpGs") +
    ylab("AUC") +
    theme_minimal() +
    scale_colour_manual(values = colour.pal.c(sets),
                        aesthetics = "colour")+
    theme(legend.position = legendPos,
          legend.title = element_blank())
}


plot_performance_compare_slope <- function(data, col, sets = 6){
  colour.pal.d8 <- c("#EA7580","#F2949A","#F6B3A1","#D8C99E","#3AB9AC","#109DB7","#0C6EA5","#172869")
  colour.pal.c <- colorRampPalette(colour.pal.d8)
  
  data %>%
    ggplot() +
    geom_line(aes(x = n,
                  y = col,
                  colour = Type)) +
    ylim(c(-1, 1)) +
    xlab("Number of input CpGs") +
    ylab("Slope") +
    theme_minimal() +
    scale_colour_manual(values = colour.pal.c(sets),
                        aesthetics = "colour")+
    theme(legend.position = "top",
          legend.title = element_blank())
}

