plot_double_roc <- function(type1, index1, type2, index2, col1 = "black", col2 = "blue", title1 = "all", title2 = "below 30", style = "default",
                            direction1 = "<",direction2 = "<", generalTitle="",textSize=2.9){
  
  require(pROC)
  require(ggsci)
  cols <- pal_lancet(palette="lanonc", alpha = 0.8)(8) 
  
  #--------------------------------#  
  # make annotation labels
  auc <- round(as.numeric(roc(type1,index1, quiet=T, ci=T, direction = direction1)$auc),digits=2)
  cil <- round(as.numeric(roc(type1,index1, quiet=T, ci=T, direction = direction1)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type1,index1, quiet=T, ci=T, direction = direction1)$ci[3]),digits=2)
  anno1 <- paste('AUC (', title1, ') = ',auc,'\n(95% CI: ', cil,'-',ciu,')',sep='')
  
  # make annotation labels
  auc2 <- round(as.numeric(roc(type2,index2, quiet=T, ci=T, direction = direction2)$auc),digits=2)
  cil2 <- round(as.numeric(roc(type2,index2, quiet=T, ci=T, direction = direction2)$ci[1]),digits=2)
  ciu2 <- round(as.numeric(roc(type2,index2, quiet=T, ci=T, direction = direction2)$ci[3]),digits=2)
  anno2 <- paste('AUC (', title2, ') = ',auc2,'\n(95% CI: ', cil2,'-',ciu2,')',sep='')
  
  if(style == "lancet"){
    anno1 <- gsub("[.]", "·", anno1)
    anno2 <- gsub("[.]", "·", anno2)
  }
  
  #--------------------------------#  
  
  roc1 <- roc(type1, index1, direction = direction1)
  roc2 <- roc(type2, index2, direction = direction2)
  title1 <- as.character(title1)
  title2 <- as.character(title2)

  
  ggplot() +
    geom_path(aes(x=1-roc1$specificities,
                  y=(roc1$sensitivities)),
              colour = col1,
              size = 0.7) +
    geom_path(aes(x=1-roc2$specificities,
                  y=(roc2$sensitivities)),
              colour = col2,
              size = 0.7) +
    annotate("segment", x = 0, y = 0,
             xend = 1, yend = 1,
             colour = "gray60") +
    xlab("1 - Specificity") +
    ylab("Sensitivity") +
    theme_minimal() +
    theme(plot.title = element_text(size=10),
          panel.grid = element_blank(),
          axis.line = element_line()) +
    annotate(geom='label',
             x=0.55,
             y=0.35,
             label=anno1,
             fill=col1,
             colour='white',
             size=textSize)  +
    annotate(geom='label',
             x=0.55,
             y=0.1,
             label=anno2,
             fill=col2,
             colour='white',
             size=textSize)  +
    ggtitle(generalTitle)  
  
}
