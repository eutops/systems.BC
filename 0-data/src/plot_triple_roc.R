plot_triple_roc <- function(type1, index1, type2, index2, type3, index3, col1 = "black", col2 = "blue", col3 = "blue", title1 = "all", title2 = "below 30", title3 = "below 30", style = "default",
                            direction1 = "<",direction2 = "<", direction3 = "<", generalTitle="",textSize=2.9, textCol1="white",textCol3="white",x1=0.55,y1=0.55,x2=0.55,y2=0.35,x3=0.55,y3=0.1){
  
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
  
  # make annotation labels
  auc3 <- round(as.numeric(roc(type3,index3, quiet=T, ci=T, direction = direction2)$auc),digits=2)
  cil3 <- round(as.numeric(roc(type3,index3, quiet=T, ci=T, direction = direction2)$ci[1]),digits=2)
  ciu3 <- round(as.numeric(roc(type3,index3, quiet=T, ci=T, direction = direction2)$ci[3]),digits=2)
  anno3 <- paste('AUC (', title3, ') = ',auc3,'\n(95% CI: ', cil3,'-',ciu3,')',sep='')
  
  if(style == "lancet"){
    anno1 <- gsub("[.]", "??", anno1)
    anno2 <- gsub("[.]", "??", anno2)
    anno3 <- gsub("[.]", "??", anno3)
  }
  
  #--------------------------------#  
  
  roc1 <- roc(type1, index1, direction = direction1)
  roc2 <- roc(type2, index2, direction = direction2)
  roc3 <- roc(type3, index3, direction = direction3)
  title1 <- as.character(title1)
  title2 <- as.character(title2)
  title3 <- as.character(title3)
  
  ggplot() +
    geom_path(aes(x=1-roc1$specificities,
                  y=(roc1$sensitivities)),
              colour = col1,
              size = 0.7) +
    geom_path(aes(x=1-roc2$specificities,
                  y=(roc2$sensitivities)),
              colour = col2,
              size = 0.7) +
    geom_path(aes(x=1-roc3$specificities,
                  y=(roc3$sensitivities)),
              colour = col3,
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
             x=x1,
             y=y1,
             label=anno1,
             fill=col1,
             colour=textCol1,
             size=textSize)  +
    annotate(geom='label',
             x=x2,
             y=y2,
             label=anno2,
             fill=col2,
             colour='white',
             size=textSize)  +
      annotate(geom='label',
               x=x3,
               y=y3,
               label=anno3,
               fill=col3,
               colour=textCol3,
               size=textSize) +
    ggtitle(generalTitle)  
  
}
