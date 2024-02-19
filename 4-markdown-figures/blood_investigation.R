# blood investigation, '* Extremes of Quadrants in Fig 6B in blood'

library(dplyr)
library(ggplot2)
library(pROC)
library(patchwork)
library(ggpubr)
library(ggtext)

here::i_am("4-markdown-figures/blood_investigation.R")

dbPath <- fs::path_expand("~/Dropbox")

# Two cases: Breast vs Control
colour.Breast <- pal.jama[4]
colour.Control <- pal.jama[3]

source(here("0-data","src", "plot_single_roc.R"))


# Load Quadrants
# Load DMPs
buccal <- readRDS(here("1-ewas", "1-buccal", "delta-beta.Rds")) |> 
  as.data.frame() |> 
  dplyr::mutate(padj = p.adjust(type, method = 'BH')) |> 
  tibble::rownames_to_column('cg')

cervical <- readRDS(here("1-ewas", "2-cervical", "delta-beta.Rds")) |> 
  as.data.frame() |> 
  dplyr::mutate(padj = p.adjust(type, method = 'BH'))|> 
  tibble::rownames_to_column('cg') 

# Merge
dmps <- buccal |> 
  dplyr::full_join(cervical, by = 'cg') |> 
  dplyr::filter(type.y < 0.05 & type.x < 0.05) # Keep those significant in both


# Two quadrants of interest: shared; extremes; # Take any with shared db > 0.02
hypo <- dmps[dmps$db_immune.x< -0.03 & dmps$db_immune.y< -0.03,]$cg
hyper <-dmps[dmps$db_immune.x>0.03 & dmps$db_immune.y>0.03,]$cg


# plot in blood data
plotBlood <- function(beta_path,
                      pheno_path,
                      title = ''){
  
  assign('beta', get(load(beta_path)))
  assign('pheno', get(load(pheno_path)))
  
  intersect <- intersect(pheno$basename, colnames(beta))
  pheno <- pheno[match(intersect, pheno$basename),]
  beta <- beta[,intersect]
  
  pheno <- pheno |> 
    dplyr::mutate(type = factor(type, levels = c("Control", 'BC case')))

  pheno$hyper <- colMeans(beta[na.omit(match(hyper, rownames(beta))),])
  pheno$hypo <- colMeans(beta[na.omit(match(hypo, rownames(beta))),])
  
  pheno2 <- pheno |> 
    tidyr::pivot_longer(c(hypo, hyper),
                        names_to = 'db_immune_cgs',
                        values_to = 'value')
  
  plot_box <- pheno2 |> 
    ggplot(aes(x = type,
               y = value)) +
    geom_boxplot(aes(fill = type),
                 alpha = 0.3,
                 outlier.shape = NA) +
    # ggbeeswarm::geom_beeswarm(aes(colour = type),
    #                           size = 0.2) +
    facet_wrap(~db_immune_cgs,
               nrow = 1,
               scales = 'free') +
    ggpubr::stat_compare_means(comparisons = list(c('Control', 'BC case'))) +
    theme_bw() +
    scale_colour_manual(values = c(colour.Control, colour.Breast),
                        aesthetics = c('fill', 'colour')) +
    labs(x = '',
         y = 'mean methylation',
         subtitle = title) +
    theme(legend.position = 'none')
  
  a1 <- plot_single_roc(pheno$type, pheno$hyper)
  a2 <- plot_single_roc(pheno$type, pheno$hypo)
  
  plot <-   plot_box / (a1|a2)
  
  return(plot)
   
  print(plot)
  
  }


beta_path = "~/Dropbox/data/3c-blood/beta_merged.Rdata"
pheno_path = here("0-data","dataframes-plotting","3c_blood_validation.Rdata") # 3C
pdf("blood_3c.pdf", width = 11, height = 6.5)
plotBlood(beta_path, pheno_path, title = '3C')
dev.off()


beta_path = "~/Dropbox/data/blood/GSE237036_wang/beta_merged.Rdata"
pheno_path = here("0-data","dataframes-plotting","blood_validation.Rdata") # blood
pdf("blood_wang.pdf", width = 11, height = 6.5)
plotBlood(beta_path, pheno_path, title = 'wang')
dev.off()
