library(dplyr)
library(ggplot2)
library(pROC)
library(patchwork)
library(ggpubr)
library(ggtext)
library(here)

here::i_am("4-markdown-figures/tcga_investigation.R")

dbPath <- fs::path_expand("~/Dropbox")
pal.jama <- ggsci::pal_jama()(5)
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

BuccalHypo_CervHyper <- dmps[dmps$dir_type.x<0 & dmps$dir_type.y>0,]$cg # BuccalHypo, CervHyper
BuccalHyper_CervHyper <- dmps[dmps$dir_type.x>0 & dmps$dir_type.y>0,]$cg # BuccalHyper, CervHyper
BuccalHypo_CervHypo <- dmps[dmps$dir_type.x<0 & dmps$dir_type.y<0,]$cg # BuccalHypo, CervHypo
BuccalHyper_CervHypo <- dmps[dmps$dir_type.x>0 & dmps$dir_type.y<0,]$cg # BuccalHyper, CervHypo


# TCGA Data
# prep pheno
load("~/Dropbox/eca/sola/1-download/2-output/data.Rdata")
data <- data |> 
  dplyr::mutate(type = factor(type, levels = c("Control", 'Cancer')))

plotTCGA <- function(beta, 
                     pheno,
                     cgs = c('BuccalHypo_CervHyper',
                             'BuccalHyper_CervHyper',
                             'BuccalHypo_CervHypo',
                             'BuccalHyper_CervHypo'),
                     project = 'BRCA'){
  
  # First, subset pheno
  intersect <- intersect(colnames(beta), pheno$barcode1)
  
  tmp1 <- beta[,match(intersect, colnames(beta))]
  tmp2 <- pheno[match(intersect, pheno$barcode1),]

  if(!identical(colnames(tmp1), tmp2$barcode1)){
    stop('names not identical')
  }
  
  # Then compute CpGs
  for(c in cgs){
    ind <- rownames(tmp1) %in% get(c)
    tmp2[[c]] <- colMeans(tmp1[ind,])
  }
  
  # Then pivot longer and plot
  tmp3 <- tmp2 |> 
    tidyr::pivot_longer(any_of(cgs),
                        names_to = 'group',
                        values_to = 'value')
  
  plot_box <- tmp3 |> 
    ggplot(aes(x = type,
               y = value)) +
    geom_boxplot(aes(fill = type),
                 alpha = 0.3,
                 outlier.shape = NA) +
    # ggbeeswarm::geom_beeswarm(aes(colour = type),
    #                           size = 0.2) +
    facet_wrap(~group,
               nrow = 1) +
    ggpubr::stat_compare_means(comparisons = list(c('Control', 'Cancer'))) +
    theme_bw() +
    scale_colour_manual(values = c(colour.Control, colour.Breast),
                        aesthetics = c('fill', 'colour')) +
    labs(x = '',
         y = 'mean methylation',
         subtitle = project) +
    theme(legend.position = 'none')

  # plot the AUCs
  a1 <- plot_single_roc(tmp2$type, tmp2$BuccalHyper_CervHyper)
  a2 <- plot_single_roc(tmp2$type, tmp2$BuccalHyper_CervHypo)
  a3 <- plot_single_roc(tmp2$type, tmp2$BuccalHypo_CervHyper)
  a4 <- plot_single_roc(tmp2$type, tmp2$BuccalHypo_CervHypo)

  
  plot <- plot_box / (a1|a2|a3|a4)
  
  print(plot)
  }



data <- data |> 
  dplyr::group_by(project) |> 
  dplyr::filter(all(c('Control', 'Cancer') %in% type)) |> 
  dplyr::ungroup()


projects <- unique(data$project)


pdf("tcga_investigation.pdf", width = 11, height = 6.5)

for (p in projects){
  
  cat('processing ', p)
  # Load 
  load(paste0("~/Dropbox/data/tcga/", p, "/beta.Rdata"))
  
  plotTCGA(beta, pheno = data,
           project = p)
  
}

dev.off()
