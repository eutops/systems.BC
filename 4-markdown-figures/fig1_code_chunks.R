# Fig 1 plots

library(patchwork)
library(ggtext)
library(ggplot2)
library(here)
library(ggpubr)
library(gtsummary)
library(dplyr)
library(gt)
library(tidyr)
library(epitools)

dbPath <- fs::path_expand("~/Dropbox")
dbPath <- sub(" Dropbox.*", " Dropbox/Bente Theeuwes", getwd())

here::i_am("4-markdown-figures/fig1_code_chunks.R")

pal.jama <- ggsci::pal_jama()(5)


# B-middle: OR plot -----------------------------

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


# Plot: OVERALL direction -----------------------------------

# BuccalHypo, CervHyper
p = nrow(dmps[dmps$dir_type.x<0 & dmps$dir_type.y>0,])/nrow(dmps)
x = nrow(dmps[dmps$dir_type.x<0 & dmps$dir_type.y>0,])

BuccalHypo_CervHyper <- ORquadrant(p, x, n = nrow(dmps))
anno_BHoCHer <- paste0('Buccal <b style="color:#0000FF">hypo</b><br>Cervical <b style="color:#FF0000">hyper</b><br>OR = ', format(round(BuccalHypo_CervHyper$or, 2), nsmall = 2), '<br>[', 
                       format(round(BuccalHypo_CervHyper$ci_lo, 2), nsmall = 2),
                       ' - ',
                       format(round(BuccalHypo_CervHyper$ci_hi, 2), nsmall = 2), ']')

# BuccalHyper, CervHyper
p = nrow(dmps[dmps$dir_type.x>0 & dmps$dir_type.y>0,])/nrow(dmps)
x = nrow(dmps[dmps$dir_type.x>0 & dmps$dir_type.y>0,])

BuccalHyper_CervHyper <- ORquadrant(p, x, n = nrow(dmps))
BuccalHyper_CervHyper
anno_BHerCHer <- paste0('Buccal <b style="color:#FF0000">hyper</b><br>Cervical <b style="color:#FF0000">hyper</b><br>OR = ', format(round(BuccalHyper_CervHyper$or, 2), nsmall = 2), '<br>[', 
                        format(round(BuccalHyper_CervHyper$ci_lo, 2), nsmall = 2),
                        ' - ',
                        format(round(BuccalHyper_CervHyper$ci_hi, 2), nsmall = 2), ']')

# BuccalHypo, CervHypo
p = nrow(dmps[dmps$dir_type.x<0 & dmps$dir_type.y<0,])/nrow(dmps)
x = nrow(dmps[dmps$dir_type.x<0 & dmps$dir_type.y<0,])

BuccalHypo_CervHypo <- ORquadrant(p, x, n = nrow(dmps))
anno_BHoCHo <- paste0('Buccal <b style="color:#0000FF">hypo</b><br>Cervical <b style="color:#0000FF">hypo</b><br>OR = ', format(round(BuccalHypo_CervHypo$or, 2), nsmall = 2), '<br>[', 
                      format(round(BuccalHypo_CervHypo$ci_lo, 2), nsmall = 2),
                      ' - ',
                      format(round(BuccalHypo_CervHypo$ci_hi, 2), nsmall = 2), ']')

# BuccalHyper, CervHypo
p = nrow(dmps[dmps$dir_type.x>0 & dmps$dir_type.y<0,])/nrow(dmps)
x = nrow(dmps[dmps$dir_type.x>0 & dmps$dir_type.y<0,])

BuccalHyper_CervHypo <- ORquadrant(p, x, n = nrow(dmps))
anno_BHerCHo <- paste0('Buccal <b style="color:#FF0000">hyper</b><br>Cervical <b style="color:#0000FF">hypo</b><br>OR = ', format(round(BuccalHyper_CervHypo$or, 2), nsmall = 2), '<br>[', 
                       format(round(BuccalHyper_CervHypo$ci_lo, 2), nsmall = 2),
                       ' - ',
                       format(round(BuccalHyper_CervHypo$ci_hi, 2), nsmall = 2), ']')


overall <- dmps |> 
  dplyr::filter(type.y < 0.05 & type.x < 0.05) |> 
  ggplot(aes(x = dir_type.x,
             y = dir_type.y)) +
  geom_hex(bins = 100) +
  labs(x = 'Buccal Δβ',
       y = 'Cervical Δβ') +
  scale_fill_viridis_c(name = "Number of CpGs", option = "plasma") +
  # ggpubr::stat_cor() +
  theme_minimal() + # Add this line to remove the grey background
  theme(legend.position = "right", legend.box = "vertical") + # Add these lines to place legend above the plot 
  annotate('richtext',
           x = 0.05, y = 0.07,
           label = anno_BHerCHer,
           label.size = 0.2,
           alpha = 0.8,
           hjust = 0
  )  +
  annotate('richtext',
           x = -0.05, y = 0.07,
           label = anno_BHoCHer,
           label.size = 0.2,
           alpha = 0.8,
           hjust = 1) +
  annotate('richtext',
           x = -0.05, y = -0.04,
           label = anno_BHoCHo,
           label.size = 0.2,
           alpha = 0.8,
           hjust = 1,
           vjust = 1) +
  annotate('richtext',
           x = 0.05, y = -0.04,
           label = anno_BHerCHo,
           label.size = 0.2,
           alpha = 0.8,
           hjust = 0,
           vjust = 1) +
  geom_vline(xintercept = 0,
             alpha = 0.6) +
  geom_hline(yintercept = 0,
             alpha = 0.6)







# eFORGE -----------------------------

plot.eFORGE <- function(path, max = 5){
  
  results<-read.table(path, header = TRUE, sep="\t") |> 
    dplyr::select(Pvalue, Cell, Tissue, Qvalue, Zscore) |> 
    dplyr::filter(Qvalue < 0.05) |> 
    dplyr::arrange(Pvalue)
  
  if(nrow(results) < max){
    results <- results |> 
      dplyr::slice(1:nrow(results))
  } else {
    results <- results |>
      dplyr::slice(1:max)
  }
  
  if(nrow(results) != 0){
  plot <- results |> 
    dplyr::mutate(label = stringr::str_wrap(Cell, width = 20)) |> 
    ggplot(aes(y = forcats::fct_reorder(label, -Pvalue),
               x = -log10(Pvalue))) +
    geom_point(size = 3,
               alpha = 0.7,
               colour = pal.jama[5]) +
    theme_bw() +
    theme(aspect.ratio = 1.75) +
    labs(x = '-log10(pvalue)',
         y = '') +
    xlim(c(0, 30))
  
  return(list(results = results, plot = plot))
  
  } else {
    cat("no significant results")
  }
}

## eFORGE buccal hypo, cerv hyper --------------------
path <- file.path(here("0-data","dataframes-plotting","eforge_hypo_hyper.tsv.gz"))

BHo_CHer <- plot.eFORGE(path)
BHo_CHer$plot

## eFORGE buccal hyper, cerv hyper --------------------
path <- file.path(here("0-data","dataframes-plotting","eforge_hyper_hyper.tsv.gz"))
BHer_CHer <- plot.eFORGE(path)
BHer_CHer$results

## eFORGE buccal hyper, cerv hypo --------------------
path <- file.path(here("0-data","dataframes-plotting","eforge_hyper_hypo.tsv.gz"))
BHer_CHo <- plot.eFORGE(path)
BHer_CHo$results
BHer_CHo$plot

## eFORGE buccal hypo, cerv hypo --------------------
path <- file.path(here("0-data","dataframes-plotting","eforge_hypo_hypo.tsv.gz"))
BHo_CHo <- plot.eFORGE(path)
BHo_CHo$results
BHo_CHo$plot


# Breast plots ------------------------------

load(here("0-data/GSE22845/beta_merged.Rdata")) # Load beta
load(here("0-data/GSE22845/pheno_trimmed.Rdata")) # Load Pheno
identical(colnames(beta_merged), pheno_trimmed$basename)

pheno <- pheno_trimmed |> 
  dplyr::mutate(type = ifelse(type == 'adj_norm', 'normal-adjacent', type),
                type = factor(type, levels = c('normal', 'normal-adjacent', 'tumor')))
beta <- beta_merged[,pheno$basename]
identical(colnames(beta), pheno$basename)

# Load quadrant list of CpGs
BuccalHypo_CervHyper <- dmps[dmps$dir_type.x<0 & dmps$dir_type.y>0,]$cg # BuccalHypo, CervHyper
BuccalHyper_CervHyper <- dmps[dmps$dir_type.x>0 & dmps$dir_type.y>0,]$cg # BuccalHyper, CervHyper
BuccalHypo_CervHypo <- dmps[dmps$dir_type.x<0 & dmps$dir_type.y<0,]$cg # BuccalHypo, CervHypo
BuccalHyper_CervHypo <- dmps[dmps$dir_type.x>0 & dmps$dir_type.y<0,]$cg # BuccalHyper, CervHypo

plot.Breast <- function(cgs, beta, pheno){
  
  ind <- rownames(beta) %in% cgs
  pheno$mean <- colMeans(beta[ind,])
  
  plot <- pheno |> 
    ggplot(aes(x = type,
               y = mean)) +
    geom_boxplot(aes(fill = type),
                 alpha = 0.2,
                 outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(size = 0.3,
                                 aes(colour = type),
                                 alpha = 0.5) +
    theme_bw() +
    theme(aspect.ratio = 1.5,
          legend.position = 'top',
          legend.box = "vertical",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +  # Remove x-axis ticks and labels
    labs(x = '', y = 'mean beta in breast tissue') +
    ggpubr::stat_compare_means(comparisons = list(c('normal', 'tumor'),
                                                  c('normal-adjacent', 'tumor'),
                                                  c('normal', 'normal-adjacent')),
                               label = "p.format") +
    scale_colour_manual(values = pal.jama[c(3,2,4)],
                        aesthetics = c("fill", 'colour')) 
  
  return(plot)
}


Breast_BHo_CHer <- plot.Breast(BuccalHypo_CervHyper, beta, pheno) # BuccalHypo, CervHyper

Breast_BHer_CHer <- plot.Breast(BuccalHyper_CervHyper, beta, pheno) # BuccalHyper, CervHyper

Breast_BHo_CHo <- plot.Breast(BuccalHypo_CervHypo, beta, pheno) # BuccalHypo, CervHypo

Breast_BHer_CHo <- plot.Breast(BuccalHyper_CervHypo, beta, pheno) # BuccalHyper, CervHypo



# Compile and save------------------------------

# 
# 
# layout <- "
# ABCD
# EEEE
# EEEE
# FGHI
# "
# 
# 
# pdf_width <- 8.3  # Width in inches (a4)
# pdf_height <- 11.7   # Height in inches (a4: 11.7)
# 
# 
# 
# 
# compiled <- BHo_CHer$plot + Breast_BHo_CHer + guide_area() + Breast_BHer_CHer  + overall + BHer_CHo$plot + Breast_BHo_CHo + BHo_CHo$plot + Breast_BHer_CHo + plot_layout(guides = 'collect')
# 
# compiled <-  compiled + plot_annotation(tag_levels = 'a') + 
#   plot_layout(design = layout)


