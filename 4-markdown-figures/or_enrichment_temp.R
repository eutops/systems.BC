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

# dbPath <- fs::path_expand("~/Dropbox")
dbPath <- sub(" Dropbox.*", " Dropbox/Bente Theeuwes", getwd())

here::i_am("4-markdown-figures/or_enrichment_temp.R")

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

# Odds ratios
ORquadrant <- function(p, x, n){
  
  # following functions here:
  # https://pages.stat.wisc.edu/~st571-1/06-tables-4.pdf
  p0 = 0.25 # Null hypothesis is equal distribution - 0.25
  or <- ((p/(1-p))/(p0/((1-p0)))) # compute OR
  
  x2 <- round(n/4)
  se <- sqrt ( (1/x) + (1/(n-x)) + (1/x2) + (1/(n-x2)) )
  
  
  df <- data.frame(or = or,
                   ci_lo = or - 1.96*se,
                   ci_hi = or + 1.96*se)
  return(df)
  
}

# Plot: DELTA BETA EPITHELIAL direction -----------------------------------
# BuccalHypo, CervHyper
p = nrow(dmps[dmps$db_epithelial.x<0 & dmps$db_epithelial.y>0,])/nrow(dmps)
x = nrow(dmps[dmps$db_epithelial.x<0 & dmps$db_epithelial.y>0,])

BuccalHypo_CervHyper <- ORquadrant(p, x, n = nrow(dmps))
anno_BHoCHer <- paste0('Buccal <b style="color:#0000FF">hypo</b><br>Cervical <b style="color:#FF0000">hyper</b><br>OR = ', format(round(BuccalHypo_CervHyper$or, 2), nsmall = 2), '<br>[', 
                        format(round(BuccalHypo_CervHyper$ci_lo, 2), nsmall = 2),
                        ' - ',
                        format(round(BuccalHypo_CervHyper$ci_hi, 2), nsmall = 2), ']')

# BuccalHyper, CervHyper
p = nrow(dmps[dmps$db_epithelial.x>0 & dmps$db_epithelial.y>0,])/nrow(dmps)
x = nrow(dmps[dmps$db_epithelial.x>0 & dmps$db_epithelial.y>0,])

BuccalHyper_CervHyper <- ORquadrant(p, x, n = nrow(dmps))
BuccalHyper_CervHyper
anno_BHerCHer <- paste0('Buccal <b style="color:#FF0000">hyper</b><br>Cervical <b style="color:#FF0000">hyper</b><br>OR = ', format(round(BuccalHyper_CervHyper$or, 2), nsmall = 2), '<br>[', 
                        format(round(BuccalHyper_CervHyper$ci_lo, 2), nsmall = 2),
                        ' - ',
                        format(round(BuccalHyper_CervHyper$ci_hi, 2), nsmall = 2), ']')

# BuccalHypo, CervHypo
p = nrow(dmps[dmps$db_epithelial.x<0 & dmps$db_epithelial.y<0,])/nrow(dmps)
x = nrow(dmps[dmps$db_epithelial.x<0 & dmps$db_epithelial.y<0,])

BuccalHypo_CervHypo <- ORquadrant(p, x, n = nrow(dmps))
anno_BHoCHo <- paste0('Buccal <b style="color:#0000FF">hypo</b><br>Cervical <b style="color:#0000FF">hypo</b><br>OR = ', format(round(BuccalHypo_CervHypo$or, 2), nsmall = 2), '<br>[', 
                       format(round(BuccalHypo_CervHypo$ci_lo, 2), nsmall = 2),
                       ' - ',
                       format(round(BuccalHypo_CervHypo$ci_hi, 2), nsmall = 2), ']')

# BuccalHyper, CervHypo
p = nrow(dmps[dmps$db_epithelial.x>0 & dmps$db_epithelial.y<0,])/nrow(dmps)
x = nrow(dmps[dmps$db_epithelial.x>0 & dmps$db_epithelial.y<0,])

BuccalHyper_CervHypo <- ORquadrant(p, x, n = nrow(dmps))
anno_BHerCHo <- paste0('Buccal <b style="color:#FF0000">hyper</b><br>Cervical <b style="color:#0000FF">hypo</b><br>OR = ', format(round(BuccalHyper_CervHypo$or, 2), nsmall = 2), '<br>[', 
                      format(round(BuccalHyper_CervHypo$ci_lo, 2), nsmall = 2),
                      ' - ',
                      format(round(BuccalHyper_CervHypo$ci_hi, 2), nsmall = 2), ']')

epi <- dmps |> 
  dplyr::filter(type.y < 0.05 & type.x < 0.05) |> 
  ggplot(aes(x = db_epithelial.x,
             y = db_epithelial.y)) +
  geom_hex(bins = 100) +
  labs(x = 'Buccal epithelial Δβ',
       y = 'Cervical epithelial Δβ') +
  scale_fill_viridis_c(name = "Number of CpGs", option = "plasma") +
  # ggpubr::stat_cor() +
  theme_minimal() + # Add this line to remove the grey background
  theme(legend.position = "bottom", legend.box = "horizontal") + # Move legend below and set it to horizontal
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
           x = -0.05, y = -0.07,
           label = anno_BHoCHo,
           label.size = 0.2,
           alpha = 0.8,
           hjust = 1,
           vjust = 1) +
  annotate('richtext',
           x = 0.05, y = -0.07,
           label = anno_BHerCHo,
           label.size = 0.2,
           alpha = 0.8,
           hjust = 0,
           vjust = 1) +
  geom_vline(xintercept = 0,
             alpha = 0.6) +
  geom_hline(yintercept = 0,
             alpha = 0.6)


# Plot: DELTA BETA IMMUNE direction -----------------------------------
# BuccalHypo, CervHyper
p = nrow(dmps[dmps$db_immune.x<0 & dmps$db_immune.y>0,])/nrow(dmps)
x = nrow(dmps[dmps$db_immune.x<0 & dmps$db_immune.y>0,])

BuccalHypo_CervHyper <- ORquadrant(p, x, n = nrow(dmps))
anno_BHoCHer <- paste0('Buccal <b style="color:#0000FF">hypo</b><br>Cervical <b style="color:#FF0000">hyper</b><br>OR = ', format(round(BuccalHypo_CervHyper$or, 2), nsmall = 2), '<br>[', 
                       format(round(BuccalHypo_CervHyper$ci_lo, 2), nsmall = 2),
                       ' - ',
                       format(round(BuccalHypo_CervHyper$ci_hi, 2), nsmall = 2), ']')

# BuccalHyper, CervHyper
p = nrow(dmps[dmps$db_immune.x>0 & dmps$db_immune.y>0,])/nrow(dmps)
x = nrow(dmps[dmps$db_immune.x>0 & dmps$db_immune.y>0,])

BuccalHyper_CervHyper <- ORquadrant(p, x, n = nrow(dmps))
BuccalHyper_CervHyper
anno_BHerCHer <- paste0('Buccal <b style="color:#FF0000">hyper</b><br>Cervical <b style="color:#FF0000">hyper</b><br>OR = ', format(round(BuccalHyper_CervHyper$or, 2), nsmall = 2), '<br>[', 
                        format(round(BuccalHyper_CervHyper$ci_lo, 2), nsmall = 2),
                        ' - ',
                        format(round(BuccalHyper_CervHyper$ci_hi, 2), nsmall = 2), ']')

# BuccalHypo, CervHypo
p = nrow(dmps[dmps$db_immune.x<0 & dmps$db_immune.y<0,])/nrow(dmps)
x = nrow(dmps[dmps$db_immune.x<0 & dmps$db_immune.y<0,])

BuccalHypo_CervHypo <- ORquadrant(p, x, n = nrow(dmps))
anno_BHoCHo <- paste0('Buccal <b style="color:#0000FF">hypo</b><br>Cervical <b style="color:#0000FF">hypo</b><br>OR = ', format(round(BuccalHypo_CervHypo$or, 2), nsmall = 2), '<br>[', 
                      format(round(BuccalHypo_CervHypo$ci_lo, 2), nsmall = 2),
                      ' - ',
                      format(round(BuccalHypo_CervHypo$ci_hi, 2), nsmall = 2), ']')

# BuccalHyper, CervHypo
p = nrow(dmps[dmps$db_immune.x>0 & dmps$db_immune.y<0,])/nrow(dmps)
x = nrow(dmps[dmps$db_immune.x>0 & dmps$db_immune.y<0,])

BuccalHyper_CervHypo <- ORquadrant(p, x, n = nrow(dmps))
anno_BHerCHo <- paste0('Buccal <b style="color:#FF0000">hyper</b><br>Cervical <b style="color:#0000FF">hypo</b><br>OR = ', format(round(BuccalHyper_CervHypo$or, 2), nsmall = 2), '<br>[', 
                       format(round(BuccalHyper_CervHypo$ci_lo, 2), nsmall = 2),
                       ' - ',
                       format(round(BuccalHyper_CervHypo$ci_hi, 2), nsmall = 2), ']')


immune <- dmps |> 
  dplyr::filter(type.y < 0.05 & type.x < 0.05) |> 
  ggplot(aes(x = db_immune.x,
             y = db_immune.y)) +
  geom_hex(bins = 100) +
  labs(x = 'Buccal immune Δβ',
       y = 'Cervical immune Δβ') +
  scale_fill_viridis_c(name = "Number of CpGs", option = "plasma") +
  # ggpubr::stat_cor() +
  theme_minimal() + # Add this line to remove the grey background
  theme(legend.position = "bottom", legend.box = "horizontal") + # Move legend below and set it to horizontal
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
           x = -0.05, y = -0.07,
           label = anno_BHoCHo,
           label.size = 0.2,
           alpha = 0.8,
           hjust = 1,
           vjust = 1) +
  annotate('richtext',
           x = 0.05, y = -0.07,
           label = anno_BHerCHo,
           label.size = 0.2,
           alpha = 0.8,
           hjust = 0,
           vjust = 1) +
  geom_vline(xintercept = 0,
             alpha = 0.6) +
  geom_hline(yintercept = 0,
             alpha = 0.6)


compiled <- epi + immune

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


dmps |> 
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
