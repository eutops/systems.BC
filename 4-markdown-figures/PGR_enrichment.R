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

# Load PGR CpGs
load("~/Dropbox/index-dev/3C/1-BC/18-TFBS-analysis/16-output/pgr_cpgs_ext3.Rdata")
head(pgr_cpgs)
nrow(pgr_cpgs)


# Buccal Hyper, Cerv Hypo
x = dmps[dmps$db_epithelial.x>0 & dmps$db_epithelial.y<0,]
head(x)
plot(x$db_epithelial.x, x$db_epithelial.y)

# Enrichment for PGRs?
                # not in PGR    # in PGR
# not in quad       A             B
# in Quad           C             D

# Remove PGR CpGs not in the arrray (buccal)
pgr_cpgs <- pgr_cpgs[rownames(pgr_cpgs) %in% buccal$cg,]

# Compute
A <- unique(buccal$cg[!buccal$cg %in% c(x$cg, rownames(pgr_cpgs))])
B <- unique(buccal$cg[buccal$cg %in% rownames(pgr_cpgs) & !buccal$cg %in% x$cg])
C <- unique(x$cg[!x$cg %in% rownames(pgr_cpgs)])
D <- unique(x$cg[x$cg %in% rownames(pgr_cpgs)])

tab <- matrix(c(length(A), length(C), length(B), length(D)), ncol = 2)
epitools::oddsratio(tab)

# Manual:
odds_quad_not_pgr <- length(C)/(sum(length(C), length(A)))
odds_quad_pgr <- length(D)/((sum(length(B), length(D))))
odds_quad_pgr/odds_quad_not_pgr

# Compared to only significant CpGs: -----
# Compute
A <- unique(dmps$cg[!dmps$cg %in% c(x$cg, rownames(pgr_cpgs))])
B <- unique(dmps$cg[dmps$cg %in% rownames(pgr_cpgs) & !dmps$cg %in% x$cg])
C <- unique(x$cg[!x$cg %in% rownames(pgr_cpgs)])
D <- unique(x$cg[x$cg %in% rownames(pgr_cpgs)])

tab <- matrix(c(length(A), length(C), length(B), length(D)), ncol = 2)
tab
epitools::oddsratio(tab)

# Manual:
odds_quad_not_pgr <- length(C)/(sum(length(C), length(A)))
odds_quad_pgr <- length(D)/((sum(length(B), length(D))))
odds_quad_pgr/odds_quad_not_pgr # underenrichment


# Odds of being overall significant if PGR -------------
A <- unique(buccal$cg[!buccal$cg %in% c(dmps$cg, rownames(pgr_cpgs))])
B <- unique(buccal$cg[buccal$cg %in% rownames(pgr_cpgs) & !buccal$cg %in% dmps$cg])
C <- unique(dmps$cg[!dmps$cg %in% rownames(pgr_cpgs)])
D <- unique(dmps$cg[dmps$cg %in% rownames(pgr_cpgs)])

tab <- matrix(c(length(A), length(C), length(B), length(D)), ncol = 2)
tab
epitools::oddsratio(tab)

# Manual:
odds_quad_not_pgr <- length(C)/(sum(length(C), length(A)))
odds_quad_pgr <- length(D)/((sum(length(B), length(D))))
odds_quad_pgr/odds_quad_not_pgr # slightly underenriched

# Buccal Hyper, Cerv Hypo (overall) -------
x = dmps[dmps$dir_type.x>0 & dmps$dir_type.y<0,]
head(x)
# Enrichment for PGRs?
# not in PGR    # in PGR
# not in quad       A             B
# in Quad           C             D

# Compute
A <- unique(buccal$cg[!buccal$cg %in% c(x$cg, rownames(pgr_cpgs))])
B <- unique(buccal$cg[buccal$cg %in% rownames(pgr_cpgs) & !buccal$cg %in% x$cg])
C <- unique(x$cg[!x$cg %in% rownames(pgr_cpgs)])
D <- unique(x$cg[x$cg %in% rownames(pgr_cpgs)])

tab <- matrix(c(length(A), length(C), length(B), length(D)), ncol = 2)
epitools::oddsratio(tab)

# Manual:
odds_quad_not_pgr <- length(C)/(sum(length(C), length(A)))
odds_quad_pgr <- length(D)/((sum(length(B), length(D))))
odds_quad_pgr/odds_quad_not_pgr
