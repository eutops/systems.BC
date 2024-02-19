manhattan <- function(dmp_buccal,
                      label_cutoff = 6,
                      sample = "buccal",
                      cols = pal.jama){
  
  require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) 
  require(ggtext)
  require(ggrepel)
  
  db <- as.data.frame(dmp_buccal) |> 
    tibble::rownames_to_column('cg') |> 
    dplyr::mutate(padj = p.adjust(type, method = 'BH'),
                  sig = ifelse(type < 0.05, "yes", NA),
                  sig_adj = ifelse(padj < 0.05, 'yes', NA))
  
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19, lociNames = db$cg)
  db$chr <- anno$chr
  db$pos <- anno$pos
  db <- db |> # Rename X to 23 for interim
    dplyr::mutate(chr = gsub("chr", "", chr),
                  chr = ifelse(chr=="X", 23, chr),
                  chr = as.numeric(chr))
  
  # Add cumulative numbers
  cumulative <- db |> 
    group_by(chr) |> 
    summarise(max_bp = max(pos)) |> 
    mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |> 
    dplyr::select(chr, bp_add) |> 
    ungroup()
  
  chrs <- sort(unique(db$chr))
  
  tmp <- db |> 
    left_join(cumulative, by = "chr") |> 
    mutate(bp_cum = as.numeric(pos) + bp_add,
           label = ifelse(-log10(type)>label_cutoff, cg, ""),
           colour = ifelse(sig_adj == 'yes' & !is.na(sig_adj), "significant",
                           ifelse(chr %in% chrs[seq(1, length(chrs), 2)], "chrs1", "chrs2")),
           colour = factor(colour, levels = c("chrs1",
                                              "chrs2",
                                              'significant')),
           size = ifelse(colour == 'significant', 1, 0))
  
  axis_set <- tmp |>  
    group_by(chr) |> 
    summarize(center = mean(bp_cum)) |> 
    ungroup() |> 
    dplyr::mutate(chr = ifelse(chr==23, "X", chr))
  
  ylim <- tmp |> 
    filter(type == min(type)) |> 
    mutate(ylim = abs(floor(log10(type))) + 2) |> 
    pull(ylim)
  
  manhplot <- tmp |> 
    ggplot(
      aes(x = bp_cum, y = -log10(type), 
          colour = colour
      )) +
    geom_point(aes(size = size)) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    scale_color_manual(values = c("black", cols[1],
                                  cols[2]),
                       name = "",
                       labels = c("",
                                  "",
                                  'significant')) +
    scale_size_continuous(range = c(0.3, 1)) +
    labs(x = NULL, 
         y = "-log<sub>10</sub>(p)",
    ) + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      plot.title = element_markdown(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    ) +
    ggrepel::geom_text_repel(aes(label = label),
                             min.segment.length = 0,
                             nudge_x = 100,
                             nudge_y = 2,
                             force = 2,
                             size = 2.6) +
    guides(size = guide_legend(show_guide = FALSE),
           colour = guide_legend(nrow = 1))
  
  return(manhplot)
}
