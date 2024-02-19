plotDMR_comparison <- function(ranges1, ranges2,
                               beta1, beta2,
                               pheno1, pheno2,
                               label1 = 'buccal', label2 = 'cervical'){
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  # Load in data + cols
  data <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
                        lociNames = intersect(rownames(beta1), rownames(beta2)))

  # Ranges
  min <- start(ranges1)
  max <- end(ranges1)
  chr <- as.character(seqnames(ranges1))
  c <- data[data$chr==chr & data$pos >= min & data$pos <= max,]
    
  # Data 1
  tmp <- beta1[match(rownames(c), rownames(beta1)),]
  tmp_pheno1 <- cbind(pheno1, t(tmp))
  tmp_pheno1 <- tmp_pheno1 |>
      tidyr::pivot_longer(rownames(tmp), names_to = "cpg", values_to = "beta") |> 
      dplyr::mutate(type2 = dplyr::case_when(type == "Control" ~ "Control",
                                             type != "Control" ~ "Case")) |> 
      dplyr::mutate(type2 = factor(type2, levels = c("Control", "Case")),
                    label = label1)
    
  # Data 2
  tmp <- beta2[match(rownames(c), rownames(beta2)),]
  tmp_pheno2 <- cbind(pheno2, t(tmp))
  tmp_pheno2 <- tmp_pheno2 |>
    tidyr::pivot_longer(rownames(tmp), names_to = "cpg", values_to = "beta") |> 
    dplyr::mutate(type2 = dplyr::case_when(type == "Control" ~ "Control",
                                           type != "Control" ~ "Case")) |> 
    dplyr::mutate(type2 = factor(type2, levels = c("Control", "Case")),
                  label = label2)
  
  tmp_pheno <- rbind(tmp_pheno1, tmp_pheno2)
  ind <- match(tmp_pheno$cpg, rownames(c))
  tmp_pheno$pos <- c$pos[ind]
  t <- tmp_pheno |> 
      dplyr::group_by(label, type2, cpg) |> 
      dplyr::reframe(mean = mean(beta),
                     pos = pos) |> 
      dplyr::distinct() 
    
    # Plot linked ranges
    (a <- t |>
        ggplot(aes(x = pos,
                   y = mean,
                   colour = type2)) +
        geom_point() +
        geom_line(aes(linetype = label)) +
        # ggtitle(paste0(chr, ":", min, "-", max, ",\n",
        #                ranges1$overlapping.genes, "\np=", signif(ranges1$Stouffer,2), " (", label1, "), ", signif(ranges2$Stouffer, 2), " (", label2, ")")) +
        theme_bw() +
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1)) +
        ylim(0, 1) +
        scale_linetype_manual(values = c('dashed', 'solid'),
                              name = '') +
        scale_colour_manual(values = c(colour.Control, colour.Breast),
                            name = "") +
        xlab(paste0(unique(c$chr))) +
        ylab(paste0("mean beta")))
    
    
    return(a)
}
