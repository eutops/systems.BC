# Delta beta 

estimateDeltaBeta <- function(beta, pheno, output, adjustment = c("age"),
                              typevar = "type", base = "Control"){
  
  # Adjustment needs a minimum of type and ic
  # need to mention what the type variable and base level are
  
  cat('Beginning delta-beta estimation script...\n\n')
  cat('Current time:',as.character(Sys.time()),'\n')
  cat('cwd:',getwd(),'\n\n')
  
  # Install packages
  for (i in c('dplyr', 'tibble', 'ggplot2', 'stringr', 'BiocManager')){
    if(length(find.package(i, quiet = T)) == 0){
      cat('Installing', i)
      install.packages(i); library(i, character.only = T)
    }
  }
  
  if(length(find.package("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quiet = T)) == 0){
    cat('Installing', i)
    BiocManager::install(i); library(i, character.only = T)
  }
  
  # Load libs
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(stringr))
  
  # Check directories and files exist
  if(!file.exists(beta)){
    stop('cannot read file at path_to_beta')
  }
  
  if(!file.exists(pheno)){
    stop('cannot read file at path_to_pheno')
  }
  
  if(!dir.exists(output)){
    dir.create(output)
    cat('output directory has been created\n')
  }
  
  # load data ------
  cat('Loading beta matrix and phenotype file...')
  if(grepl(".Rdata", beta)){
    beta_tmp <- load(file = beta)
    beta <- get(beta_tmp);rm(beta_tmp)
    pheno_tmp <- load(file = pheno)
    pheno <- get(pheno_tmp);rm(pheno_tmp)
  } else {
    beta <- readRDS(beta)
    pheno <- readRDS(pheno)
  }
  cat(' done\n\n')
  
  # Get overlapping ids in beta and pheno
  intersect <- intersect(colnames(beta), rownames(pheno))
  beta <- beta[,intersect];gc()
  pheno <- pheno[intersect,]
  # ind <- match(colnames(beta), rownames(pheno))
  # pheno <- pheno[ind,]
  
  cat('Statistical adjustment will be made for:\n')
  
  adjustment <- c("type", "ic", adjustment)
  adjustment <- unique(adjustment[adjustment!=""])
  for(i in adjustment){cat(i,'\n')}
  cat('\n')
  
  # Set factor levels correctly
  pheno <- pheno |> 
    dplyr::mutate(type = case_when(get(typevar) == base ~ "Control",
                                   get(typevar) != base ~ "Case"),
                  type = factor(type, levels = c("Control", "Case")),
                  ic = as.numeric(ic)) |> 
    dplyr::select(all_of(adjustment))
  
  if("age" %in% adjustment){
    pheno <- pheno |> 
      mutate(age = as.numeric(age))
  }
  
  # fit linear models ----
  cat('Estimating delta beta...\n')
  db <- matrix(NA, nrow=nrow(beta),ncol=3+ncol(pheno))
  rownames(db) <- rownames(beta)
  colnames(db) <- c('db_epithelial','db_immune',colnames(pheno), "dir_type")
  
  ic.control <- pheno$ic[which(pheno$type=='Control')]
  ic.case <- pheno$ic[which(pheno$type!='Control')]
  
  pB <- txtProgressBar(min=1,max=nrow(beta), width =50L, style = 3)
  for (i in 1:nrow(beta)){
    setTxtProgressBar(pB, i)
    
    ldat <- data.frame(beta=as.numeric(beta[i,]),
                       pheno)
    
    lfit <- lm(beta ~ ., data=ldat)
    db[i,3:(ncol(db)-1)] <- summary(lfit)$coefficients[2:(ncol(pheno)+1),4] # pval
    db[i, ncol(db)] <- summary(lfit)$coefficients[2,1]  # overall directionality
    
    beta.control <- as.numeric(beta[i,which(pheno$type=='Control')])
    beta.case <- as.numeric(beta[i,which(pheno$type!='Control')])
    
    dat.control <- data.frame(beta=beta.control, ic=ic.control)
    dat.case <- data.frame(beta=beta.case, ic=ic.case)
    
    fit.control <- lm(beta ~ ic, data=dat.control)
    fit.case <- lm(beta ~ ic, data=dat.case)
    
    delta_beta_epithelial <- round(fit.case$coefficients[1]-fit.control$coefficients[1],digits=3)
    delta_beta_immune <- round(-fit.control$coefficients[1] - fit.control$coefficients[2]
                               +fit.case$coefficients[1] + fit.case$coefficients[2],digits=3)
    db[i,1] <- delta_beta_epithelial
    db[i,2] <- delta_beta_immune
    
  }
  close(pB)
  cat('done\n\n')
  
  saveRDS(db, file=paste(output,"/delta-beta.Rds",sep=""))
  
  # plot top 100 epithelial delta-betas -----  
  cols <- c("#83A58C", "#DF8A5A")
  
  cat('Plotting top 100 epithelial CpGs...')
  
  tmp <- db |> 
    as.data.frame() |> 
    arrange(desc(abs(db_epithelial))) |> 
    dplyr::slice(1:100) |> 
    tibble::rownames_to_column("cg") |> 
    dplyr::rename(p_type = type,
                  p_ic = ic)
  
  beta_tmp <- data.frame(t(beta[tmp$cg,match(rownames(pheno), colnames(beta))]))
  # identical(colnames(beta_tmp), tmp$cg)
  # identical(rownames(beta_tmp), rownames(pheno))
  
  pheno_tmp <- cbind(pheno, beta_tmp) |> 
    tidyr::pivot_longer(all_of(tmp$cg),
                        names_to = "cg",
                        values_to = "beta") |> 
    dplyr::full_join(tmp, by = "cg") 
  
  # Annotate gene
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, lociNames = tmp$cg) |> 
    data.frame() |> 
    mutate(gene = stringr::str_split(UCSC_RefGene_Name, ";", simplify = T)[,1]) |> 
    tibble::rownames_to_column("cg") |> 
    select(cg, gene)
  
  pheno_tmp <- pheno_tmp |> 
    full_join(anno) |> 
    mutate(cg_label = case_when(gene != "" ~ paste0(cg, " (", gene,
                                                    ")\ndb_imm = ", signif(db_epithelial, 3),
                                                    ",\np = ", signif(p_type, 3)),
                                TRUE ~ paste0(cg,
                                              "\ndb_imm = ", signif(db_epithelial, 3),
                                              ",\np = ", signif(p_type, 3))))
  
  pdf(file=paste(output,"/top-epithelial-db.pdf",sep=""), width=4,height=4.5)
  
  for (i in 1:100){
    plot <- pheno_tmp |>
      dplyr::filter(cg == tmp$cg[i]) |> 
      ggplot(aes(x = ic,
                 y = beta,
                 colour = type)) +
      geom_point(size = 0.75) +
      geom_smooth(method = "lm",
                  se = F) +
      facet_wrap(~ cg_label) +
      scale_colour_manual(values = cols,
                          name = "") +
      theme_bw() +
      theme(strip.text = element_text(hjust = 0))
    
    if(i == 1){
      plot <- plot +
        theme(legend.position = "top")
    } else {
      plot <- plot +
        theme(legend.position = "none")
    }
    
    print(plot)
  }
  
  dev.off()
  
  cat(' done\n')
  
  # plot top 100 immune delta-betas -----  
  
  cat('Plotting top 100 immune CpGs...')
  
  tmp <- db |> 
    as.data.frame() |> 
    arrange(desc(abs(db_immune))) |> 
    dplyr::slice(1:100) |> 
    tibble::rownames_to_column("cg") |> 
    dplyr::rename(p_type = type,
                  p_ic = ic)
  
  beta_tmp <- data.frame(t(beta[tmp$cg,match(rownames(pheno), colnames(beta))]))
  # identical(colnames(beta_tmp), tmp$cg)
  # identical(rownames(beta_tmp), rownames(pheno))
  
  pheno_tmp <- cbind(pheno, beta_tmp) |> 
    tidyr::pivot_longer(all_of(tmp$cg),
                        names_to = "cg",
                        values_to = "beta") |> 
    dplyr::full_join(tmp, by = "cg") 
  
  # Annotate gene
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, lociNames = tmp$cg) |> 
    data.frame() |> 
    mutate(gene = stringr::str_split(UCSC_RefGene_Name, ";", simplify = T)[,1]) |> 
    tibble::rownames_to_column("cg") |> 
    select(cg, gene)
  
  pheno_tmp <- pheno_tmp |> 
    full_join(anno) |> 
    mutate(cg_label = case_when(gene != "" ~ paste0(cg, " (", gene,
                                                    ")\ndb_imm = ", signif(db_immune, 3),
                                                    ",\np = ", signif(p_type, 3)),
                                TRUE ~ paste0(cg,
                                              "\ndb_imm = ", signif(db_immune, 3),
                                              ",\np = ", signif(p_type, 3))))
  
  pdf(file=paste(output,"/top-immune-db.pdf",sep=""), width=4,height=4.5)
  
  for (i in 1:100){
    plot <- pheno_tmp |>
      dplyr::filter(cg == tmp$cg[i]) |> 
      ggplot(aes(x = ic,
                 y = beta,
                 colour = type)) +
      geom_point(size = 0.75) +
      geom_smooth(method = "lm",
                  se = F) +
      facet_wrap(~ cg_label) +
      scale_colour_manual(values = cols,
                          name = "") +
      theme_bw() +
      theme(strip.text = element_text(hjust = 0))
    
    if(i == 1){
      plot <- plot +
        theme(legend.position = "top")
    } else {
      plot <- plot +
        theme(legend.position = "none")
    }
    
    print(plot)
  }
  
  dev.off()
  
  # plot top p-ranked ----------
  cat('Plotting top 100 significant CpGs...')
  
  tmp <- db |> 
    as.data.frame() |> 
    arrange(type) |> 
    dplyr::slice(1:100) |> 
    tibble::rownames_to_column("cg") |> 
    dplyr::rename(p_type = type,
                  p_ic = ic)
  
  beta_tmp <- data.frame(t(beta[tmp$cg,match(rownames(pheno), colnames(beta))]))
  # identical(colnames(beta_tmp), tmp$cg)
  # identical(rownames(beta_tmp), rownames(pheno))
  
  pheno_tmp <- cbind(pheno, beta_tmp) |> 
    tidyr::pivot_longer(all_of(tmp$cg),
                        names_to = "cg",
                        values_to = "beta") |> 
    dplyr::full_join(tmp, by = "cg") 
  
  # Annotate gene
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, lociNames = tmp$cg) |> 
    data.frame() |> 
    mutate(gene = stringr::str_split(UCSC_RefGene_Name, ";", simplify = T)[,1]) |> 
    tibble::rownames_to_column("cg") |> 
    select(cg, gene)
  
  pheno_tmp <- pheno_tmp |> 
    full_join(anno) |> 
    mutate(cg_label = case_when(gene != "" ~ paste0(cg, " (", gene,
                                                    ")\np = ", signif(p_type, 3)),
                                TRUE ~ paste0(cg,
                                              "\np = ", signif(p_type, 3))))
  
  pdf(file=paste(output,"/top-p.pdf",sep=""), width=4,height=4.5)
  
  for (i in 1:100){
    plot <- pheno_tmp |>
      dplyr::filter(cg == tmp$cg[i]) |> 
      ggplot(aes(x = ic,
                 y = beta,
                 colour = type)) +
      geom_point(size = 0.75) +
      geom_smooth(method = "lm",
                  se = F) +
      facet_wrap(~ cg_label) +
      scale_colour_manual(values = cols,
                          name = "") +
      theme_bw() +
      theme(strip.text = element_text(hjust = 0))
    
    if(i == 1){
      plot <- plot +
        theme(legend.position = "top")
    } else {
      plot <- plot +
        theme(legend.position = "none")
    }
    
    print(plot)
  }
  
  dev.off()
  
  
  cat(' done\n\n')
  cat('Session info:\n\n')
  sessionInfo()
}
