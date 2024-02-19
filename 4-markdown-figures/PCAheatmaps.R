library(ComplexHeatmap)

colour.green <- "#1BA68C"
colour.blue <- "#3B5C8D"
colour.red <- "#A61B35"

# Function to subset and get top variable CpGs
topVarCpGs <- function(pheno, beta, top_sites = 0.10, vars = features){
  
  tmp <- pheno[,vars]
  beta <- as.matrix(beta)
  
  # sd
  sd <- matrixStats::rowSds(beta)
  names(sd) <- rownames(beta)
  sd <- sd[order(sd, decreasing = T)]
  sd_top <- sd[1:round(length(sd)*top_sites)]
  beta_top <- beta[names(sd_top),]
  
  # pc
  pc <- prcomp(t(beta_top), scale. = T, center = T)
  
  # append info
  dat <- cbind(tmp, pc$x[,1:10])
  dat$variance <- variance_explained
  
  return(dat)
}
pcHeatmap <- function(pc_dat){
  
  pc_dat$set <-NULL
  pc_dat$sampletype <- NULL
  pc_dat$OCP_ever <- NULL
  pc_dat$cancer_stage <- NULL
  pc_dat$time_to_event <- NULL
  
  pcs <- pc_dat |> 
    dplyr::select(PC1:PC10)
  
  tmp <- pc_dat |> 
    dplyr::select(any_of(features), 'ic', 'hepidish_Neutro') |> 
    dplyr::mutate_if(is.character, as.factor)
  
  mat <- matrix(ncol = ncol(pcs),
                nrow = ncol(tmp))
  colnames(mat) <- colnames(pcs)
  rownames(mat) <- colnames(tmp)
  
  for (i in colnames(mat)){
    print(i)
    for (j in rownames(mat)){
      print(j)
      
      if(is.numeric(tmp[,j])){
        mat[j,i] <- cor.test(pcs[,i], tmp[,j])$p.value
        # pmat[j,i] <- cor.test(pcs[,i], tmp[,j])
      } else {
        mat[j,i] <- kruskal.test(pcs[,i], tmp[,j])$p.value
      }
      
    }
  }
  
  return(mat)
}


## BUCCAL

# Variables
setname = 'discovery_buccal'
sampletypename = 'buccal'

# Load data
load(here("0-data", "data.Rdata"))
load(here(db_path, "data", "3c-buccal", "beta_merged.Rdata"))


# Pheno buccal
pheno_buccal <- droplevels(data[data$set == setname  & data$sampletype == sampletypename,]) |> 
  dplyr::mutate(type = factor(type, levels = c("Control", "BC case")))

features <-  colnames(pheno_buccal)

# Beta buccal
beta_buccal <- beta_merged[,pheno_buccal$basename]
pc_buccal <- topVarCpGs(pheno_buccal, beta_buccal, vars = features)


## CERVICAL

setname = 'discovery_cervical'
sampletypename = 'cervical'

# Load data
load(here("0-data", "data.Rdata"))
load(here(db_path, "data", "3c", "beta_merged.Rdata"))

# Pheno cervical
pheno_cervical <- droplevels(data[data$set == setname & data$sampletype == sampletypename,]) |> 
  dplyr::mutate(type = factor(type, levels = c("Control", "BC case")))

features <-  colnames(pheno_cervical)

# Beta cervical
beta_cervical <- beta_merged[,pheno_cervical$basename]
pc_cervical <- topVarCpGs(pheno_cervical, beta_cervical, vars = features)


## BLOOD
setname = 'discovery_blood'
sampletypename = 'blood'

# Load data
load(here("0-data", "data.Rdata"))
load(here(db_path, "data", "3c-blood", "beta_merged.Rdata"))

# Pheno cervical
pheno_blood <- droplevels(data[data$set == setname & data$sampletype == sampletypename,]) |> 
  dplyr::mutate(type = factor(type, levels = c("Control", "BC case")))

features <-  colnames(pheno_blood)

# Beta blood
beta_blood <- beta_merged[,pheno_blood$basename]
pc_blood <- topVarCpGs(pheno_blood, beta_blood, vars = features)


#=====

pcHeatmap_buccal <- pcHeatmap(pc_buccal)  
pcHeatmap_cervical <- pcHeatmap(pc_cervical)  
pcHeatmap_blood <- pcHeatmap(pc_blood)  


combineHeatmaps <- function(pcHeatmap_blood, pcHeatmap_buccal, pcHeatmap_cervical){
  cols <- c("#1b69a1", "#48a0af", "#71b5a9", "#ec6669", "#f39668", "#bd647d", "#832c9b", "#5f70a8")
  map <- cbind(pcHeatmap_blood,
               pcHeatmap_buccal,
               pcHeatmap_cervical)
  groups <- rep(c("blood", "buccal", "cervical"), each = 10)
  map <- apply(map, 2, function(t) ifelse(t < 0.05, t, NA))
  
  # Selection of features to show
  
  variables <- c("age", "type", "ic", "age_at_menarche", "menopause", "OCP_current", "HRT_total_years", "hepidish_Neutro", "stage_N")
  labs <- c("Age", "Type (BC / Control)", "Immune Cell proportion", "Age at menarche", "Menopause (pre/post)", "OCP current", "HRT (total years)", "Neutrophil proportion", "Nodal stage")
  
  map <- map[variables,]

  
  

  p <- Heatmap(-log10(map),
               row_labels = labs,
               cluster_columns = F,
               # cell_fun = function(j, i, x, y, width, height, fill) {
               #   grid.text(pmap[i, j], x, y, gp = gpar(fontsize = 10))
               # },
               row_names_side = 'left',
               cluster_rows = F,
               column_split = groups,
               column_title = NULL,
               show_row_dend = F, show_column_dend = F,
               #name = '-log10(pvalue)',
               na_col = 'white',
               col = circlize::colorRamp2(breaks = seq(40, 1.3, length.out = 5),
                                          # colors = rev(viridis::viridis(5)),
                                          # colors = rev(color("batlow")(5)),
                                          colors = rev(cols[c(2,5,4, 6, 7)])
               ),
               row_names_gp = grid::gpar(fontsize = 9),
               column_names_gp = grid::gpar(fontsize = 9),
               border_gp = gpar(lwd = 0.5),
               border = T,
               top_annotation = HeatmapAnnotation(grp = anno_block(gp = gpar(fill = c(colour.red,colour.green,colour.blue),
                                                                             lwd = 0.5),
                                                                   labels = c("Blood", "Buccal", "Cervical"),
                                                                   labels_gp = gpar(col = "white",
                                                                                    fontsize = 10,
                                                                                    fontface = "bold"))))
  
  return(p)
}

p <- combineHeatmaps(pcHeatmap_blood, pcHeatmap_buccal, pcHeatmap_cervical)
draw(p)




pca_blood <- ggplot(data = pc_blood, aes(x = PC1, y = PC2, color = hepidish_Neutro, shape = type)) +
  geom_point(alpha = 1, size = 2) +  # Plot points
  labs(
    y = "PC2",
    x = "PC1",
    color = "Neutrophil\nproportion",  # Rename color legend title
    shape = "",  # Rename shape legend title
    title = 'Blood Discovery Set'
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +  # Remove legend
  theme(panel.background = element_blank()) +
  scale_color_gradient(low = "#00E3FF", high = "#FF1C00") + 
  scale_shape_manual(values = c(1, 16))



pca_buccal<- ggplot(data = pc_buccal, aes(x = PC1, y = PC2, color = ic, shape = type)) +
  geom_point(alpha = 1, size = 2) +  # Plot points
  labs(
    y = "PC2",
    x = "PC1",
    color = "Immune Cell\nproportion",  # Rename color legend title
    shape = "",  # Rename shape legend title
    title = 'Buccal Discovery Set'
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +  # Remove legend
  theme(panel.background = element_blank()) +
  scale_color_gradient(low = "#00E3FF", high = "#FF1C00", limits = c(0, 1)) + 
  scale_shape_manual(values = c(1, 16))

pca_buccal

pca_cervical<- ggplot(data = pc_cervical, aes(x = PC1, y = PC2, color = ic, shape = type)) +
  geom_point(alpha = 1, size = 2) +  # Plot points
  labs(
    y = "PC2",
    x = "PC1",
    color = "Immune Cell\nproportion",  # Rename color legend title
    shape = "",  # Rename shape legend title
    title = 'Cervical Discovery Set'
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +  # Remove legend
  theme(panel.background = element_blank()) +
  scale_color_gradient(low = "#00E3FF", high = "#FF1C00", limits = c(0, 1)) + 
  scale_shape_manual(values = c(1, 16))

pca_cervical

#==
# Plot
layout <- "
AAA
AAA
BCD
EEE
"


pdf_width <- 8.3  # Width in inches (a4)
pdf_height <- 8   # Height in inches (a4: 11.7)

grob = grid.grabExpr(draw(p)) 
compiled <- wrap_elements(grob)  + pca_blood +pca_buccal + pca_cervical + guide_area()

compiled <-  compiled + plot_annotation(tag_levels = 'a') + 
  plot_layout(design = layout, guides = 'collect')

ggsave(file.path(outPath,"ED3.pdf"), width = pdf_width, height = pdf_height, compiled)
ggsave(file.path(jpgPath,"ED3.jpg"), width = pdf_width, height = pdf_height, compiled)
