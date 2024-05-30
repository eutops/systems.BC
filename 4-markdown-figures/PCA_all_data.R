# Buccal discovery
# Buccal validation 
# BRCA buccal
# Cervical discovery
# Cervical validation
# BRCA cervical
# Blood discovery
# Blood validation (Wang GSE237036)
# BRCA blood

# Breast data GSE225845 CHECKED
# Breast data tissue at risk CHECKED
# Breast data TCGA CHECKED

# Colour by dataset, shape by tissue type/sample. We could try the top 5% or 10% of variable and overlapping CpGs between EPIC + 450k?

library(here)
library(ggplot2)
library(EpiDISH)

# Functions ---------------------------------------------------------------

# Check overlaps basename - returns beta
subsetData <- function(beta_merged,data){
  intersect <- intersect(colnames(beta_merged), data$basename)
  
  beta <- beta_merged[,intersect]
  return(beta)
}

# Check overlaps basename - returns pheno
subsetPheno <- function(beta_merged,data){
  intersect <- intersect(colnames(beta_merged), data$basename)
  
  pheno <- data[match(intersect, data$basename),]
  return(pheno)
}



# Function to compute top x variable CpGs
topVarCpGs <- function(pheno, beta, top_sites = 0.20){
  start_time <- Sys.time()
  beta <- as.matrix(beta)
  
  end_time <- Sys.time()
  cat("Time elapsed for dataframe to matrix operation:\n")
  print(end_time - start_time)
  
  # sd
  start_time <- Sys.time()
  sd <- matrixStats::rowSds(beta)
  names(sd) <- rownames(beta)
  sd <- sd[order(sd, decreasing = T)]
  sd_top <- sd[1:round(length(sd)*top_sites)]
  beta_top <- beta[names(sd_top),]
  
  end_time <- Sys.time()
  cat("Time elapsed for sd calculation:\n")
  print(end_time - start_time)
  
  # pc
  start_time <- Sys.time()
  pc <- prcomp(t(beta_top), scale. = T, center = T)
  end_time <- Sys.time()
  cat("Time elapsed for PCA:\n")
  print(end_time - start_time)
  
  # append info
  start_time <- Sys.time()
  dat <- cbind(pheno, pc$x[,1:10])
  end_time <- Sys.time()
  cat("Time elapsed for binding PC rows:\n")
  print(end_time - start_time)
  
  summ <- summary(pc)
  
  dat$variance <- summ$importance[2,]
  
  return(dat)
}


# Load data ---------------------------------------------------------------

load(here("0-data", "data.Rdata"))
db_path <- "C:/Users/elisa/EUTOPS Team Dropbox/Bente Theeuwes"
dataDir <- "C:/Users/elisa/Documents/data local" # Folder where data is stored locally

## 1 - Buccal validation -----------------------------------------------------------------------

load(here(db_path, "data", "4c-buccal", "beta_merged.Rdata"))
beta.buccal_validation <- subsetData(beta_merged,data)
rm(beta_merged)

## 2 - Cervical validation -----------------------------------------------------------------------

load(here(db_path, "data", "3c", "beta_merged.Rdata"))
beta.cervical_validation <- subsetData(beta_merged,data)

rm(beta_merged)
## 3 - Blood validation (Wang; GSE237036) -----------------------------------------------------------------------

load(here(db_path, "data", "blood", "GSE237036_wang", "beta_merged.Rdata"))
load(here(db_path, "data", "blood", "GSE237036_wang", "GSE237036_pheno_trimmed.Rdata"))
beta.blood_validation <- subsetData(beta_merged,pheno_trimmed)

pheno_trimmed$subset <- "GSE237036"
pheno_trimmed$sampletype <- "blood"

load(file.path(db_path,"eca/bc-blood/Manuscript/0-data/cent12CT.m.rda")) # epiDISH reference matrix

# EpiDISH 
data(cent12CT.m)
frac.m <- hepidish(beta.m = beta_merged,
                   ref1.m = centEpiFibIC.m,
                   ref2.m = cent12CT.m,
                   h.CT.idx = 3,
                   method = 'RPC',
                   maxit = 500)

pheno_trimmed <- cbind(pheno_trimmed, frac.m)


data <-merge(data,pheno_trimmed,by = intersect(names(data), names(pheno_trimmed)),all=T)


rm(beta_merged)
## 4 - BRCA -----------------------------------------------------------------------

load(here(db_path, "data", "BRCA-DS1", "beta_merged.Rdata"))
beta.BRCA_validation <- subsetData(beta_merged,data)

rm(beta_merged)
## 5 - Tissue at Risk -----------------------------------------------------------------------

load(here(db_path, "data", "TISSUE-AT-RISK", "beta_merged.Rdata"))
beta.tar <- subsetData(beta_merged,data)

rm(beta_merged)
## 6 - TCGA -----------------------------------------------------------------------


load(here(db_path, "data", "tcga", "TCGA-BRCA", "beta.Rdata"))

# Load TCGA pheno
load(here(db_path, "eca", "sola", "1-download", "2-output", "data.Rdata"))
data.tcga <- data
load(here("0-data", "data.Rdata"))
intersect <- intersect(data.tcga$barcode1, colnames(beta))
beta <- beta[,match(intersect, colnames(beta))]
data.tcga <- data.tcga[match(intersect, data.tcga$barcode1),]
data.tcga$subset <- "TCGA"
data.tcga$sampletype <- "breast"
data.tcga$basename <- data.tcga$barcode1
beta.tcga <- subsetData(beta,data.tcga)

data<-merge(data,data.tcga,by = intersect(names(data), names(data.tcga)),all=T)

rm(beta)


## 7 - Blood discovery ------------------------------------------------------

load(here(db_path, "data", "3c-blood", "beta_merged.Rdata"))
beta.blood_discovery <- subsetData(beta_merged,data)

rm(beta_merged)


## 8- Breast - GSE225845 ------------------------------------------------------

load(file.path(dataDir,"GSE225845","pheno_trimmed.Rdata"))
load(file.path(dataDir,"GSE225845","qc","beta_merged.Rdata"))

pheno_trimmed$subset <- "GSE225845"
pheno_trimmed$sampletype <- "breast"

beta.GSE225845 <- subsetData(beta_merged,pheno_trimmed)

data <-merge(data,pheno_trimmed,by = intersect(names(data), names(pheno_trimmed)),all=T)

rm(beta_merged)


## 9- Buccal discovery --------------------------------------------------------

load(here(db_path, "data", "3c-buccal", "beta_merged.Rdata"))

beta.buccal_discovery <- subsetData(beta_merged,data)

rm(beta_merged)

# Merge beta matrices by common probes -----------------------------------------------------------------------
merged <- merge(beta.buccal_validation, beta.cervical_validation, by='row.names')
rownames(merged) <- merged[,1]
merged[,1] <- NULL

merged <- merge(merged, beta.blood_validation, by='row.names')
rownames(merged) <- merged[,1]
merged[,1] <- NULL

merged <- merge(merged, beta.BRCA_validation, by='row.names')
rownames(merged) <- merged[,1]
merged[,1] <- NULL

merged <- merge(merged, beta.tar, by='row.names')
rownames(merged) <- merged[,1]
merged[,1] <- NULL

merged <- merge(merged, beta.tcga, by='row.names')
rownames(merged) <- merged[,1]
merged[,1] <- NULL

merged <- merge(merged, beta.blood_discovery, by='row.names')
rownames(merged) <- merged[,1]
merged[,1] <- NULL

merged <- merge(merged, beta.GSE225845, by='row.names')
rownames(merged) <- merged[,1]
merged[,1] <- NULL

merged <- merge(merged, beta.buccal_discovery, by='row.names')
rownames(merged) <- merged[,1]
merged[,1] <- NULL

rm(beta.buccal_validation,beta.cervical_validation, beta.blood_validation,beta.BRCA_validation,beta.tar, beta.tcga,beta.blood_discovery, beta.GSE225845, beta.buccal_discovery) 



# Create master pheno file-------------------------------------------------------------------------
pheno <- subsetPheno(merged,data)


# PCA on all samples ------------------------------------------------------

# Top x CpGs and PCA
pcs <- topVarCpGs(pheno,merged)

# Save dataframe for plotting
save(pcs,file=file.path(getwd(),"0-data/dataframes-plotting/pca_all.Rdata"))


# PCA per sample type -----------------------------------------------------

## Buccal -----------------------------------------------------
pheno.buccal <- pheno[pheno$sampletype=="buccal",]
beta.buccal <- merged[,pheno$sampletype=="buccal"]

# Top x CpGs and PCA
pcs <- topVarCpGs(pheno.buccal,beta.buccal)



# Save dataframe for plotting
save(pcs,file=file.path(getwd(),"0-data/dataframes-plotting/pca_all_buccal.Rdata"))

rm(beta.buccal, pheno.buccal)

# Cervical
pheno.cervical <- pheno[pheno$sampletype=="cervical",]
beta.cervical <- merged[,pheno$sampletype=="cervical"]

# Top x CpGs and PCA
pcs <- topVarCpGs(pheno.cervical,beta.cervical)

# Save dataframe for plotting
save(pcs,file=file.path(getwd(),"0-data/dataframes-plotting/pca_all_cervical.Rdata"))

rm(beta.cervical, pheno.cervical)

# Blood
pheno.blood <- pheno[pheno$sampletype=="blood",]
beta.blood <- merged[,pheno$sampletype=="blood"]

# Top x CpGs and PCA
pcs <- topVarCpGs(pheno.blood,beta.blood)

# Save dataframe for plotting
save(pcs,file=file.path(getwd(),"0-data/dataframes-plotting/pca_all_blood.Rdata"))

rm(beta.blood, pheno.blood)

# Breast
pheno.breast <- pheno[pheno$sampletype=="breast",]
beta.breast <- merged[,pheno$sampletype=="breast"]

# Top x CpGs and PCA
pcs <- topVarCpGs(pheno.breast,beta.breast)

# Save dataframe for plotting
save(pcs,file=file.path(getwd(),"0-data/dataframes-plotting/pca_all_breast.Rdata"))

rm(beta.breast, pheno.breast)




# Visualisation -----------------------------------------------------------
pca_plot <- ggplot(data = pcs, aes(x = PC1, y = PC2, color = sampletype, shape = set)) +
  geom_point(alpha = 1, size = 2) +  # Plot points
 # geom_rug(col=type,alpha=0.1, size=1.5) + 
  labs(
    y = "PC2",
    x = "PC1",
    color = "Sample type",  # Rename color legend title
    shape = "Dataset",  # Rename shape legend title
    title = 'Merge try out'
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +  # Remove legend
  theme(panel.background = element_blank()) +
  scale_shape_manual(values = c(1, 16, 13, 14, 10,11,12,15))


