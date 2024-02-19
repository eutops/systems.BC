source('C:/Users/elisa/EUTOPS Team Dropbox/Bente Theeuwes/eca/Buccal BC index/repository/WID-buccal/0-data/src/estimateDeltaBeta.R')

estimateDeltaBeta(pheno = 'C:/Users/elisa/Documents/data local/GSE225845/temp.Rdata',
beta = 'C:/Users/elisa/Documents/data local/GSE225845/qc/beta_merged.Rdata',
output = 'C:/Users/elisa/EUTOPS Team Dropbox/Bente Theeuwes/eca/Buccal BC index/repository/WID-buccal/1-ewas/4-GSE225845',
adjustment = c('age', 'ic'))
