source('/Users/chiara/EUTOPS Team Dropbox/Chiara Herzog/eca/Buccal BC index/repository/WID-buccal/0-data/src/estimateDeltaBeta.R')

estimateDeltaBeta(pheno = '/Users/chiara/EUTOPS Team Dropbox/Chiara Herzog/eca/Buccal BC index/repository/WID-buccal/1-ewas/1-buccal/pheno.Rdata',
beta = '/Users/chiara/Dropbox/data/3c-buccal/beta_merged.Rdata',
output = '/Users/chiara/EUTOPS Team Dropbox/Chiara Herzog/eca/Buccal BC index/repository/WID-buccal/1-ewas/1-buccal',
adjustment = c('age', 'ic'))
