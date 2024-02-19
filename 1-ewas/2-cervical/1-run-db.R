source('/Users/chiara/EUTOPS Team Dropbox/Chiara Herzog/eca/Buccal BC index/repository/WID-buccal/0-data/src/estimateDeltaBeta.R')

estimateDeltaBeta(pheno = '/Users/chiara/EUTOPS Team Dropbox/Chiara Herzog/eca/Buccal BC index/repository/WID-buccal/1-ewas/2-cervical/pheno.Rdata',
beta = '/Users/chiara/Dropbox/data/3c/beta_merged.Rdata',
output = '/Users/chiara/EUTOPS Team Dropbox/Chiara Herzog/eca/Buccal BC index/repository/WID-buccal/1-ewas/2-cervical',
adjustment = c('age', 'ic'))
