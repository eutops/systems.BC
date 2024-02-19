# SOX2 index
# Author: Chiara Herzog
# Date: 14 Dec 23
# Computes z Score SOX2 index

index_sox2 <- function(beta_input){
  
  load(here("0-data","src", "weights.Rdata"))
  
  
  # Compute intersect of CpGs.
  intersect <- intersect(rownames(beta_input), rownames(beta_signif))
  mu <- mu[match(intersect, names(mu))]
  sigma <- sigma[match(intersect, names(sigma))]
  wx <- wx[match(intersect, names(wx))]
  
  # Subset beta input
  tmp <- beta_input[match(intersect, rownames(beta_input)),]
  
  idx <- apply(tmp, 2, function(x){
    return(sum(wx*((x-mu)/sigma)))
  })
  
  return(idx)
  
}
