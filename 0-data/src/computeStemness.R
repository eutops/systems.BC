#' Compute Stemness Index
#' 
#' This function computes the stemness index as defined in the paper 'Machine Learning Identifies Stemness Features Associated with Oncogenic Dedifferentiation'.
#'
#' The stemness index is calculated based on the provided beta matrix and phenotype data.
#'
#' @param beta A matrix of beta values representing DNA methylation data.
#' @param pheno A data frame containing phenotype information, with a column named 'basename' used for matching with column names of beta matrix.
#' @return A modified phenotype data frame with an additional column 'stemness' containing computed stemness index.
#' @details This function performs the following steps:
#'   1. Identifies the common samples between the beta matrix and phenotype data.
#'   2. Extracts the relevant columns from beta matrix and updates the phenotype data accordingly.
#'   3. Loads pre-computed weights (w) from the 'pcbc-stemsig.p219.Rda' file.
#'   4. Checks overlaps between weight names and beta matrix rows.
#'   5. Computes the stemness scores based on the intersection of weights and beta matrix rows.
#'   6. Scales the stemness scores to be between 0 and 1.
#'   7. Returns the stemness scores vector.
#'
#' @references
#' PanCanStem Workflow: https://github.com/BioinformaticsFMRP/PanCanStem_Web/blob/master/Rmardown/mDNAsi.rmd
#'
#' @examples
#' result <- computeStemness(beta, pheno)
#'
#' @seealso
#' \url{https://github.com/BioinformaticsFMRP/PanCanStem_Web/blob/master/Rmardown/mDNAsi.rmd}
#' 
#' @author B. Theeuwes
#'
#' @export
computeStemness <- function(beta, pheno){

  intersect <- intersect(colnames(beta), pheno$basename)
  beta <- beta[,intersect]
  pheno <- pheno[match(intersect, pheno$basename),]
  
  load(here("0-data", "src", "pcbc-stemsig.p219.Rda")) # w
  w <- mm$w
  
  # Check overlaps with coefficients
  intersect <- intersect(names(w), rownames(beta))
  
  # how many are present?
  present <- sum(names(w) %in% rownames(beta))
  percPresent <- (present/length(w))
  
  # Display info
  if(!identical(intersect, names(w))){
    cat("[Stemness index] ", round(percPresent*100, 2), "% of CpGs present (", present, "/", length(w), ")\n", sep = "")
  }
  
  X <- beta[match(intersect, rownames(beta)),]
  w <- w[names(w) %in% rownames(beta)] 
  X <- as.matrix(X)
  ss <- t(w) %*% X
  
  ## Scale the scores to be between 0 and 1
  ss <- ss - min(ss)
  ss <- ss / max(ss)

  # Result
  stemness <- as.numeric(ss)

  return(stemness)
}




