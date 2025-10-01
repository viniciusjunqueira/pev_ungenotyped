# ============================================
# Simulation Functions for PEV Analysis
# ============================================

#' Build pedigree-based relationship matrix (A)
#'
#' @param ped Data frame with columns: id, sire, dam
#' @return Additive relationship matrix
build_A_matrix <- function(ped) {
  n <- nrow(ped)
  A <- matrix(0, n, n)
  
  for(i in 1:n) {
    sire <- ped$sire[i]
    dam <- ped$dam[i]
    
    # Diagonal element
    if(sire == 0 && dam == 0) {
      A[i,i] <- 1.0
    } else if(sire == 0 || dam == 0) {
      known_parent <- max(sire, dam)
      A[i,i] <- 1.0 + 0.5 * A[known_parent, known_parent] / 2
    } else {
      A[i,i] <- 1.0 + 0.5 * A[sire, dam]
    }
    
    # Off-diagonal elements
    if(i > 1) {
      for(j in 1:(i-1)) {
        if(sire > 0 && dam > 0) {
          A[i,j] <- 0.5 * (A[j, sire] + A[j, dam])
        } else if(sire > 0) {
          A[i,j] <- 0.5 * A[j, sire]
        } else if(dam > 0) {
          A[i,j] <- 0.5 * A[j, dam]
        }
        A[j,i] <- A[i,j]
      }
    }
  }
  
  return(A)
}

#' Simulate genotypes with Mendelian inheritance
#'
#' @param ped Pedigree data frame
#' @param genotyped_ids Vector of IDs to genotype
#' @param n_markers Number of SNP markers
#' @param n_founders Number of founder animals
#' @return Matrix of genotypes (0, 1, 2)
simulate_genotypes <- function(ped, genotyped_ids, n_markers, n_founders) {
  M <- matrix(0, length(genotyped_ids), n_markers)
  
  for(i in seq_along(genotyped_ids)) {
    id <- genotyped_ids[i]
    
    if(id <= n_founders) {
      # Founders: random genotypes
      M[i,] <- rbinom(n_markers, 2, 0.5)
    } else {
      # Progeny: Mendelian sampling
      sire <- ped$sire[id]
      dam <- ped$dam[id]
      
      sire_allele <- if(sire > 0 && sire %in% genotyped_ids) {
        sire_idx <- which(genotyped_ids == sire)
        rbinom(n_markers, 1, M[sire_idx,] / 2)
      } else {
        rbinom(n_markers, 1, 0.5)
      }
      
      dam_allele <- if(dam > 0 && dam %in% genotyped_ids) {
        dam_idx <- which(genotyped_ids == dam)
        rbinom(n_markers, 1, M[dam_idx,] / 2)
      } else {
        rbinom(n_markers, 1, 0.5)
      }
      
      M[i,] <- sire_allele + dam_allele
    }
  }
  
  return(M)
}

#' Build pedigree with diverse family structures
#'
#' @param n_founders Number of founder animals
#' @param n_sires Number of sires
#' @param n_dams Number of dams
#' @param n_progeny Number of progeny
#' @return Data frame with pedigree structure
build_pedigree <- function(n_founders, n_sires, n_dams, n_progeny) {
  n_total <- n_founders + n_progeny
  
  ped <- data.frame(
    id = 1:n_total,
    sire = 0,
    dam = 0
  )
  
  sire_ids <- 1:n_sires
  dam_ids <- (n_sires + 1):n_founders
  current_id <- n_founders + 1
  
  # Full-sibs (25%)
  n_fullsibs <- floor(n_progeny * 0.25)
  for(i in seq_len(n_fullsibs)) {
    ped$sire[current_id] <- sample(sire_ids, 1)
    ped$dam[current_id] <- sample(dam_ids, 1)
    current_id <- current_id + 1
  }
  
  # Paternal half-sibs (35%)
  n_paternal <- floor(n_progeny * 0.35)
  for(i in seq_len(n_paternal)) {
    ped$sire[current_id] <- sample(sire_ids, 1)
    ped$dam[current_id] <- sample(dam_ids, 1)
    current_id <- current_id + 1
  }
  
  # Maternal half-sibs (25%)
  n_maternal <- floor(n_progeny * 0.25)
  for(i in seq_len(n_maternal)) {
    ped$sire[current_id] <- sample(sire_ids, 1)
    ped$dam[current_id] <- sample(dam_ids, 1)
    current_id <- current_id + 1
  }
  
  # One parent unknown (remaining)
  while(current_id <= n_total) {
    if(runif(1) < 0.5) {
      ped$sire[current_id] <- sample(sire_ids, 1)
    } else {
      ped$dam[current_id] <- sample(dam_ids, 1)
    }
    current_id <- current_id + 1
  }
  
  return(ped)
}

#' Build genomic relationship matrix
#'
#' @param M Centered marker matrix
#' @return Genomic relationship matrix
build_G_matrix <- function(M) {
  p <- colMeans(M) / 2
  M_centered <- scale(M, center = 2*p, scale = FALSE)
  scaling_factor <- 2 * sum(p * (1 - p))
  G <- (M_centered %*% t(M_centered)) / scaling_factor
  
  # Ensure positive definite
  min_eval <- min(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
  if(min_eval < 0.01) {
    G <- G + diag(0.01 - min_eval, nrow(G))
  }
  
  return(G)
}

#' Build H-inverse for single-step GBLUP
#'
#' @param A Pedigree relationship matrix
#' @param G Genomic relationship matrix (blended)
#' @param genotyped_ids IDs of genotyped animals
#' @return H-inverse matrix
build_H_inverse <- function(A, G, genotyped_ids) {
  A_inv <- solve(A)
  A22 <- A[genotyped_ids, genotyped_ids]
  
  G_inv <- solve(G)
  A22_inv <- solve(A22)
  correction <- G_inv - A22_inv
  
  H_inv <- A_inv
  H_inv[genotyped_ids, genotyped_ids] <- 
    H_inv[genotyped_ids, genotyped_ids] + correction
  
  return(H_inv)
}