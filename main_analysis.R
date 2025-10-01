# ============================================
# Main Analysis: PEV Comparison
# ============================================

source("config.R")
source("simulation_functions.R")

library(MASS)

#' Run PEV comparison analysis
#'
#' @param problem_sizes Data frame with problem size specifications
#' @param n_replicates Number of replicates per size
#' @param seed Random seed
#' @param sigma2_e Residual variance
#' @param sigma2_u Genetic variance
#' @param alpha Weight on genomic relationship
#' @param beta Weight on pedigree relationship
#' @return List with timing and PEV results
run_pev_analysis <- function(problem_sizes, n_replicates, seed, 
                             sigma2_e, sigma2_u, alpha, beta) {
  
  set.seed(seed)
  
  all_timing_results <- list()
  all_pev_results <- list()
  individual_pev_list <- list()
  
  for(size_idx in seq_len(nrow(problem_sizes))) {
    size_params <- problem_sizes[size_idx, ]
    
    message(sprintf("\nAnalyzing %s problem size...", size_params$size_label))
    message(sprintf("  Genotyped: %d, Young: %d, Markers: %d", 
                    size_params$n_genotyped, size_params$n_young, 
                    size_params$n_markers))
    
    timing_reps <- data.frame()
    pev_comparison <- data.frame()
    
    for(rep in seq_len(n_replicates)) {
      message(sprintf("  Replicate %d/%d", rep, n_replicates))
      
      set.seed(seed + size_idx * 100 + rep)
      
      tryCatch({
        # Build pedigree and relationships
        ped <- build_pedigree(size_params$n_founders, size_params$n_sires,
                              size_params$n_dams, size_params$n_progeny)
        n_total <- nrow(ped)
        A <- build_A_matrix(ped)
        
        # Define populations
        genotyped_ids <- c(
          1:size_params$n_founders,
          sample((size_params$n_founders + 1):(n_total - size_params$n_young), 
                 size_params$n_genotyped - size_params$n_founders)
        )
        genotyped_ids <- sort(genotyped_ids)
        young_ids <- (n_total - size_params$n_young + 1):n_total
        training_ids <- genotyped_ids
        n_training <- length(training_ids)
        
        # Simulate genotypes and build G
        M <- simulate_genotypes(ped, genotyped_ids, size_params$n_markers, 
                                size_params$n_founders)
        G_genomic <- build_G_matrix(M)
        A22 <- A[genotyped_ids, genotyped_ids]
        G_blended <- alpha * G_genomic + beta * A22
        
        # Build H-inverse
        H_inv <- build_H_inverse(A, G_blended, genotyped_ids)
        lambda <- sigma2_e / sigma2_u
        
        # Partition
        G_inv_tt <- H_inv[training_ids, training_ids] * lambda
        G_inv_ty <- H_inv[training_ids, young_ids] * lambda
        G_inv_yt <- H_inv[young_ids, training_ids] * lambda
        G_inv_yy <- H_inv[young_ids, young_ids] * lambda
        
        # Design matrices
        Z_t <- diag(n_training)
        R_inv <- diag(n_training) / sigma2_e
        X <- matrix(1, n_training, 1)
        
        # Method 1: Direct MME
        time_direct <- system.time({
          C_full <- rbind(
            cbind(t(X) %*% R_inv %*% X, t(X) %*% R_inv %*% Z_t, 
                  matrix(0, ncol(X), size_params$n_young)),
            cbind(t(Z_t) %*% R_inv %*% X, t(Z_t) %*% R_inv %*% Z_t + G_inv_tt, 
                  G_inv_ty),
            cbind(matrix(0, size_params$n_young, ncol(X)), G_inv_yt, G_inv_yy)
          )
          C_inv <- solve(C_full)
          
          # Correct indexing for young animals block
          idx_start <- ncol(X) + n_training + 1
          idx_end <- ncol(X) + n_training + size_params$n_young
          PEV_direct <- sigma2_e * C_inv[idx_start:idx_end, idx_start:idx_end]
        })[3]
        
        # Method 2: Schur complement
        time_schur <- system.time({
          B <- rbind(
            cbind(t(X) %*% R_inv %*% X, t(X) %*% R_inv %*% Z_t),
            cbind(t(Z_t) %*% R_inv %*% X, t(Z_t) %*% R_inv %*% Z_t + G_inv_tt)
          )
          B_inv <- solve(B)
          
          # Correct indexing for training animals block in B_inv
          idx_train_start <- ncol(X) + 1
          idx_train_end <- ncol(X) + n_training
          B22_inv <- B_inv[idx_train_start:idx_train_end, idx_train_start:idx_train_end]
          
          S <- G_inv_yy - G_inv_yt %*% B22_inv %*% G_inv_ty
          PEV_schur <- sigma2_e * solve(S)
        })[3]
        
        # Method 3: Parent average
        PEV_training <- sigma2_e * diag(B_inv)[idx_train_start:idx_train_end]
        acc_training <- sqrt(pmax(0, 1 - PEV_training / sigma2_u))
        acc_map <- numeric(n_total)
        acc_map[training_ids] <- acc_training
        
        PEV_PA <- vapply(young_ids, function(id) {
          acc_sire <- if(ped$sire[id] > 0) acc_map[ped$sire[id]] else 0
          acc_dam <- if(ped$dam[id] > 0) acc_map[ped$dam[id]] else 0
          acc_PA <- sqrt((acc_sire^2 + acc_dam^2) / 4)
          sigma2_u * (1 - acc_PA^2)
        }, numeric(1))
        
        # Store results
        timing_reps <- rbind(timing_reps, data.frame(
          replicate = rep,
          method = c("Direct", "Schur"),
          time_seconds = c(time_direct, time_schur)
        ))
        
        pev_comparison <- rbind(pev_comparison, data.frame(
          replicate = rep,
          mean_pev_direct = mean(diag(PEV_direct)),
          mean_pev_schur = mean(diag(PEV_schur)),
          mean_pev_PA = mean(PEV_PA),
          max_diff_direct_schur = max(abs(PEV_direct - PEV_schur)),
          cor_direct_PA = cor(diag(PEV_direct), PEV_PA),
          cor_schur_PA = cor(diag(PEV_schur), PEV_PA)
        ))
        
        if(rep == 1) {
          individual_pev_list[[size_params$size_label]] <- data.frame(
            animal = seq_len(size_params$n_young),
            PEV_direct = diag(PEV_direct),
            PEV_schur = diag(PEV_schur),
            PEV_PA = PEV_PA,
            size = size_params$size_label
          )
        }
        
      }, error = function(e) {
        warning(sprintf("Error in size %s, replicate %d: %s", 
                        size_params$size_label, rep, e$message))
      })
    }
    
    all_timing_results[[size_params$size_label]] <- timing_reps
    all_pev_results[[size_params$size_label]] <- pev_comparison
  }
  
  return(list(
    timing = all_timing_results,
    pev = all_pev_results,
    individual = individual_pev_list,
    problem_sizes = problem_sizes
  ))
}