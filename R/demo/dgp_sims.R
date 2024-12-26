# Load necessary library
library(MASS) # For generating multivariate normal data

# 1. Define Stability Check Function
check_stability <- function(A) {
  eigenvalues <- eigen(A)$values
  return(all(Mod(eigenvalues) < 1)) # Stability condition
}

# 2. Define Function to Generate Stable VAR(1) Matrix
generate_var_matrix <- function(m, mu_D, sigma_D, mu_O, sigma_O, off_diag_prob) {
  stable <- FALSE
  while (!stable) {
    # Initialize matrix
    A <- matrix(0, nrow = m, ncol = m)
    
    # Generate diagonal elements
    diag(A) <- rnorm(m, mean = mu_D, sd = sigma_D)
    
    # Generate off-diagonal elements
    for (i in 1:m) {
      for (j in 1:m) {
        if (i != j && runif(1) < off_diag_prob) {
          A[i, j] <- rnorm(1, mean = mu_O, sd = sigma_O)
        }
      }
    }
    
    # Check stability
    stable <- check_stability(A)
  }
  return(A)
}

# 3. Simulate VAR(1) Time Series Data
simulate_var1 <- function(T, m, A, mu_I, sigma_I) {
  # Generate intercepts
  intercept <- rnorm(m, mean = mu_I, sd = sigma_I)
  
  # Initialize time series data
  Y <- matrix(0, nrow = T, ncol = m)
  
  # Generate the first observation
  Y[1, ] <- rnorm(m, mean = 0, sd = 1)
  
  # Generate subsequent observations
  for (t in 2:T) {
    Y[t, ] <- intercept + A %*% Y[t - 1, ] + rnorm(m, mean = 0, sd = 0.001)
  }
  
  return(Y)
}

# 4. Set Simulation Parameters
T_values <- c(50, 100, 150, 200, 250)
m_values <- c(10, 20, 50, 100)
off_diag_prob <- c(0.01, 0.1, 0.8) # Sparse, intermediate, dense probabilities

# Means and standard deviations for intercepts and matrix entries
mu_I <- 0.01
sigma_I <- 0.01
mu_D <- 0.15
sigma_D <- 0.15

# Define densities
sparse <- list(mu_O = 0.3, sigma_O = 0.3, mu_D = 0.3, sigma_D = 0.3)
intermediate <- list(mu_O = 0.1, sigma_O = 0.1, mu_D = 0.15, sigma_D = 0.15)
dense <- list(mu_O = 0.01, sigma_O = 0.01, mu_D = 0.15, sigma_D = 0.15)

# 5. Main Simulation Function
simulate_data <- function(T_values, m_values, density, off_diag_prob) {
  results <- list()
  for (T in T_values) {
    for (m in m_values) {
      # Generate stable VAR(1) coefficient matrix
      A <- generate_var_matrix(m, density$mu_D, density$sigma_D, density$mu_O, density$sigma_O, off_diag_prob)
      
      # Simulate VAR(1) time series data
      Y <- simulate_var1(T, m, A, mu_I, sigma_I)
      
      # Store results
      results[[paste0("T_", T, "_m_", m)]] <- list(A = A, Y = Y)
    }
  }
  return(results)
}

# 6. Simulate Data for Sparse, Intermediate, and Dense Configurations
set.seed(123) # For reproducibility
sparse_results <- simulate_data(T_values, m_values, sparse, off_diag_prob[1])
intermediate_results <- simulate_data(T_values, m_values, intermediate, off_diag_prob[2])
dense_results <- simulate_data(T_values, m_values, dense, off_diag_prob[3])

# 7. Save or Analyze Results
# Example: Save results to an RData file
# save(sparse_results, intermediate_results, dense_results, file = "var1_simulation_results.RData")

sparse_results[["T_50_m_10"]] # Print A and Y for sparse case with T=50 and m=10

heatmap(sparse_results[["T_50_m_10"]]$A)

# Basic heatmap using base R
heatmap(sparse_results[["T_50_m_10"]]$A, 
        Colv = NA, 
        Rowv = NA, 
        scale = "none", 
        col = colorRampPalette(c("blue", "white", "red"))(100),
        main = "Matrix Heatmap")
