# MetaPAM Function ----

# load package(s) ----
library(tidyverse)
library(tidymodels)
library(here)
library(gower)
library(WeightedCluster)

# load data
McLouth_Data <- readRDS(here("data/McLouth_data.rds"))


# MetaPAM function
MetaPAM <- function(data, study_col, yi_col, vi_col, mods, rho, num_clusters, ...) {
  
  dat <- data %>%
    select(-yi_col, -vi_col, -study_col)
  
  # Trace function
  tr <- function(M) sum(diag(M))
  
  # Cleaning
  
  s <- data[[study_col]]
  if (!is.factor(s)) s <- factor(s)
  y <- as.numeric(data[[yi_col]])
  v <- as.numeric(data[[vi_col]])
  X <- model.matrix(mods, data = dat)
  p <- ncol(X)
  m <- nlevels(s)
  
  # Split by study
  split_idx <- split(seq_len(nrow(data)), s)
  
  
  S_TWT <- 0
  S_TWX <- matrix(0, 1, p)
  S_XWX <- matrix(0, p, p)
  
  S_w_over_k_XX <- matrix(0, p, p)
  S_w_over_k_D  <- matrix(0, p, p)
  S_w2_XJX      <- matrix(0, p, p)
  sum_k_w       <- 0
  
  for (idx in split_idx) {
    k_j <- length(idx)
    y_j <- matrix(y[idx], ncol = 1)           # k_j x 1
    X_j <- X[idx, , drop = FALSE]             # k_j x p
    J_j <- matrix(1, k_j, k_j)
    
    vbar_j <- mean(v[idx], na.rm = TRUE)
    w_j    <- 1 / (k_j * vbar_j)              # your working weight (study-constant)
    W_j    <- diag(w_j, k_j, k_j)
    
    S_TWT <- S_TWT + t(y_j) %*% W_j %*% y_j
    S_TWX <- S_TWX + t(y_j) %*% W_j %*% X_j
    S_XWX <- S_XWX + t(X_j) %*% W_j %*% X_j
    
    sum_k_w       <- sum_k_w + k_j * w_j
    S_w_over_k_XX <- S_w_over_k_XX + (w_j / k_j) * (t(X_j) %*% X_j)
    S_w_over_k_D  <- S_w_over_k_D  + (w_j / k_j) * ((t(X_j) %*% J_j %*% X_j) - (t(X_j) %*% X_j))
    S_w2_XJX      <- S_w2_XJX      + (w_j^2)     *  (t(X_j) %*% J_j %*% X_j)
  }
  
  # Q_E calculation
  XWX_inv <- solve(S_XWX)
  Q_E     <- as.numeric(S_TWT - S_TWX %*% XWX_inv %*% t(S_TWX))
  
  # method of moments estimator fortau^2 split by numerator & denominator
  V   <- XWX_inv
  num <- Q_E - m + tr(V %*% S_w_over_k_XX) + rho * tr(V %*% S_w_over_k_D)
  den <- sum_k_w - tr(V %*% S_w2_XJX)
  
  tau2_hat <- max(0, as.numeric(num / den))
  
  # Calculate total Variance variance
  tot.variance <- v + tau2_hat
  # Make inverse variance weights
  inverse_variance_mv <- 1 / tot.variance
  
  gower.dist_mv <- daisy(dat, metric = "gower")
  
  
  MetaPAM_result <- wcKMedoids(gower.dist_mv, k = num_clusters, weights = inverse_variance_mv,
                               method="PAMonce",...)
  
  list(results = MetaPAM_result, tau2 = tau2_hat, distances = gower.dist_mv)
  
}


# Function call
MetaPAM_clusters <- MetaPAM(
  data   = McLouth_Data,
  study = "study_id",
  yi    = "g",
  vi    = "Vg",
  mods  = ~ 1,
  rho   = 0.8,
  num_clusters = 9,
  cluster.only = FALSE,
  debuglevel=0,
  npass = 1,
  initialclust=NULL
)

# Function results
MetaPAM_clusters$tau2
MetaPAM_clusters$results
MetaPAM_clusters$distances
