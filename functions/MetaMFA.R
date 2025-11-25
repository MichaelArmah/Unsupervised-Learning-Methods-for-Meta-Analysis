# MetaMFA Function and Example Data Preparation ----

# load package(s) ----
library(tidyverse)
library(tidymodels)
library(here)
library(FactoMineR)


# load data
McLouth_Data <- readRDS(here("data/McLouth_data.rds"))


# Variable grouping for MUTOS
# M – Methods
methods_vars <- c(
  "ROB_seq_gen_clean", "ROB_alloc_clean",
  "ROB_participant_blinding_clean", "ROB_assessor_blinding_clean",
  "ROB_outcome_reporting_clean", "attrition_overall", "attrition_differential"
)

# U – Units
units_vars <- c(
  "age_M_cent", "female_pct_cent", "nonwhite_pct_cent", "college_pct"
)

# T – Treatment
treatment_vars <- c(
  "trt_arm", "ctl_arm", "intervention_detailed", "intervention_class",
  "control_class", "sessions", "session_length", "total_hours",
  "delivery_format", "in_person", "care_provider", "discipline"
)

# O – Outcomes
outcome_vars <- c(
  "outcome", "outcome_group", "outcome_category", "facit"
)

# S – Settings
settings_vars <- c(
  "cancer_type", "type", "stage", "phase", "focus", "n_weeks", "weeks_cent", "weeks_post"
)

# Ordered MUTOS variables
MUTOS_vars <- c(methods_vars, units_vars, treatment_vars, outcome_vars, settings_vars)

# Group sizes
group_sizes_MUTOS <- c(
  length(methods_vars),
  length(units_vars),
  length(treatment_vars),
  length(outcome_vars),
  length(settings_vars)
)


# MetaMFA function
MetaMFA <- function(data, study_col, yi_col, vi_col, mods, rho, num_pcs, group_assignments, group_sizes, group_type, group_names,...) {
  
  exclude_in_dat <- c("study_col", "yi_col", "vi_col")
  dat <- data %>%
    select(-any_of(exclude_in_dat))%>%
    select(all_of(group_assignments))
  
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
  
  # Calculate total variance
  tot.variance <- v + tau2_hat
  # Make inverse variance weights
  inverse_variance_mv <- 1 / tot.variance
  
  MetaMFA_result <- MFA(
    dat,
    group = group_sizes,
    type = group_type,
    name.group = group_names,
    ncp = num_pcs,
    row.w = inverse_variance_mv,
    ...
  )
  
  list(results = MetaMFA_result, tau2 = tau2_hat)
  
}


# Example function call
MetaMFA_model <- MetaMFA(
  data   = McLouth_Data,
  study = "study_id",
  yi    = "g",
  vi    = "Vg",
  mods  = ~ 1,
  rho   = 0.8,
  num_pcs = 100,
  group_assignments = MUTOS_vars,
  group_sizes = group_sizes_MUTOS,
  group_type = c("m", "s", "m","n", "m"),
  group_names = c("Methods", "Units", "Treatment", "Outcomes", "Settings"),
  graph = FALSE
)

# Example function results
MetaMFA_model$tau2        
MetaMFA_model$results 



