# Code to produce figures used in the paper ----

# load package(s) ----
library(tidyverse)
library(tidymodels)
library(here)
library(factoextra)
library(FactoMineR)
library(gower)
library(patchwork)


McLouth_Data <- readRDS(here("data/McLouth_data.rds"))
McLouth_Covariates <- McLouth_Data %>%
  select(-g, -study_id, -Vg)

# Meta FAMD function call
MetaFAMD_model <- MetaFAMD(
  data   = McLouth_Data,
  study = "study_id",
  yi    = "g",
  vi    = "Vg",
  mods  = ~ 1,
  rho   = 0.8,
  num_pcs = 100,
  graph = FALSE
)

# MetaMFA function call
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

# Function results
MetaFAMD_results <- MetaFAMD_model$results
MetaMFA_results <- MetaMFA_model$results 



# MetaFAMD Scree Plot
fviz_screeplot(MetaFAMD_results) + 
  labs(
    title = "MetaFAMD Scree Plot",
    x = "Principal Components") +
  theme(
    axis.text = element_text(size = 16),        # Axis labels
    axis.title = element_text(size = 18),       # Axis titles
    plot.title = element_text(size = 20, face = "bold"),  # Plot title
    legend.text = element_text(size = 14),      # Legend text
    legend.title = element_text(size = 16)      # Legend title
  )



# MetaFAMD - Top 10 Contributing Variables
fviz_famd_var(
  MetaFAMD_results,
  choice = "var",
  col.var = "cos2",
  select.var = list(contrib = 10),
  repel = TRUE
) +
  labs(
    title = "MetaFAMD - Top 10 Contributing Variables",
    x = "Principal Component 1 (11.14% Variance Explained)",
    y = "Principal Component 2 (9.53% Variance Explained)"
  ) +
  theme(
    axis.text = element_text(size = 16),        # Axis tick labels
    axis.title = element_text(size = 18),       # Axis titles
    plot.title = element_text(size = 20, face = "bold"),  # Plot title
    legend.text = element_text(size = 14),      # Legend text
    legend.title = element_text(size = 16)      # Legend title
  )


# MetaMFA - MUTOS Variable Groups by Contribution
fviz_mfa_var(
  MetaMFA_results,
  choice = "group",
  col.var = "contrib",
  repel = TRUE
  )  +
  labs(
    title = "MetaMFA - MUTOS Variable Groups by Contribution",
    x = "Principal Component 1 (13.52% Variance Explained)",
    y = "Principal Component 2 (10.29% Variance Explained)"
  )



# Hierarchical Cluster Dendrogram on Raw Features (Gower Distance)
gower_dist_matrix <- daisy(McLouth_Covariates, metric = "gower") # Creates distance matrix using Gower's Distance
hclust_result <- hclust(gower_dist_matrix, method = "ward.D2") # Performs hierarchical clustering using the Gower's distance matrix

# Plots the result in hclust_result
fviz_dend(
  hclust_result, k = 9, 
  rect = TRUE, 
  rect_fill = TRUE, 
  show_labels = TRUE, 
  cex = 1.2, 
  main = "Hierarchical Cluster Dendrogram on Raw Features (Gower Distance)"
  ) +
  xlab("Individual Studies Identified by Order in Dataset") +
  ylab("Dissimilarity Distance of Studies")  +
  theme(
    axis.text = element_text(size = 16),        # Axis tick labels
    axis.title = element_text(size = 18),       # Axis titles
    plot.title = element_text(size = 20, face = "bold"),  # Plot title
    legend.text = element_text(size = 14),      # Legend text
    legend.title = element_text(size = 16)      # Legend title
  )



# MetaFAMD + Hierarchical Clustering
HCPC_result <- HCPC(MetaFAMD_results, nb.clust = -1)  # Performs hierarchical clustering on PCs using an automatically chosen nuber of clusters

# Plots the result in HCPC_result
fviz_cluster(HCPC_result, geom = c("point", "text"), main = "MetaFAMD + Hierarchical Clustering")  +
  xlab("Principal Component 1 (11.14% Variance Explained)") +
  ylab("Principal Component 2 (9.53% Variance Explained)") +
  theme(
    axis.text = element_text(size = 16),        # Axis tick labels
    axis.title = element_text(size = 18),       # Axis titles
    plot.title = element_text(size = 20, face = "bold"),  # Plot title
    legend.text = element_text(size = 14),      # Legend text
    legend.title = element_text(size = 16)      # Legend title
  )





