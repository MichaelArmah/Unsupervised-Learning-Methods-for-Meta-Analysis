# Code to produce figures used in the example section of the paper ----

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


MetaFAMD_model_reduced <- MetaFAMD(
  data   = McLouth_Data,
  study = "study_id",
  yi    = "g",
  vi    = "Vg",
  mods  = ~ 1,
  rho   = 0.8,
  num_pcs = 11,
  graph = FALSE
)


# Function results
MetaFAMD_results <- MetaFAMD_model$results
MetaFAMD_model_reduced_results <- MetaFAMD_model_reduced$results





# Contribution to the first dimension
p1 <- fviz_contrib(MetaFAMD_results, "var", axes = 2) +
  ggtitle("Contributions of Variables to PC2")
# Contribution to the second dimension
p2 <- fviz_contrib(MetaFAMD_results, "var", axes = 7) +
  ggtitle("Contributions of Variables to PC7")
# Contribution to the ninth dimension
p3 <- fviz_contrib(MetaFAMD_results, "var", axes = 9) +
  ggtitle("Contributions of Variables to PC9")
# Contribution to the tenth dimension
p4 <- fviz_contrib(MetaFAMD_results, "var", axes = 10) +
  ggtitle("Contributions of Variables to PC10")


p1 + p2

p3 + p4


##############################################################
##############################################################
##############################################################

# Old plots for the example





# MetaFAMD - Top 10 Contributing Variables for PC2 adn PC7
p1 <- fviz_famd_var(
  MetaFAMD_results,
  choice = "var",
  col.var = "cos2",
  axes = c(2,7),
  select.var = list(contrib = 10),
  repel = TRUE,
  geom = c("point", "text"),
  labelsize = 8,
  pointsize = 4
) +
  labs(
    title = "MetaFAMD - Top 10 Contributing Variables",
    x = "Principal Component 2 (9.5% Variance Explained)",
    y = "Principal Component 7 (5.9% Variance Explained)",
    color = "cos^2"
  ) +
  theme(
    axis.text = element_text(size = 23),        # Axis tick labels
    axis.title = element_text(size = 23),       # Axis title
    plot.title = element_text(size = 29, face = "bold"),  # Plot title
    legend.text = element_text(size =23),      # Legend text
    legend.title = element_text(size = 23)      # Legend title
  )

# MetaFAMD - Top 10 Contributing Variables for PC9 adn PC10
p2 <- fviz_famd_var(
  MetaFAMD_results,
  choice = "var",
  col.var = "cos2",
  axes = c(9,10),
  select.var = list(contrib = 10),
  repel = TRUE,
  geom = c("point", "text"),
  labelsize = 8,
  pointsize = 4
) +
  labs(
    title = "         ",
    x = "Principal Component 9 (4.5% Variance Explained)",
    y = "Principal Component 10 (3.9% Variance Explained)",
    color = "cos^2"
  ) +
  theme(
    axis.text = element_text(size = 23),        # Axis tick labels
    axis.title = element_text(size = 23),       # Axis titles
    plot.title = element_text(size = 26, face = "bold"),  # Plot title
    legend.text = element_text(size = 23),      # Legend text
    legend.title = element_text(size = 23)      # Legend title
  )

# plots the two biplots side by side using the library "patchwork"
p1 + p2


# Hierarchical clustering on reduced dimension space (Dimension spaced reduced to the first 11 PCs)
Example_HCPC_results <- HCPC(MetaFAMD_model_reduced_results, nb.clust = -1)  # Performs hierarchical clustering on PCs using an automatically chosen nuber of clusters

fviz_cluster(Example_HCPC_results, axes = c(2, 9), geom = c("point", "text"), labelsize = 15,
             pointsize = 2, main = "MetaFAMD + Hierarchical Clustering")  +
  xlab("Principal Component 2 (9.5% Variance Explained)") +
  ylab("Principal Component 9 (4.5% Variance Explained)") +
  theme(
    axis.text = element_text(size = 23),        # Axis tick labels
    axis.title = element_text(size = 23),       # Axis titles
    plot.title = element_text(size = 26, face = "bold"),  # Plot title
    legend.text = element_text(size = 23),      # Legend text
    legend.title = element_text(size = 23)      # Legend title
  )





