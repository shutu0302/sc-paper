#---------------------- dyno ---------------

library(dyno)
library(tidyverse)
library(SummarizedExperiment)

DefaultAssay(cells2) <- "RNA"

# expr matrix after normalization

gene_exp <- cells2@assays[["RNA"]]@data
gene_counts <- cells2@assays[["RNA"]]@counts

# 3. cell info

cell_info <- rownames_to_column(cells2@meta.data, 
                                var = "cell_id")

###  wrap_expression 

dynob_cells2 <- wrap_expression(
  expression = t(gene_exp), 
  counts = t(gene_counts), 
  cell_info = cell_info
)

#dynob <- add_prior_information(
#dynob,
# start_id = c("Sox2+","DC/PCs")
#)

dynob_cells2 <- add_prior_information(
  dataset = dynob_cells2,
  start_id = colnames(cells2)[cells2$cell_type == "Sox2+"]
)

dynob_cells2 <- add_grouping(
  dynob_cells2,
  cell_info$cell_type
)
head(dynob$grouping)

## TI method

guidelines_shiny(dataset = dynob_cells2)


# Reproduces the guidelines as created in the shiny app
answers <- dynguidelines::answer_questions(
  multiple_disconnected = FALSE, 
  expect_topology = TRUE, 
  expected_topology = "bifurcation", 
  n_cells = 1618, 
  n_features = 18629, 
  memory = "2GB", 
  prior_information = "start_id", 
  docker = TRUE
)
guidelines <- dynguidelines::guidelines(answers = answers)

methods_selected <- guidelines$methods_selected
methods_selected

model <- infer_trajectory(dynob_cells2, 
                          method = "paga_tree", 
                          give_priors = c("start_id"),docker = TRUE)

plot_dimred(model, 
            grouping = dynob_cells2$grouping) + 
  ggtitle("Cell grouping")

plot_dimred(model, 
            grouping = dynob_cells2[["cell_info"]][["sample"]]) + 
  ggtitle("Group by sampe")

plot_dimred(model, 
            feature_oi = "Slc17a8", 
            expression_source = dynob_cells2) + 
  ggtitle("VGluT3")                                  


plot_dimred(model, "pseudotime", pseudotime = calculate_pseudotime(model)) 
+ ggtitle("Pseudotime")

