# scRNAseq of mouse cochleae
Analysis pipeline for the data reported in "Systematic single cell RNA sequencing analysis reveals unique transcriptional regulatory networks of Atoh1-mediated hair cell conversion in adult mouse cochleae"

# Contents
- Scanpy analysis (Figure 1)\
\
  Contains code used to reanalyse data from GSE173217. Used Jupyter notebook.
  
- SCENIC and WGCNA (Figure 2)\
\
  Contains R scripts used to analyse Gene regulatory network (GRN) inference using WGCNA and SCENIC. It also has script for GO analysis of WGCNA modules in the supplementary figure using clusterprofiler. The ARACNE results used to find the GRNs of each modules are uploaded in the preprocessed data folder. Requires the R packages WGCNA (V1.70-3), SCENIC (v1.1.2-2) and clusterprofiler (v4.0.2).
  
- TI method and SCT integration (Figure 3)\
\
  Contains scripts used to analyse PAGA trajectory, seurat analysis of retriving and combining GSE137299 and GSE173217, and the SCTransform integration of two datasets. In addition to that, extracted sct data used for the correlation plot in Figure 3E was added to the preprocessed data folder below. Requires R package Seurat (v4.0.1) and dyno (v0.1.2).
  
- Viper (Figure 4)\
\
  Scripts used for viper analysis in Figure 4. It also applied for the viper analysis used in plotting the supplementary figure. Requires R package 1.22.0.
  
  
# Downloads
Preprocessed data and other files necessary for this analysis can be downloaded from https://drive.google.com/drive/folders/1RYiemRLYgBm2D8WaYobz3GTPgZGnVo8C?usp=share_link
