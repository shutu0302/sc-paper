## GSE137299  scRNAseq  - reanalysis

#####------------------------------------------

## data from GSE137299 was retrived and then processed using Cellranger

library(Seurat, lib.loc = "/pub/software/R_Lib")
library(dplyr, lib.loc = "/pub/software/R_Lib")
library(patchwork, lib.loc = "/pub/software/R_Lib")


seob[["percent.mt"]] <- PercentageFeatureSet(seob, pattern = "^mt-")

# Visualize QC metrics as a violin plot
Idents(seob) <- "Time"
VlnPlot(seob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seob <- subset(seob, 
               subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 40 & nCount_RNA < 20000)

seob <- NormalizeData(seob, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)
seob <- FindVariableFeatures(seob, selection.method = "vst", nfeatures = 3000)


all.genes <- rownames(seob)
seob <- ScaleData(seob, features = all.genes)

seob <- RunPCA(seob, features = VariableFeatures(object = seob))

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
seob <- JackStraw(seob, num.replicate = 100 ,dims = 70)
seob <- ScoreJackStraw(seob, dims = 1:50)

JackStrawPlot(seob, dims = 1:50)

seob <- FindNeighbors(seob, dims = 1:50)
seob <- FindClusters(seob, resolution = 0.8)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
seob <- RunUMAP(seob, dims = 1:50)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(seob, reduction = "umap")

seob <- RunTSNE(seob, dims = 1:50)
DimPlot(seob, reduction = "tsne", label = T) + NoLegend()
DimPlot(seob, reduction = "tsne", label = T, split.by = "Time") + NoLegend()

## Expression plot

FeaturePlot(seob, features = c("Myo7a", "Pou4f3"), reduction = "tsne")

FeaturePlot(seob, features = c("Slc26a5", "Ocm"), reduction = "tsne")
FeaturePlot(seob, features = c("Slc17a8", "Otof"), reduction = "tsne")

# HC markers

FeaturePlot(seob, 
            features = c("Myo6", "Myo7a", "Rasd2", "Chrna9", "Pvalb", "Pou4f3", "Chrna10"), 
            reduction = "tsne")


# --------

##iHC_E14 L.PsC_E14  iOHC_E16 L.PsC_E16 DC/OPC_P1    IPC_P1    OHC_P1    OHC_P7 
##45       123       214       540       284       151       490       105 
##DC/PC_P7 
##187 


# subset object contains cHCs1-3 and OHCs from Yamashita et al

Idents(subset) <- factor(Idents(subset), 
                         levels = c("iHC_E14", "iOHC_E16", "OHC_P1","OHC_P7","OHCs",
                                    "cHCs1","cHCs2","cHCs3"))
DoHeatmap(subset, features = top30$gene
) + NoLegend()

seob_SCT <- merge(seob, y = subset)

# NatCom PlosGen 
# 2433     204 


DefaultAssay(seob_SCT) <- "RNA"
Idents(seob_SCT) <- "ident"

seob_list <- SplitObject(seob_SCT, split.by = "ident")

for(i in 1:length(seob_list)){
  seob_list[[i]] <- SCTransform(
    seob_list[[i]], 
    variable.features.n = 3000,
    # vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
    verbose = FALSE)
}


# 2. Integration
## select genes
features <- SelectIntegrationFeatures(object.list = seob_list,
                                      nfeatures = 3000)
## integrate
seob_list <- PrepSCTIntegration(object.list = seob_list, 
                                anchor.features = features)
## find anchors，15-30 min
anchors <- FindIntegrationAnchors(object.list = seob_list, 
                                  # reference = 3 can speed up
                                  normalization.method = "SCT",
                                  anchor.features = features
)
## 
seob_SCT <- IntegrateData(anchorset = anchors, 
                       normalization.method = "SCT")






