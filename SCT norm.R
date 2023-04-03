seob3[['SCT']] <- NULL

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
seob3 <- IntegrateData(anchorset = anchors, 
                      normalization.method = "SCT")
DefaultAssay(seob3) <- "integrated"





