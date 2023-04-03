#---------------- SCENIC --------------------

#for a Seurat object
# select cells
cells.use <- WhichCells(object = seob, ident = 
                          c("cHCs1", "cHCs2", "cHCs3", "OHCs","IHCs", "SCs_P12", "SCs_P33")) 
# obtain matrix
expr <- GetAssayData(object = cells, assay= "RNA", slot = "data")
# format transform 

expr <- as(Class = 'matrix', object = expr)
exprMat <- expr
cellInfo <- data.frame(seuratCluster=Idents(cells),
                       nCount = cells@meta.data[["nCount_RNA"]],
                       nGene = cells@meta.data[["nFeature_RNA"]])


cellInfo <- data.frame(cellInfo)
cellTypeColumn <- "seuratCluster"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"

head(cellInfo)

#create folder
dir.create("SCENIC analysis")
setwd("SCENIC analysis")
#
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

#assign colors

colVars <- list(CellType=c("cHCs1"="forestgreen", 
                           "cHCs2"="darkorange", 
                           "cHCs3"="magenta4", 
                           "OHCs"="hotpink", 
                           "IHCs"="red3", 
                           "SCs_P12"="skyblue",
                           "SCs_P33"="darkblue"
))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")

plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

## ---------- SCENIC settings -------------

library(SCENIC)

mgi_dbs <- list('500bp'= 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather', 
                '10kb' = 'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
db_mcVersion <- 'v9'

scenicOptions <- initializeScenic(org='mgi',
                                  dbs = mgi_dbs,
                                  dbDir="cisTarget_databases",
                                  datasetTitle='Scenic Analysis LA_LAL', 
                                  nCores=10) 

scenicOptions@settings$db_mcVersion <- db_mcVersion
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#----------------------

#Filter by the total number of reads per gene
#Filter by the number of cells in which the gene is detected
#only the genes that are available in RcisTarget databases will be kept

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.005*ncol(exprMat),
                           minSamples=ncol(exprMat)*.005)

exprMat_filtered <- exprMat[genesKept, ]

interestingGenes <- c("Sox2", "Atoh1-HA", "Myo7a","Kdm1a")
interestingGenes[which(!interestingGenes %in% genesKept)]

## [1] 770 200

# delete
rm(exprMat)

# ---------
runCorrelation(exprMat, scenicOptions)


# To be run
runGenie3(exprMat_filtered, scenicOptions)

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) #** Only for toy run!!
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat)

saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status