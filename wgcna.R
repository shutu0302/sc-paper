## -----------------------------------------------------------
library(WGCNA)
library(tidyverse)


gene_exp <- read.csv(file = 'data/wgcna/gene_exp.csv',
                     row.names = 1)


##-----------------------------

## find variable genes
seurat <- FindVariableFeatures(seob, nfeatures =4000)
DefaultAssay(seurat) <- "RNA"
seurat <- NormalizeData(seurat, 
                        normalization.method = "LogNormalize", 
                        scale.factor = 10000) 

seurat <- ScaleData(seurat, vars.to.regress = c("percent.mt", "UMI"))

gene_exp <- as.matrix(seurat@assays$RNA@data)
datTraits <- seurat@meta.data

#colnames(gene_exp)[grepl("[12]_Cel",colnames(gene_exp))]
gene_exp <- gene_exp[intersect(Seurat::VariableFeatures(seurat),rownames(gene_exp)),]

# 
datExpr <- t(gene_exp)

##corrlate two dataframe (same row order)
datExpr <- datExpr[match(rownames(datTraits), rownames(datExpr)), ]

## -----------------------------------------------------------
{
  
  # Gene filter if needed
  # gsg <- goodSamplesGenes(
    #datExpr0, 
    #minFraction = 1/2 
    #)
  #datExpr <- datExpr0[gsg$goodSamples, gsg$goodGenes]
  
  
  # clustering dendro
  plot(hclust(dist(datExpr)), 
       cex.lab = 1.5,
       cex.axis = 1.5, 
       cex.main = 2)
  

#   library(genefilter)
#   var_gene_exp <- varFilter(
#     as.matrix(t(datExpr)),
#     var.func = IQR,
#     var.cutoff = 0.05, # below 0.5
#     filterByQuantile = TRUE)
  # 
  #datExpr <- t(var_gene_exp)
  datExpr[1:4,1:4]
}



## ----message=FALSE------------------------------------------
##find beta
{
  library(WGCNA)
  # multi threads
  enableWGCNAThreads(nThreads = 6)
  ## disableWGCNAThreads()
    
  # define power
  sft <- pickSoftThreshold(
    datExpr, 
    powerVector = 1:20, # 
    networkType = "signed" 
    )
}


## -----------------------------------------------------------
{
  # plot
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(ggthemes)
  
  fig_power1 <- ggplot(data = sft$fitIndices,
         aes(x = Power,
             y = SFT.R.sq)) +
    geom_point(color = 'red') +
    geom_text_repel(aes(label = Power)) +
    geom_hline(aes(yintercept = 0.85), color = 'red') +
    labs(title = 'Scale independence',
         x = 'Soft Threshold (power)',
         y = 'Scale Free Topology Model Fit,signed R^2') +
    theme_few() +
    theme(plot.title = element_text(hjust = 0.5))
     
  fig_power2 <- ggplot(data = sft$fitIndices,
         aes(x = Power,
             y = mean.k.)) +
    geom_point(color = 'red') +
    geom_text_repel(aes(label = Power)) +
    labs(title = 'Mean connectivity',
         x = 'Soft Threshold (power)',
         y = 'Mean Connectivity') +
    theme_few()+
    theme(plot.title = element_text(hjust = 0.5))
    
  plot_grid(fig_power1, fig_power2)
}


## ----message=FALSE------------------------------------------
##network
{
  net <- blockwiseModules(
    # 0.input data
    datExpr, 
    
    # 1. correlation
    corType = "pearson", # 相关系数算法，pearson|bicor
      
    # 2. matrix
    power = 6, # 前面得到的 soft power
    networkType = "signed", # unsigned | signed | signed hybrid
      
    # 3. TOM 
    TOMType = "unsigned", # none | unsigned | signed
    saveTOMs = TRUE,
    saveTOMFileBase = "blockwiseTOM",
  
    # 4. 
    deepSplit = 6, # 0|1|2|3|4, 值越大得到的模块就越多越小
    minModuleSize = 20,
    
    # 5. 
    mergeCutHeight = 0.25, 
  
    # others
    numericLabels = FALSE, 
    nThreads = 0, 
    maxBlockSize = 100000 
    )
  # check gene numbers
  table(net$colors) 
}


## -----------------------------------------------------------
{
  library(tidyverse)
  wgcna_result <- data.frame(gene_id = names(net$colors),
                 module = net$colors) #%>%
    #left_join(gene_info, by = c('gene_id' = 'gene_id')) 
  head(wgcna_result)
}


## -----------------------------------------------------------
{
  # plot
  plotDendroAndColors(
    dendro = net$dendrograms[[1]], 
    colors = net$colors,
    groupLabels = "Module colors",
    dendroLabels = FALSE, 
    addGuide = TRUE)
}


## ----message=FALSE------------------------------------------
{
  # 
  moduleTraitCor <- cor(
    net$MEs,
    datTraits,
    use = "p",
    method = 'spearman' # 注意相关系数计算方式
    )
  
  #  Pvalue
  moduleTraitPvalue <- corPvalueStudent(
    moduleTraitCor, 
    nrow(datExpr))
}


## -----------------------------------------------------------
{
  #  heatmap
  sizeGrWindow(10,6)
  
  #  pvalue
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) <- dim(moduleTraitCor)
  
  
  # heatmap 
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(net$MEs),
                 ySymbols = names(net$MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(500),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
}


## -----------------------------------------------------------
{
  MET <- orderMEs(cbind(net$MEs, dplyr::select(datTraits,CellType)))
  
  plotEigengeneNetworks(
    multiME = MET, 
    setLabels = "Eigengene dendrogram", 
    plotDendrograms = TRUE, 
    plotHeatmaps = FALSE,
    colorLabels = TRUE,
    marHeatmap = c(3,4,2,2))
}

#sample dendrogram
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",
                    dendroLabels = F)
#color legend
colVars <- list(CellType=c("1"= "#FFFFFF", 
                  "2"="#FFCDC1", 
                        "3"="#FF977D", 
                         "4"="#FF653E", 
                         "5"="#FF3300"))

colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), datTraits$CellType)]
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))


## -----------------------------------------------------------
{
  plotEigengeneNetworks(
    multiME = MET, 
    setLabels = "Eigengene dendrogram", 
    plotDendrograms = FALSE, 
    plotHeatmaps = TRUE,
    colorLabels = TRUE,
    marHeatmap = c(8,8,2,2))
}


##cytoscape
moduleLabels = net$colors
moduleColors <- dplyr::pull(wgcna_result, module)

# 。
ADJ1=abs(cor(dataExpr,use="p"))^softPower 
# kWithin，kOut=kTotal-kWithin， kDiff=kIn-kOut=2*kIN-kTotal 
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors) 

module = "turquoise";
# Select module probes
probes = colnames(dataExpr) 
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
modProbes

##LOADING TOM
load(net$TOMFiles, verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM

##NETWORK
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.2,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
)
## -----------------------------------------------------------
{
  # select module
  my_modules <- c('turquoise') 
  
  # retrive matrxi
  m_wgcna_result <- filter(wgcna_result, module %in% my_modules)
  m_datExpr <- datExpr[, m_wgcna_result$gene_id]
  
  # TOM 
  m_TOM <- TOMsimilarityFromExpr(
    m_datExpr,
    power = 3,
    networkType = "unsigned",
    TOMType = "unsigned")
    
  dimnames(m_TOM) <- list(colnames(m_datExpr), colnames(m_datExpr))
    
  #  Cytoscape output
  cyt <- exportNetworkToCytoscape(
    m_TOM,
    edgeFile = "CytoscapeInput-network2.txt",
    weighted = TRUE,
    threshold = 0.0005) #0.2
}

##################other graphs (plotMat)
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
#(1) another plot for realtionship between module eigengenes
plotMEpairs(MEs,y=datTraits$CellType)
#(2) Diagnostics: heatmap plots of module expression
sizeGrWindow(8,9)
#par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
# for black module
which.module="turquoise";
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
#（3） Diagnostics: displaying module heatmap and the eigengene
sizeGrWindow(8,7);
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="MPP")

##########################
module = which.module
geneTree = hclust(as.dist(dissTOM), method = "average")
plotMat(t(scale(datExpr[,moduleColors==module ])),nrgcols=30,rlabels=F,rcols=module,main=module, cex.main=2)

#Correlate the module eigengenes with the trait
y=datTraits$CellType

#plot Gene signicance (y-axis) vs. intramodular connectivity (x-axis) 
colorlevels=unique(moduleColors)
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels))) 
{
  whichmodule=colorlevels[[i]];
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(k.in$kWithin[restrict1],
                     GeneSignificance[restrict1], col=moduleColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}

#Generalizing intramodular connectivity for all genes on the array
datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
# Display the first few rows of the data frame
head(datKME)



#Output of the results of network screening analysis
NS1=networkScreening(y=y, datME=MEs, datExpr=datExpr,oddPower=3, blockSize=10000, minimumSampleSize=4,addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)
GeneResultsNetworkScreening=data.frame(GeneName=row.names(NS1), NS1)
#write.table(GeneResultsNetworkScreening, file="GeneResultsNetworkScreening.csv",row.names=F,sep=",")
MEsy = data.frame(y, MEs)
eigengeneSignificance = cor(MEsy, y);
eigengeneSignificance[1,1] = (1+max(eigengeneSignificance[-1, 1]))/2
eigengeneSignificance.pvalue = corPvalueStudent(eigengeneSignificance, nSamples = length(y))
namesME=names(MEsy)
# Form a summary data frame
out1=data.frame(t(data.frame(eigengeneSignificance,eigengeneSignificance.pvalue, namesME, t(MEsy))))
# Set appropriate row names
dimnames(out1)[[1]][1]="EigengeneSignificance"
dimnames(out1)[[1]][2]="EigengeneSignificancePvalue"
dimnames(out1)[[1]][3]="ModuleEigengeneName"
dimnames(out1)[[1]][-c(1:3)]=dimnames(datExpr)[[1]]
# Write the data frame into a file
write.table(out1, file="MEResultsNetworkScreening.csv", row.names=TRUE, col.names = TRUE, sep=",")
GeneName=dimnames(datExpr)[[2]]
GeneSummary=data.frame(GeneName, moduleColors, NS1)
write.table(GeneSummary, file="GeneSummaryTutorial.csv", row.names=F,sep=",")
datTraits=data.frame(ArrayName, MEsy)
dimnames(datTraits)[[2]][2:length(namesME)]=paste("Trait",dimnames(datTraits)[[2]][2:length(namesME)],sep=".")
write.table(datTraits, file="TraitsTutorial.csv", row.names=F,sep=",")

#Relationships among the top 30 most signicant genes,correlation heatmaps for signed network:
names(datExpr)=colnames(datExpr)
#??????????????????????
sizeGrWindow(7,7)
NS1=networkScreening(y=y, datME=MEs, datExpr=datExpr,oddPower=3, blockSize=100000, minimumSampleSize=4,addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)
topList=rank(NS1$p.Weighted,ties.method="first")<=100
gene.names= names(datExpr)[topList]
# The following shows the correlations between the top genes
plotNetworkHeatmap(datExpr, plotGenes = gene.names,networkType="unsigned", useTOM=FALSE,power=1, main="signed correlations")



GS_spe <-  as.numeric(apply(datExpr[[1]]$data,2,function(x){
  biserial.cor(x,datTraits[,trait_minP$trait[i]])
}))
GeneSignificance_spe <- abs(GS_spe)
#
NS1=networkScreening(y=allTraits[,trait_minP$trait[[i]]], 
                     datME=datME, 
                     datExpr=multiExpr[[1]]$data, 
                     oddPower=3, 
                     blockSize=1000, 
                     minimumSampleSize=4, 
                     addMEy=TRUE, 
                     removeDiag=FALSE, 
                     weightESy=0.5) 
rownames(NS1) <- colnames(multiExpr[[1]]$data)

# biserial.cor, select genes
FilterGenes_spe = ((GeneSignificance_spe > 0.2) & (abs(datKME[paste("MM.",Freq_MS_max$GS_color,sep="")])>0.8) & (NS1$q.Weighted < 0.01) ) 
table(FilterGenes_spe)
# 
trait_hubGenes_spe <- colnames(multiExpr[[1]]$data)[FilterGenes_spe] 

## --------------find hub gene-----------------
ADJ1=abs(cor(datExpr,use="p"))^4 
# 
# 
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors) 
head(Alldegrees1)

datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
datKME[1:4,1:4]

table(moduleColors)
module = "red"
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
modProbes

## ---------------plot module significance-------------
GS1=as.numeric(cor(y,datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, moduleColors, mean, na.rm=T)

sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,moduleColors)

