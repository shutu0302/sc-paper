#---------------------clusterprofile------------------
library(clusterProfiler)
rvcheck::check_bioc("clusterProfiler")
library(tidyverse)

{
  tib_wgcna_result <- as.tibble(wgcns_result)
  yellow_gene <- filter(tib_wgcna_result, module == "yellow") %>%
    select(gene_id) %>%
    pull(gene_id)
  yellow_gene <- as.character(yellow_gene)
  
  new_ids_yellow <- bitr(yellow_gene, fromType = 'SYMBOL', 
                         toType = c('SYMBOL', 'ENTREZID'), 
                         OrgDb = 'org.Hs.eg.db')
  yellow_ID <- dplyr::select(new_ids_yellow,ENTREZID) %>%
    pull(ENTREZID)
}

## run KEGG
kk <- enrichKEGG(gene         = red_ID,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
head(kk)

## KEGG module
mkk <- enrichMKEGG(gene = red_ID,
                   organism = 'mmu')

head(mkk)

## ------------GO ANALYSIS----------

## ENRICH
yellow_go <- enrichGO(gene         = new_ids_yellow$ENTREZID,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = 'ENTREZID',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
head(yellow_go)

yellow_go <- setReadable(yellow_go, OrgDb = org.Mm.eg.db)

yellow_go_simplify <- clusterProfiler::simplify(yellow_go, cutoff=0.7, by="p.adjust")

## visulization
library(enrichplot)
upsetplot(turquoise_go)

cnetplot(turquoise_go, node_label="category") 

heatplot(yellow_go_simplify, showCategory = 15) + ggplot2::ylab('Annotations of Biological Processes')+
  ggplot2::xlab('Gene Symbols')+
  ggplot2::ggtitle('yellow module')