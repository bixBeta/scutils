library(tibble)
library(dplyr)
library(networkD3)
library(tidyr)
library(reshape2)
library(viridis)
sv4.raw.t.meta = sv4.raw.t@meta.data
v4t.obj.added.clusts.meta <- readRDS("/home/MAIN2/February_03_CITEseq_fixed_AB_annot/SEURAT_OBJECT_V4_T_Added_Clusters/v4t.obj.added.clusts.meta.RDS")


v4raw.v.v4t = as.matrix(table(sv4.raw.t.meta$SCT_snn_res.0.6, v4t.obj.added.clusts.meta$SCT_snn_res.0.8))
colnames(v4raw.v.v4t) = paste0("res0.8_T_cluster", colnames(v4raw.v.v4t))
rownames(v4raw.v.v4t) = paste0("res0.6_SV4_Cluster", rownames(v4raw.v.v4t))

getSankey <- function(matrix, floor){
  
  x = as.data.frame(matrix) 
  colnames(x) = c("source", "target", "value")
  
  x = x %>% filter(value > floor)
  
  nodes = as.data.frame(c(as.character(x$source), as.character(x$target)) %>% unique())
  
  colnames(nodes) = "name"
  
  x$IDsource=match(x$source, nodes$name)-1 
  x$IDtarget=match(x$target, nodes$name)-1
  
  y = sankeyNetwork(Links = x, Nodes = nodes,
                    Source = "IDsource", Target = "IDtarget",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=12, nodePadding=20, iterations = 0)
  
  return(y)
  
}
getSankey(matrix = v4raw.v.v4t, floor = 10)
