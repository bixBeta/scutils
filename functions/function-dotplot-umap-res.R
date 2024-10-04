for (i in levels(sv4.meta$SCT_snn_res.0.6)) {
  png(filename = paste0("highlight-cells-sv4-res0.6/clust", i, ".png"), width = 1000, height = 600, res = 100)
  highlight.clusters(sobj = sv4, res = 0.6, clust = i)
  dev.off()
}
sv4.meta = sv4@meta.data

lv.feats = c("CD3E","CD247","CD4","CD8A","TRBC1","TRGC1","TRDC","CCR7","CD28","SELL","CD69","FAS","PDCD1","CCL5","IFNG","EOMES","TBX21","CXCR3","IL2","GATA3","IL4R","CCR4","PTGDR2","IRF4","SPI1","KLRB1","CCR6","RORC","AHR","FOXO4","CCR10","CTLA4","FOXP3","IL2RA","LAG3","IL10","IKZF2","NT5E","NRP1","NCAM1","FCGR3A","KLRK1","NCR3","LILRB1","TRPM3","CD14","LGALS2","S100A8","PF4","MS4A1","CD1C","HLA-DRA","ITGAX","CD34")
Idents(sv4) <- sv4$SCT_snn_res.0.6

png("feature_plot-sv4-1.png", width = 1200, height = 1000, res = 100)
FeaturePlot(sv4, features = lv.feats[1:20], raster = F)
dev.off()

png("feature_plot-sv4-2.png", width = 1200, height = 1800, res = 100)
FeaturePlot(sv4, features = lv.feats[21:54], raster = F)
dev.off()




getUmapDotPlots = function(sobj, res, feats){
  
  title=deparse(substitute(sobj))
  Idents(object = sobj) <- sobj@meta.data[,paste0("SCT_snn_res.", res)]
  uplot = DimPlot(sobj , group.by = paste0("SCT_snn_res.", res), label = T, label.size = 6, repel = T, raster = F) + ggtitle(paste0(title, " UMAP -- res ", res))
  
  dplot = DotPlot(object = sobj, assay="RNA", features=feats, col.min=-1.5,
                  group.by = paste0("SCT_snn_res.", res), cols="Spectral", dot.scale=3, cluster.idents = T, scale.by = 'size') + 
    RotatedAxis() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggtitle(paste0(title, " DotPlot -- res ", res))
  
  cplot = uplot + dplot
  
  
  png(filename = paste0("sv4__umap__dotPlot__res__",res, ".png"), width = 3000, height = 900, res = 100)
  print(cplot)
  #Sys.sleep(3)
  dev.off()
  
  # return(cplot)
  
}


getUmapDotPlots(sobj = sv4, res = 0.6, feats = lv.feats)
