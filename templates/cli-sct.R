#> input sobj.by.batch.list.unfiltered.RDS
#> old universal filtering criteria:
#>          sobj.filtered <- subset(sobj.raw, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20 & log10GenesPerUMI > 0.80 & nFeature_ADT > 30)

#> 


.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.2")
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(reticulate)
library(glmGamPoi)
library(dplyr)
library(plotly)



sobj.by.batch <- readRDS("/workdir/TCELL_VDJ_ADT_ALL_h5_August13_2024/NEW/RDS/02_sobj_by_batch_list_unfilltered.RDS")
mdata.raw = read.csv("new-mdata.raw.csv", header = T, row.names = 1)

processSOBJ <- function(sobj, filter = TRUE, nfeatMin = 500, nfeatMax = 6000, mt = 20, noveltyScore = 0.80, nfeatADTMin = 30, prefix = "SOBJ"){
  
  if(filter){
    
    pdf(file = paste0("figures/", prefix, "_PRE_FILTER_RAW_VLN_PLOT.pdf"), width = 14, height = 18)
    a = VlnPlot(sobj, features = colnames(mdata.raw)[c(2:7)],
                pt.size = 0, group.by = "enid", ncol = 3, raster = F, split.by = "group")
    b = VlnPlot(sobj, features = colnames(mdata.raw)[c(2:7)],
                pt.size = 0, group.by = "orig.ident", ncol = 3, raster = F)
    print(a /b) 
    dev.off()

    png(filename =  paste0("figures/", prefix, "_PRE_FILTER_RAW_VLN_PLOT.byENID.png"), width = 1920, height = 1080, res= 150)
    c = VlnPlot(sobj, features = colnames(mdata.raw)[c(2:7)],
                pt.size = 0, group.by = "enid", ncol = 3, raster = F, split.by = "group")
    print(c)
    dev.off()
    
    png(filename =  paste0("figures/", prefix, "_PRE_FILTER_RAW_VLN_PLOT.bySample.png"), width = 1920, height = 1080, res= 150)
    d = VlnPlot(sobj, features = colnames(mdata.raw)[c(2:7)],
                pt.size = 0, group.by = "orig.ident", ncol = 3, raster = F)
    
    print(d) 
    dev.off()
    
    
    
    
    sobj.filtered <- subset(sobj, subset = nFeature_RNA > nfeatMin & 
                              nFeature_RNA < nfeatMax & 
                                percent.mt < mt & 
                                  log10GenesPerUMI > noveltyScore & 
                                    nFeature_ADT > nfeatADTMin)
    
    
    pdf(file = paste0("figures/", prefix, "_POST_FILTER_RAW_VLN_PLOT.pdf"), width = 14, height = 18)
    a = VlnPlot(sobj.filtered, features = colnames(mdata.raw)[c(2:7)],
                pt.size = 0, group.by = "enid", ncol = 3, raster = F, split.by = "group")
    b = VlnPlot(sobj.filtered, features = colnames(mdata.raw)[c(2:7)],
                pt.size = 0, group.by = "orig.ident", ncol = 3, raster = F)
    print(a /b) 
    dev.off()
    
    
    
    png(filename =  paste0("figures/", prefix, "_POST_FILTER_RAW_VLN_PLOT.byENID.png"), width = 1920, height = 1080, res= 150)
    c = VlnPlot(sobj.filtered, features = colnames(mdata.raw)[c(2:7)],
                pt.size = 0, group.by = "enid", ncol = 3, raster = F, split.by = "group")
    print(c)
    dev.off()
    
    
    png(filename =  paste0("figures/", prefix, "_POST_FILTER_RAW_VLN_PLOT.bySample.png"), width = 1920, height = 1080, res= 150)
    
    d = VlnPlot(sobj.filtered, features = colnames(mdata.raw)[c(2:7)],
                pt.size = 0, group.by = "orig.ident", ncol = 3, raster = F)
    print(d) 
    dev.off()
    
    return(sobj.filtered)
    
    } 
}

sobj.by.batch.filtered = list()
for (i in 1:length(sobj.by.batch)) {
  sobj.by.batch.filtered[[i]] <- processSOBJ(sobj = sobj.by.batch[[i]], filter = TRUE, prefix = paste0("sobj_", names(sobj.by.batch[i])))
  names(sobj.by.batch.filtered)[[i]] <- names(sobj.by.batch)[i]
}

saveRDS(sobj.by.batch.filtered, "RDS/03_sobj_by_batch_list_filltered.RDS")


# split layers by orig.ident ----

for (i in 1:length(sobj.by.batch.filtered)) {
  
  sobj.by.batch.filtered[[i]][["RNA"]] <- split(sobj.by.batch.filtered[[i]][["RNA"]], f = sobj.by.batch.filtered[[i]]$orig.ident)
}


# run sctransform ---- 
sobj.by.batch.filtered <- lapply(sobj.by.batch.filtered, function(x){
  
  x <- SCTransform(object = x , assay = "RNA", vars.to.regress = "percent.mt", vst.flavor = "v2", return.only.var.genes = T, verbose = T)
  
})

saveRDS(sobj.by.batch.filtered, file = "RDS/04_SCTransformed_sobjs_by_batch_list_unIntegrated.RDS")
 

for (i in 1:length(sobj.by.batch.filtered)) {
  saveRDS(object = sobj.by.batch.filtered[[i]], file = paste0("RDS/sobj_scTransformed_unintegrated_", names(sobj.by.batch.filtered)[i], ".RDS" ))
}




