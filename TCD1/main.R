library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(DESeq2)
library(magrittr)
library(stringr)
library(purrr)
library(aplot)
library(paletteer)

# Import ----
# --------------------------------------------------------------------------------------------------------------
setwd("/workdir/TCELL_PRED1_DoubletFilteredSCTHarmony/gg03-session")
sobj <- readRDS("/workdir/TCELL_PRED1_DoubletFilteredSCTHarmony/work/cb/d907ff995a5a44f6886c4a8fbdab45/preD1_dblFiltered_Harmony-by-sample-sobj.RDS")

sobj.meta = sobj@meta.data


# CellProps ----
# --------------------------------------------------------------------------------------------------------------
mdata.list = list(all_batches = sobj.meta)

saveNcells(list_ = mdata.list, rez_ = 0.8)
saveNcells(list_ = mdata.list, rez_ = 1)
saveNcells(list_ = mdata.list, rez_ = 1.5)
saveNcells(list_ = mdata.list, rez_ = 2)



# PseudoBulks ----
# --------------------------------------------------------------------------------------------------------------
res1.5_pb = getPseudoBulk(sobj_ = sobj, res_ = 1.5)
res1_pb = getPseudoBulk(sobj_ = sobj, res_ = 1)

saveRDS(res1.5_pb, file = "RDS/preD1_dblFiltered_res1.5_pb.RDS")
saveRDS(res1_pb, file = "RDS/preD1_dblFiltered_res1_pb.RDS")



# DDS ----
# --------------------------------------------------------------------------------------------------------------

metadata = sobj@meta.data

box.meta = metadata |> select(c(orig.ident, sample_name, batch, group, enid, day))
box.meta = unique(box.meta)
rownames(box.meta) <- box.meta$orig.ident
metadata = box.meta 
metadata$id_10x = metadata$orig.ident

saveRDS(box.meta, "RDS/box.meta.RDS")


all.pbs = readRDS("RDS/preD1_dblFiltered_res1_pb.RDS")
# filter to keep only first 20 clusters
all.pbs2 = all.pbs[1:20]


all.dds <- lapply(all.pbs2, function(x){
  setDDS(l_ = x)
})


all.dds.deseq = lapply(all.dds, function(x){
  
  DESeq(object = x, minReplicatesForReplace = Inf)
  
})

saveRDS(all.dds.deseq, "RDS/preD1_dblFiltered_Top20_clusters_postDESeq_Call_dds_res.1.RDS")


all.dds.vst = lapply(all.dds.deseq, function(x){
  
  varianceStabilizingTransformation(object = x)
})

plotPCA(all.dds.vst$cluster__1__res1, intgroup = "group", returnData = F) 

saveRDS(all.dds.vst, "RDS/preD1_dblFiltered_Top20_clusters_VST_res.1.RDS")


pca.plots.list = lapply(X = all.dds.vst, FUN = function(x){getPCAsbatch(vst_ = x)})
pca.only = lapply(pca.plots.list, function(x){
  pluck(x,"plot")
})

png(filename = "figures/preD1_dblFiltered_PCA_PLOTS_res1_batch.png", res = 150, height = 4000, width = 5000)
plot_list(gglist =pca.only, ncol = 5, byrow = T, widths = 800, heights = 400)
dev.off()

# UCELL ----
# --------------------------------------------------------------------------------------------------------------


library(UCell)


markers <- list()
markers$Tcell_CD4 <- c("CD4", "CD40LG")
markers$Tcell_CD8 <- c("CD8A", "CD8B")
markers$Tcell_Treg <- c("FOXP3", "IL2RA")
markers$Tcell_MAIT <- c("KLRB1", "SLC4A10", "NCR3")
markers$Tcell_gd <- c("TRDC", "TRGC1", "TRGC2", "TRDV1", "TRAC-", "TRBC1-", "TRBC2-")
markers$Tcell_NK <- c("FGFBP2", "SPON2", "KLRF1", "FCGR3A", "KLRD1", "TRDC", "CD3E-",
                      "CD3G-")


Idents(sobj) <- sobj$harmony_clusters_res1.5

sobj <- AddModuleScore_UCell(sobj, features = markers, assay = "SCT")

sobj.meta = sobj@meta.data

saveRDS(sobj.meta, "RDS/preD1_dblFiltered_UCELL-meta.RDS")

signature.names <- paste0(names(markers), "_UCell")


png("figures/preD1_dblFiltered_ucell.png", width = 2400, height = 1600, res = 200)

VlnPlot(sobj, features = signature.names, group.by = "harmony_clusters_res1.5", stack = T, flip = T, cols = turbo(8), sort = "increasing")

dev.off()





# SCPLOTTER ----
# --------------------------------------------------------------------------------------------------------------


library(scplotter)

png("preD1_dblFiltered_celldim_stat-d1.png", width = 2400, height = 1600, res = 150)

CellDimPlot(sobj,
            group_by = "harmony_clusters_res1.5", reduction = "umap.harmony", stat_by = "group",
            stat_plot_type = "ring", stat_plot_label = TRUE, stat_plot_size = 0.05, label_repel = T, raster = F)

dev.off()


png("preD1_dblFiltered_cell_stat_sankey2-d1.png", width = 8400, height = 8600, res = 150)

CellStatPlot(sobj, plot_type = "sankey", alpha = .6,
             group_by ="group")

dev.off()



png("preD1_dblFiltered_cell_stat_circos-d1.png", width = 1200, height = 1200, res = 100)

CellStatPlot(sobj, plot_type = "circos", alpha = .6,
             group_by = "group", ident = "harmony_clusters_res1.5")

dev.off()


png("preD1_dblFiltered_Clustree-d1.png", width = 1200, height = 1200, res = 100)
ClustreePlot(sobj, prefix = "harmony_clusters_res")
dev.off()


FeatureStatPlot(sobj, features = markers$Tcell_CD4,
            ident = "harmony_clusters_res1.5", plot_type = "heatmap")

p1 = FeatureStatPlot(sobj, markers$Tcell_CD4, reduction = "umap.harmony",
                lower_cutoff = 1, upper_cutoff = 4, plot_type = "dim")

p2 = FeatureStatPlot(sobj, setdiff(m2,x1), reduction = "umap.harmony",
                     lower_cutoff = 1, upper_cutoff = 4, plot_type = "dim")

png("preD1_dblFiltered_TCELLmarkersUMAP.png", width = 2800, height = 2800, res = 200)
FeatureStatPlot(sobj, setdiff(m2,x1), reduction = "umap.harmony",
                lower_cutoff = 1, upper_cutoff = 4, plot_type = "dim", add_density = TRUE, density_color = "gray22")
dev.off()

png("preD1_dblFiltered_TCELLmarkers.png", width = 2800, height = 3200, res = 150)
FeatureStatPlot(sobj, features = setdiff(m2,x1),
                ident = "harmony_clusters_res1.5", 
                add_bg = TRUE, stack = TRUE,
                legend.position = "bottom", legend.direction = "horizontal")

dev.off()

