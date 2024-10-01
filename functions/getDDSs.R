#> input = pseudobulk list at a given res; current available res1 and res1.5
#> 
#> 

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(DESeq2)
library(magrittr)
library(stringr)
library(purrr)
library(aplot)

# mdata.paths = read.csv(file = "metadata.paths.csv", header = T)
# mdata.list = list()
# 
# for (i in 1:nrow(mdata.paths)) {
#   mdata.list[[i]] <- read.csv(mdata.paths$mata.data[1])
#   names(mdata.list)[[i]] <- strsplit(basename(mdata.paths[i,1]), split = "_metadata.csv")[[1]][1]
# }

mdata.list = readRDS("doublet_filtered_meta_sobjs.RDS")
metadata = read.csv("/workdir/TCELL_VDJ_ADT_ALL_h5_August13_2024/NEW/new-metadata.csv", header = T)

metadata = metadata |> filter(day == "D1pre")

setDDS = function(l_){
  
  count.matrix = l_
  count.matrix = magrittr::set_colnames(count.matrix, value = str_split_i(string = colnames(count.matrix), pattern = "_",i = 1))
  
  target = metadata |> filter(id_10x %in% colnames(count.matrix)) |> as.data.frame()
  
  rownames(target) <- target$id_10x
  
  count.matrix = count.matrix |> select(all_of(target$id_10x))
  
  dds = DESeqDataSetFromMatrix(countData = count.matrix,
                               colData = target, design =  ~ group)

  return(dds)

}

all.pbs = readRDS("Pseudobulk_RES.1.5.RDS")

# filter to keep only first 20 clusters
all.pbs2 = lapply(all.pbs, FUN = function(x){
  x[1:20]
})

all.dds = list()
batches = names(all.pbs2)

for (i in 1:length(batches)) {
  all.dds[[i]] <- lapply(all.pbs2[[batches[i]]], function(x){
    setDDS(l_ = x)
  })
  names(all.dds)[[i]] <- names(all.pbs2)[i]
}

all.dds.deseq = list()
for (i in 1:length(batches)) {

  all.dds.deseq[[i]] <- lapply(all.dds[[batches[i]]], function(x){
    
    DESeq(object = x, minReplicatesForReplace = Inf)
  })
    
  names(all.dds.deseq)[[i]] <- names(all.dds)[i]
}

saveRDS(all.dds.deseq, "Top20_clusters_postDESeq_Call_dds_res.1.5.RDS")

all.dds.vst = list()

for (i in 1:length(batches)) {
  
  all.dds.vst[[i]] <- lapply(all.dds.deseq[[batches[i]]], function(x){
    
    varianceStabilizingTransformation(object = x)
  })
  
  names(all.dds.vst)[[i]] <- names(all.dds.deseq)[i]
}

saveRDS(all.dds.vst, "Top20_clusters_VST_res.1.5.RDS")


plotPCA(all.dds.vst$batch_1$cluster__0__res1, intgroup = "group", returnData = T) 



pca.plots.list = list()

# for (i in 1:length(batches)) {
#  
#   pca.plots.list[[i]] <- lapply(all.dds.vst[[batches[i]]], function(x){
#     
#     plotPCA(x, intgroup = "group") + theme_light() + theme(legend.position = "bottom")  + 
#       scale_color_manual(values = c("green2", "purple"))  + geom_point(size=7 , aes(shape=group))
#     
#   })
#   
#   names(pca.plots.list)[[i]] <- names(all.dds.vst)[i]
#    
#   
# }

# library(aplot)
# png(filename = "PCA_PLOTS_res1.5_all_batches.png", width = 4000, height = 4000, res = 150)
# plot_list(gglist = pca.plots.list$batch_1, ncol = 5, byrow = T, widths = 200, heights = 100)
# dev.off()


# for (i in 1:length(batches)) {
#   png(filename = paste0("pca_res1/PCA_PLOTS_res1_", batches[i], ".png"), width = 4000, height = 4000, res = 150)
#   print(plot_list(gglist = pca.plots.list[[batches[i]]], ncol = 5, byrow = T, widths = 200, heights = 100))
#   dev.off()
# }


