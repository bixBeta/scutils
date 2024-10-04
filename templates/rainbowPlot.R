
sv3@meta.data$seurat_clusters <- sv3@meta.data$SCT_snn_res.0.6
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# get number of cells per cluster per sample
n_cells_per_cluster_per_sample3 <- FetchData(sv3, 
                                            vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

write.csv(n_cells_per_cluster_per_sample3, "SEURAT_OBJECT_V3/n_cells_per_cluster_per_sample_SV3.csv")


rownames(n_cells_per_cluster_per_sample3) <- n_cells_per_cluster_per_sample3$orig.ident

n_cells.matrix = n_cells_per_cluster_per_sample3 %>% select(-1)
n_cells.matrix[is.na(n_cells.matrix)] <- 0

n_cells.matrix = n_cells.matrix / rowSums(n_cells.matrix) * 100
n_cells.matrix$sample = n_cells_per_cluster_per_sample3$orig.ident

n_cells.matrix.joined = left_join(n_cells.matrix, sv3.meta, by = c("sample" = "orig.ident")) 

gg3.join = n_cells.matrix.joined %>% select(1:29,batch10x, sample)
gg3.melt = reshape2::melt(gg3.join)

colnames(gg3.melt) <- c("batch10x","sample", "cluster.0.6", "value")

png("SEURAT_OBJECT_V3/distibution_of_cells_per_cluster_per_sample_SV3.png", width = 4000, height = 10000, res = 200)


ggplot(gg3.melt,
       aes(fill=batch10x, y=value, x=sample, label = batch10x)) + 
  geom_bar(position="dodge", stat="identity") + theme(legend.position = "none", 
                                                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_viridis(discrete = T) +
  facet_grid(cluster.0.6 ~ batch10x , scales = "free") + theme(axis.line=element_line())


dev.off()
