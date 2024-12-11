# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# CellProps ----
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------

getCounts = function(sobjmeta_, res_) {
  
  meta = sobjmeta_
  #meta = meta |> filter(scDblFinder.class == "singlet")
  
  ncells.table = dplyr::count(meta, get(noquote(paste0("harmony_clusters_res", res_))), orig.ident)
  
  colnames(ncells.table)[1] <- "harmony_cluster"
  
  
  max_cluster = length(levels(meta[[paste0("harmony_clusters_res", res_)]]))
  ncells.table$factor = factor(ncells.table$harmony_cluster, levels = 0:max_cluster)
  
  total_cells_per_cluster = ncells.table |>
    group_by(orig.ident) |>
    summarise(Sum = sum(n)) |>
    ungroup()
  
  
  ncells.table2 = left_join(ncells.table, total_cells_per_cluster, by = "orig.ident")
  
  ncells.table2$norm = (ncells.table2$n / ncells.table2$Sum) * 100
  
  ns1 = ggplot(ncells.table2, aes(x=harmony_cluster, y=norm, fill=orig.ident)) +
    geom_bar(position="stack", stat="identity") + 
    scale_color_manual(values = colors, aesthetics = c("colour", "fill")) +
    theme_linedraw() + ggtitle(paste0("Ncells Per Sample Per Cluster -- Res = "), res_) + ylab("norm.by.sample")
  
  ns2 = ggplot(ncells.table2, aes(x=factor, y=norm, fill=orig.ident)) +
    geom_bar(position="fill", stat="identity") +  
    scale_color_manual(values = colors, aesthetics = c("colour", "fill")) + 
    theme_linedraw() + ggtitle(paste0("Percentage Contribution of each Sample Per Cluster -- Res = ", res_)) + xlab("harmony_cluster") +
    ylab("norm.by.sample")
  
  
  #print(ns2)
  
  return(list(ncells.table, ncells.table2, plot = ns2))
  
  
  
}


saveNcells = function(list_, rez_){
  
  t.list = lapply(list_, function(x){
    getCounts(sobjmeta_ = x, res_ = rez_)}
  )
  
  for (i in 1:length(t.list)) {
    
    png(filename = paste0("figures/ncells_", names(t.list)[i], "_res_", rez_ ,  ".png"), width = 2200, height = 800, res = 150)
    
    print(t.list |> pluck(i, "plot"))
    
    dev.off()
    
    
  }
  
  
  
  
}
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# PseudoBulks ----
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------

getByCluster = function(bulk_, cluster_){
  
  a = bulk_ |> select(matches(paste0("_", cluster_, "$"))) 
  return(a)
  
  
}

getPseudoBulk = function(sobj_, res_){
  
  bulk <- AggregateExpression(sobj_, group.by = c("orig.ident", paste0("harmony_clusters_res", res_)), return.seurat = F)
  bulk.matrix = as.data.frame(bulk$RNA) 
  
  nclusters = length(unlist(map(strsplit(colnames(bulk.matrix), split = "_"), ~pluck(.,2))) |> unique())
  
  pb.list = list()
  
  for (i in 1:nclusters) {
    message(paste("Getting Cluster:", i-1))
    pb.list[[i]] <-   getByCluster(bulk_ = bulk.matrix, cluster_ = i-1)
    names(pb.list)[[i]] <- paste0("cluster__", i-1, "__res", res_)
    
  }
  
  
  return(pb.list)
  
  
}



# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# DDS DESEQ2 ---- 
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------


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


getPCAsbatch= function(vst_){
  
  
  meta = as.data.frame(colData(vst_))
  rv <- rowVars(assay(vst_))
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  
  pca <- prcomp(t(assay(vst_)[select,]))
  
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  pVar.df <- as.data.frame(percentVar)
  pVar.df$x = as.factor(paste0("PC",rownames(pVar.df)))
  
  pVar.df = pVar.df[ , order(names(pVar.df))]
  pVar.df$percentVar = pVar.df$percentVar * 100
  pVar.df$percentVar = round(pVar.df$percentVar, digits = 2)
  
  
  d <- data.frame(pca$x, label=rownames(pca$x))
  d2 <- left_join(d, meta, by = c("label" = "id_10x"))
  
  
  
  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(ggplot2))
  
  pc1 = ggplot(d2, aes(x=PC1, y=PC2, color = batch)) +
    geom_point(size=5, aes(shape=group)) +
    geom_label_repel(aes(label = label),
                     box.padding   = 0.8, 
                     point.padding = 0.5,
                     segment.color = 'grey55', show.legend = F) + 
    xlab(paste0(pVar.df$x[1], "  ", pVar.df$percentVar[1], "%") ) +
    ylab(paste0(pVar.df$x[2], "  ", pVar.df$percentVar[2], "%") ) + 
    theme_light() + theme(legend.position = "bottom")  + 
    scale_color_manual(values = colors2) 
    # geom_mark_hull(aes(fill = batch),
    #                concavity = 20, alpha = 0.08
    #)
  
  # png(paste0(pin, "_PC1_PC2.png"), width = 1200, height = 1200, res = 150)
  # pc1
  # dev.off()
  
  
  
  
  return(list(
    prcomp.out = pca,
    Variance.df    = pVar.df,
    colData    = meta,
    PCA.df      = d2,
    plot = pc1
  ))
  
}








# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# GLOBAL ---- 
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------

colors <- (c(paletteer_d("ggsci::springfield_simpsons"), 
             paletteer_d("ggsci::hallmarks_light_cosmic"),
             paletteer_d("tidyquant::tq_light"),
             paletteer_d("ggsci::dark_uchicago"),
             paletteer_d("ggsci::default_nejm"),
             paletteer_d("ggsci::lanonc_lancet"),
             paletteer_d("ggthemes::Tableau_20")
             
)
)

colors2 <-  c("#EF8A62",
            "#1f78b4",
            "#1b9e77",
            "purple3",
            "khaki4",
            "#E9A3C9",
            "#A1D76A",
            "red",
            "grey")

