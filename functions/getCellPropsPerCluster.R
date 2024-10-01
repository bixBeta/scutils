dplyr::count(sobjmeta, harmony_clusters_res1, orig.ident)


getCounts = function(sobjmeta_, res_) {
  
  meta = sobjmeta_
  meta = meta |> filter(scDblFinder.class == "singlet")

  ncells.table = dplyr::count(meta, get(noquote(paste0("harmony_clusters_res", res_))), orig.ident)
  
  colnames(ncells.table)[1] <- "harmony_cluster"
  
  
  max_cluster = max(ncells.table$harmony_cluster)
  ncells.table$factor = factor(ncells.table$harmony_cluster, levels = 0:max_cluster)
  
  total_cells_per_cluster = ncells.table |>
                                group_by(orig.ident) |>
                                summarise(Sum = sum(n)) |>
                                ungroup()
  
  
  ncells.table2 = left_join(ncells.table, total_cells_per_cluster, by = "orig.ident")
  
  ncells.table2$norm = (ncells.table2$n / ncells.table2$Sum) * 100
  
  ns1 = ggplot(ncells.table2, aes(x=harmony_cluster, y=norm, fill=orig.ident)) +
    geom_bar(position="stack", stat="identity") + 
    scale_color_manual(values = colors2, aesthetics = c("colour", "fill")) +
    theme_linedraw() + ggtitle(paste0("Ncells Per Sample Per Cluster -- Res = "), res_) + ylab("norm.by.sample")
  
  ns2 = ggplot(ncells.table2, aes(x=factor, y=norm, fill=orig.ident)) +
    geom_bar(position="fill", stat="identity") +  
    scale_color_manual(values = colors2, aesthetics = c("colour", "fill")) + 
    theme_linedraw() + ggtitle(paste0("Percentage of Cells Per Sample Per Cluster -- Res = ", res_)) + xlab("harmony_cluster") +
    ylab("norm.by.sample")
  
  
  #print(ns2)
  
  return(list(ncells.table, ncells.table2, plot = ns2))
  
  
  
}

mdata.paths = read.csv("metadata.paths.csv", header = F, col.names = "meta.data.path")

mdata.list = list()
for (i in 1:nrow(mdata.paths)) {
  mdata.list[[i]] = read.csv(mdata.paths$meta.data.path[i], row.names = 1)
  names(mdata.list)[[i]] <- gsub(pattern = "_scDBL_doublet_identities.csv", replacement = "", x = basename(mdata.paths$meta.data.path[i]))

}


saveNcells = function(list_, rez_){
  
  t.list = lapply(list_, function(x){
    getCounts(sobjmeta_ = x, res_ = rez_)}
  )
  
  for (i in 1:length(t.list)) {
    
    png(filename = paste0("ncells_figures/ncells_", names(t.list)[i], "_res_", rez_ ,  ".png"), width = 1920, height = 600, res = 150)
    
    print(t.list |> pluck(i, "plot"))
    
    dev.off()
    
    
  }
  
  
  
  
}


saveNcells(list_ = mdata.list, rez_ = 0.8)
saveNcells(list_ = mdata.list, rez_ = 1)
saveNcells(list_ = mdata.list, rez_ = 1.5)
saveNcells(list_ = mdata.list, rez_ = 2)



getDblPercents <- function(sobjmeta_){
  
  meta = sobjmeta_
  
  x = table(meta$scDblFinder.class)
  x["doublet.percentage"] <- round((x["doublet"] / (sum(x))) * 100, digits = 2) 
  
  x.df = stack(x)
  
  colnames(x.df) <- c("scDblCounts", "scDblClass")
  
  return(x.df)
  
}

dblPercents.list = lapply(mdata.list, function(x){
  getDblPercents(sobjmeta_ = x)
})

saveRDS(dblPercents.list, "doublet_percentages.RDS")


sink(file = "dblPercents.log")
lapply(dblPercents.list, function(x){
  x
})
sink()


doubletFiltered.sobj.meta.list <- lapply(mdata.list, function(x){
  
  fmeta = x |> filter(scDblFinder.class == "singlet")
  return(fmeta)
})

saveRDS(doubletFiltered.sobj.meta.list, "doublet_filtered_meta_sobjs.RDS")
