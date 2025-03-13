library(decontX)
library(Seurat)
library(ggplot2)
library(patchwork)
library(reticulate)
library(dplyr)
library(scplotter)
library(dplyr)
library(tibble)
library(tidyr)
library(networkD3)

sobj.meta2 = read.csv("/SSD/copy-from-gg03/preD1_dblFiltered_metadata.csv", header = T, row.names = 1)
s35 = readRDS("cli_rds/high2_delta_high_bg_sd_sample35.RDS")
s35.sobj = s35$sobj_decontX
s35.meta = s35.sobj@meta.data



#ref_set = sobj.meta2 |> filter(orig.ident == "24A035" & harmony_clusters_res1 %in% c(2,22,23,29,35))
ref_set = sobj.meta2 |> filter(orig.ident == "24A035")

ref_set$bcf = rownames(ref_set)
ref_set = separate(ref_set, col = "bcf", into = c("bc", "metaN"), sep = "_")

query_set = s35.meta |> filter(decontXcounts_snn_res.0.5 %in% c(13,17,11,10,9,16,2,8))
query_set$bc = rownames(query_set)

str(ref_set$harmony_clusters_res1)
str(query_set$decontXcounts_snn_res.0.5)

query_set$decontXcounts_snn_res.0.5 = as.integer(as.numeric(as.character(query_set$decontXcounts_snn_res.0.5)))
str(query_set$decontXcounts_snn_res.0.5)


# ---- 
getSankey <- function(ref_set, query_set, cluster_col_ref, bc_col_ref, cluster_col_query, bc_col_query){
  
  df1 <- ref_set |> select(!!sym(bc_col_ref), !!sym(cluster_col_ref))
  df2 <- query_set |> select(!!sym(bc_col_query), !!sym(cluster_col_query))
  
  # Rename columns for clarity
  colnames(df1) <- c("barcode", "cluster_ref")
  df1$cluster_ref <- paste0("reference_", df1$cluster_ref)
  
  colnames(df2) <- c("barcode", "cluster_query")
  df2$cluster_query <- paste0("query_", df2$cluster_query)
  
  # Merge dataframes based on barcode
  merged_df <- merge(df1, df2, by = "barcode")
  merged_df$cluster_ref = as.character(merged_df$cluster_ref)
  merged_df$cluster_query = as.character(merged_df$cluster_query)
  
  # Prepare nodes for Sankey plot
  nodes <- unique(c(merged_df$cluster_ref, merged_df$cluster_query))
  nodes_df <- data.frame(name = nodes)
  
  # Prepare links for Sankey plot
  links <- merged_df %>%
    mutate(source = match(cluster_ref, nodes) - 1,
           target = match(cluster_query, nodes) - 1,
           value = 1) %>%
    group_by(source, target) %>%
    summarise(value = sum(value)) %>%
    ungroup()
  
  # Create Sankey plot
  sankey_plot <- sankeyNetwork(Links = links, Nodes = nodes_df,
                               Source = "source", Target = "target",
                               Value = "value", NodeID = "name",
                               fontSize = 24, nodeWidth = 60, nodePadding=20, sinksRight=FALSE, height = 1080,
                               width = 1920, iterations = 0)
  
  return(list(sankey = sankey_plot,
              links = links, 
              nodes_df = nodes_df,
              merged_df = merged_df))
  
}

s35q_sobj_ref <- getSankey(ref_set = ref_set , 
                          query_set = query_set,
                          cluster_col_ref = "harmony_clusters_res1",
                          bc_col_ref = "bc", 
                          cluster_col_query = "decontXcounts_snn_res.0.5",
                          bc_col_query = "bc" )

saveRDS(s35q_sobj_ref, "s35q_sobj_ref.RDS")
