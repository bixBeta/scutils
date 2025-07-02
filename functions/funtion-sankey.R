library(dplyr)
library(tibble)
library(networkD3)
library(htmlwidgets)
# Import ----
meta = readRDS("data/adt_harmony_meta.RDS")

# Clean UP ----

meta = meta |> rownames_to_column("bc") |> select(bc, matches("harmony"))
gex_res = 0.8
adt_res = 0.2

getSankey <- function(gex_res, adt_res, floor = 10) {
  prefixg = paste0("GEX_")
  prefixa = paste0("ADT_")
  message(paste0(prefixg, gex_res))
  message(paste0(prefixa, adt_res))

  meta.df = meta |>
    select(
      bc,
      paste0("harmony_clusters_res", gex_res),
      paste0("harmony_decont_clusters_res", adt_res)
    )
  meta.df[[2]] <- as.character(meta.df[[2]])
  meta.df[[3]] <- as.character(meta.df[[3]])

  matrix = as.matrix(table(meta.df[[2]], meta.df[[3]]))

  colnames(matrix) <- paste0(prefixa, "_", colnames(matrix))
  rownames(matrix) <- paste0(prefixg, "_", rownames(matrix))

  df = as.data.frame(matrix)

  colnames(df) <- c("source", "target", "value")

  df = df %>% filter(value > floor)
  df = df |> arrange(desc(value))

  nodes <- data.frame(
    name = c(as.character(df$source), as.character(df$target)) %>% unique()
  )
  df$IDsource = match(df$source, nodes$name) - 1
  df$IDtarget = match(df$target, nodes$name) - 1
  ColourScal = 'd3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'

  p = # Make the Network
    sankeyNetwork(
      Links = df,
      Nodes = nodes,
      Source = "IDsource",
      Target = "IDtarget",
      Value = "value",
      NodeID = "name",
      sinksRight = FALSE,
      colourScale = ColourScal,
      nodeWidth = 50,
      fontSize = 13,
      nodePadding = 20
    )

  return(p)
}


g8a2 = getSankey(gex_res = 0.8, adt_res = 0.2)
g1a3 = getSankey(gex_res = 1, adt_res = 0.3)
g15a4 = getSankey(gex_res = 1.5, adt_res = 0.4)
g2a5 = getSankey(gex_res = 2, adt_res = 0.5)

saveWidget(g8a2, file = "figures/gex_0.8_adt_0.2.html", selfcontained = TRUE)
saveWidget(g1a3, file = "figures/gex_1_adt_0.3.html", selfcontained = TRUE)
saveWidget(g15a4, file = "figures/gex_1.5_adt_0.4.html", selfcontained = TRUE)
saveWidget(g2a5, file = "figures/gex_2_adt_0.5.html", selfcontained = TRUE)
