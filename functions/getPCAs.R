library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(DESeq2)
library(magrittr)
library(stringr)
library(purrr)
library(aplot)

getPCAs= function(vst_){
  
  
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
 
   pc1 = ggplot(d2, aes(x=PC1, y=PC2, color = group)) +
    geom_point(size=8, aes(shape=group)) +
    geom_label_repel(aes(label = label),
                     box.padding   = 0.8, 
                     point.padding = 0.5,
                     segment.color = 'grey55', show.legend = F) + 
    xlab(paste0(pVar.df$x[1], "  ", pVar.df$percentVar[1], "%") ) +
    ylab(paste0(pVar.df$x[2], "  ", pVar.df$percentVar[2], "%") ) + 
     theme_light() + theme(legend.position = "bottom")  + 
     scale_color_manual(values = c("darkgreen", "purple"))
  
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
getPCAs(vst_ = all.dds.vst$batch_10$cluster__0__res1)


pca.plots.list = list()
  
for (i in 1:length(batches)) {

  pca.plots.list[i] <- lapply(all.dds.vst[[batches[i]]], function(x){

    lapply(X = all.dds.vst[[batches[i]]], FUN = function(x){getPCAs(vst_ = x)})
    

  })

  names(pca.plots.list)[i] <- names(all.dds.vst)[i]


}
saveRDS(pca.plots.list, "PCAs_res1.RDS")

pca.plots.list = readRDS("PCAs_res1.RDS")

pca.plots.list |> pluck(1, 1,'plot')


# map(mylist, ~pluck(.,1,1) %>% map(pluck, "value_I_want")) 

pca.only = lapply(pca.plots.list, function(x){
  
  map(x, ~ pluck(.,"plot")) 
  
})

png(filename = "test.png", res = 100, height = 4000, width = 4000)
plot_list(gglist =pca.only$batch_10, ncol = 5, byrow = T, widths = 800, heights = 400)
dev.off()


for (i in 1:length(batches)) {
  png(filename = paste0("PCAS_RES1_New//PCA_PLOTS_res1_", batches[i], ".png"), width = 4000, height = 3000, res = 110)
  print(plot_list(gglist = pca.only[[batches[i]]], ncol = 5, byrow = T, widths = 400, heights = 400))
  dev.off()
}
