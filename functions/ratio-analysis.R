# load("/Users/faraz/Desktop/4052-4061/rawCounts/4052D-cpet/4052D-cpet.RData")
# save.image(file = "data/4052.dds_SARtools.RData")

library(ggplot2)
library(DESeq2)
library(dplyr)
library(tibble)
library(scales)

load("data/4052.dds_SARtools.RData")
dds = out.DESeq2$dds


extract.filter.normCounts = function(dds){
  
  # extract coldata
  meta.data = colData(dds) %>% data.frame()
  
  # add enid.id column for grepping later
  
  # meta.data$orig.ident = paste0("X", meta.data$label)
  # meta.data$orig.ident = gsub(pattern = "-", replacement = ".", x = meta.data$orig.ident)
  colnames(meta.data)[which(colnames(meta.data)== 'label')] <- 'orig.ident'
  meta.data$enid.id = paste0(meta.data$orig.ident, "_", meta.data$enid)
  # extract normalized counts
  norm.counts = counts(dds, normalized = T) %>% data.frame()
  
  # add row medians to norm.counts
  norm.counts$row.medians = apply(X = norm.counts, MARGIN = 1, FUN = median)
  
  # get the quantile profile of the median
  q.tile = quantile(norm.counts$row.medians)
  
  # bins.quantiles(norm.counts$row.medians, target.bins = 4, max.breaks = 15)$binct
  q.tile[4] <- 5
  q75.filtered.norm.counts = norm.counts[norm.counts$row.medians > q.tile[4], ]
  
  
  
  # plot norm counts pre and post filter
  
  counts.df.no.filter = norm.counts %>% select(-row.medians)
  counts.df.stacked.no.filter = stack(counts.df.no.filter)
  gg.df.no.filter = left_join(counts.df.stacked.no.filter, meta.data, by = c("ind" = "orig.ident"))
  
  p.pre.filter <- ggplot(gg.df.no.filter, aes(x=.data$values+1)) +
    stat_density(aes(group=.data$ind, color=.data$day), position="identity", geom="line", show.legend=TRUE) +
    scale_x_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(~10^.x))) +
    labs(color="") +
    xlab(paste0("normCounts_", deparse(substitute(dds)))) +
    ylab("Density") +
    ggtitle("Density of counts distribution") +
    theme_gray() + facet_wrap("enid")
  
  
  counts.df.post.filter = q75.filtered.norm.counts %>% select(-row.medians)
  counts.df.stacked.post.filter = stack(counts.df.post.filter)
  gg.df.post.filter = left_join(counts.df.stacked.post.filter, meta.data, by = c("ind" = "orig.ident"))
  
  p.post.filter <- ggplot(gg.df.post.filter, aes(x=.data$values+1)) +
    stat_density(aes(group=.data$ind, color=.data$day), position="identity", geom="line", show.legend=TRUE) +
    scale_x_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(~10^.x))) +
    labs(color="") +
    xlab(paste0("normCounts_", deparse(substitute(dds)), "__", "Row_median >", round(q.tile[4], 2))) +
    ylab("Density") +
    ggtitle("Density of counts distribution") +
    theme_gray() + facet_wrap("ind")
  
  
  
  #return(p.pre.filter)
  return(list(colData = meta.data,
              norm.counts.No.filter = norm.counts,
              Quantile = q.tile,
              q75.filtered.norm.counts  = q75.filtered.norm.counts,
              ggplots = 
                list(ggplot.pre.filter = p.pre.filter,
                     ggplot.post.filter = p.post.filter)
              
  )
  )
  
}

xl = extract.filter.normCounts(dds = out.DESeq2$dds)
saveRDS(xl, file = "data/xl__quantile.bin.filtered.normalizedCounts.RDS")

# original beta function
getRatios = function(cluster){
  
  metadata = cluster[["colData"]] %>% data.frame()
  metadata$enid.id =  paste0(metadata$orig.ident,"_",metadata$ENID, "_",metadata$day)
  
  filtered.counts = cluster[["q75.filtered.norm.counts"]] %>% data.frame() %>% select(-row.medians)
  
  # return(filtered.counts %>% head()) 
  colnames(filtered.counts) <- metadata$enid.id
  
  meta.filter.day1 = metadata %>% filter(day == "D1") %>% select(orig.ident)
  meta.filter.day2 = metadata %>% filter(day == "D2") %>% select(orig.ident)
  
  d1 = unname(unlist(meta.filter.day1))
  d2 = unname(unlist(meta.filter.day2))
  
  counts.day1 = filtered.counts %>% rownames_to_column("gene") %>% select(matches("D1"), "gene") %>% column_to_rownames("gene")
  counts.day2 = filtered.counts %>% rownames_to_column("gene") %>% select(matches("D2"), "gene") %>% column_to_rownames("gene")
  
  counts.day1 = counts.day1 + 0.1
  counts.day2 = counts.day2 + 0.1
  
  enids = unique(metadata$ENID)
  ratios = list()
  
  for (i in 1:length(enids)) {
    
    ratios[[i]] <- as.data.frame(counts.day2[,grep(enids[i], colnames(counts.day2))]  / counts.day1[,grep(enids[i], colnames(counts.day1))])
    
    names(ratios)[[i]] <- enids[i]
    
    rownames(ratios[[i]]) <- rownames(counts.day2)
    
    colnames(ratios[[i]]) <- enids[i]
    
  }
  
  ratios.matrix =  as.matrix(do.call(cbind, ratios))
  log10.ratios.matrix = log10(ratios.matrix)
  
  return(list(D1 = counts.day1,
              D2 = counts.day2, 
              raw.ratios.list = ratios,
              non.log.ratios = ratios.matrix,
              log.10.ratios = log10.ratios.matrix,
              colData = metadata))
}
