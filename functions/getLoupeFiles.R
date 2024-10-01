x = read.delim("list.log", header =T)

processSOBJ <- function(sobj, batch_){
  
  Idents(sobj) <- sobj$harmony_clusters_res1
  assay <- sobj[["SCT"]]
  counts <- counts_matrix_from_assay(assay)
  
  create_loupe(
    counts,
    clusters = select_clusters(sobj),
    projections = select_projections(sobj),
    output_name = paste0("loupe/sobj_batch_", batch_, "_res1_SCT")
  )
  
}

for (i in 1:1) {
  
  batch_ = (strsplit(x[i,1], split = "_")[[1]][2])
  assign(x = paste0("sobj.batch.", batch_), value = readRDS(x[i,1]), envir = .GlobalEnv)
  
  processSOBJ(sobj = get(paste0("sobj.batch.", batch_)), batch_ = batch_)
  
  rm(list = c(paste0("sobj.batch.", batch_)))
  
}



